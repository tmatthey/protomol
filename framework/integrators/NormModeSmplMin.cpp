#include "NormModeSmplMin.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ReducedHessAngle.h"
#include "reducedHessBond.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#if defined(HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
extern "C"{
        void SimTK_about_SimTKlapack(const char *, int, char *);
	    void SimTK_version_SimTKlapack(int *, int *, int *);
	}
#endif

using namespace std;


using namespace ProtoMol::Report;

using std::string;
using std::vector;


#define FIXCOFM	1

namespace ProtoMol {
  //__________________________________________________ NormModeSmplMin

  const string NormModeSmplMin::keyword( "NormModeSmplMin" );

  NormModeSmplMin::NormModeSmplMin() : STSIntegrator() 
		{}

  NormModeSmplMin::NormModeSmplMin(Real timestep, Real minimlim, bool masswt, bool multts, bool rforce, bool rediag,
					 ForceGroup *overloadedForces) 
    : STSIntegrator(timestep,overloadedForces), minLim(minimlim), massWeight(masswt), mts(multts), randforce(rforce), reDiag(rediag)
		{}


  NormModeSmplMin::~NormModeSmplMin() 
  {  
  }

  void NormModeSmplMin::initialize(GenericTopology *topo,
				      Vector3DBlock   *positions,
				      Vector3DBlock   *velocities,
				      ScalarStructure *energies){
    //cout << "VALUE CHECK MIN: " << mhQ[0] << endl;
    STSIntegrator::initialize(topo, positions, velocities, energies);
    initializeForces();
	myPrevMTS  = (NormModeInt*)myPreviousIntegrator;
	//calc sizes/initialize
    _N = (int)myPositions->size(); //get N
	pk.resize(_N);
	randStp = 0.0; //fix random step for first loop
	//***********************************************************
  }

  //Steepest decent minimizer for all modes outside subspace.
  //PK update respects 'c' subspace positions.
  int NormModeSmplMin::updtConstrainedModes(Real peLim) {
	int in, itr;
	Real oldPot, lambda, lambda1, lambdaSlp1, lambdaRef, lastDiff; 
	int rsCG;

	rsCG = 0;
	lastDiff = 5; //start at 5 kcal mol^{-1}
	totGamma = 0.0;
	numGam = 0;
	//
	//find new forces
	itr = 0;
	calculateForces();
	forceCalc++;
	//Set first value of /lambda to be 1/eigval 
	lambdaRef = 1.0 / myPrevMTS->maxEigval;	
	if(!massWeight) lambdaRef /= 2.0;
	lambda = lambdaRef;
	if(mts){
		Real ts = getTimestep();
		report <<debug(5)<<"[NormModeSmplMin::updtConstrainedModes] MTS steps= "<<ts<<endl;
		for(in=0;in<(int)ts;in++){
			itr++;
			if(in>0){
				calculateForces();
				forceCalc++;
			}
			//find forces in compliment space.
			myPrevMTS->nonSubspaceForce(myForces, myForces);
			//sift so position move is in compliment space
			for(int i=0;i<_N;i++) (*myForces)[i] /= myTopo->atoms[i].scaledMass;
			(*myPositions).intoWeightedAdd(lambda,*myForces);
		}
	}else{
		for(in=0;in<100;in++){
			itr++;
			//
			report.precision(10);
			report <<debug(6)<<"[NormModeSmplMin::updtConstrainedModes] PE= "<<myEnergies->potentialEnergy()<<endr;
			//****find search direction vector pk
			//find forces in compliment space.
			myPrevMTS->nonSubspaceForce(myForces, myForces);
			//sift so position move is in compliment space
			if(massWeight){
				for(int i=0;i<_N;i++) (*myForces)[i] /= myTopo->atoms[i].scaledMass;
			}else{ 
				myPrevMTS->nonSubspacePosition(myForces, myForces);
			}
			//set pk
			pk = *myForces;
			//find slope of original PE with /lambda=0 here. If +ve then no solution so reset CG
			lambdaSlp1 = 0.0;
			for( int k = 0; k < _N; k++ ) lambdaSlp1 -= pk[k].dot((*myForces)[k]);
			//save PE at /lambda=0
			oldPot = myEnergies->potentialEnergy();
			report <<debug(7)<<"[NormModeSmplMin::updtConstrainedModes] lambd= "<<lambda<<endl;
			//find force at new position pos+/lambda*pk
			(*myPositions).intoWeightedAdd(lambda,pk);
			//update total gamma
			totGamma += lambda;
			numGam++;
			calculateForces();
			forceCalc++;
			//test for end, too large lambda test first
			if((oldPot - myEnergies->potentialEnergy()) < 0){
				if(rsCG>4) report << error << "[NormModeSmplMin::updtConstrainedModes] Minimization failed. Rediagonalization maybe required for the current conformation. Aborting."<<endr;
				else{
					if(!reDiag){  //allow minimization if mode at angle to sub-space
						//calc optimum lambda from first slope 
						Real a1;
						a1 = (myEnergies->potentialEnergy() - oldPot - lambdaSlp1 * lambda) / (lambda * lambda);
						lambda1 = -lambdaSlp1 / (2 * a1);
						//Test that the quadratic solution gives a predicted PE value
						//where the difference from the old PE value is bounded by the difference
						//from the last succesful step, else solve quadratic for the last difference.  
						Real calcPE = a1*lambda1*lambda1+lambdaSlp1*lambda1+oldPot;
						if(oldPot - calcPE > lastDiff && lambdaSlp1 != 0.0){
							lambda1 = (-lastDiff * 2.0) / lambdaSlp1;
						}
						(*myPositions).intoWeightedAdd(-lambda,pk);		//reset positions
						totGamma -= lambda;
						numGam--;
						calculateForces();
						forceCalc++;
						if(lambda1 > 0.0 && lambda1 < lambda) lambda = lambda1;
						else lambda /= 2.0;
						rsCG++;
						report <<debug(1)<<"[NormModeSmplMin::updtConstrainedModes] Reset CG, PE fail. Cycle= "<<rsCG<<" lambda= "<<lambda<<endl;
					}else{
						(*myPositions).intoWeightedAdd(-lambda,pk);		//reset positions
						myPrevMTS->mhQ = NULL;	//force rediagonalization
						return -1;				//flag aborted
					}
				}
			}else{
				rsCG = 0;
				lambda = lambdaRef; //0.5 / myPrevMTS->maxEigval;	//Factor of 2 to be conservative
			}
			if((oldPot - myEnergies->potentialEnergy()) < peLim && !rsCG) break;
			if(!rsCG) lastDiff = oldPot - myEnergies->potentialEnergy();
			//
		}
	}
	return itr;
  }

  void NormModeSmplMin::run(int numTimesteps) {
	if( numTimesteps != -1 && numTimesteps != -2 ){
		preStepModify();
		postStepModify();
		return;
	}else{
		preStepModify();
		//restore previous pre-purturbed positions
		if(randforce) (*myPositions).intoWeightedAdd(-randStp,pk);
		//do minimization with local forces
		forceCalc = 0;
		myPrevMTS->itrs = updtConstrainedModes(minLim);
		myPrevMTS->minForceCalc = forceCalc;
		if(randforce && myPrevMTS->itrs != -1){	//add random force, but not if rediagonalizing
			//add random force
			if(totGamma > 0 && numGam > 0){
				randStp = sqrt(totGamma / (Real)numGam);
				//randStp = sqrt((getTimestep() * totGamma * Constant::TIMEFACTOR * Constant::TIMEFACTOR)  / (Real)numGam) * Constant::INV_TIMEFACTOR * Constant::INV_TIMEFACTOR;
				pk = myPrevMTS->gaussRandCoordM;	//save random purtubation
				(*myPositions).intoWeightedAdd(randStp,myPrevMTS->gaussRandCoordM);
			}
		}else{
			randStp = 0.0;
		}
		//
		postStepModify();
	}
  }  

  void NormModeSmplMin::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
	parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));
	parameters.push_back(Parameter("massweight",Value(massWeight,ConstraintValueType::NoConstraints()),true,Text("Position move based on mass weighted forces.")));
    parameters.push_back(Parameter("mts",Value(mts,ConstraintValueType::NoConstraints()),false,Text("Multiple time step mode")));
    parameters.push_back(Parameter("randforce",Value(randforce,ConstraintValueType::NoConstraints()),true,Text("Add random force")));
    parameters.push_back(Parameter("rediag",Value(reDiag,ConstraintValueType::NoConstraints()),false,Text("Force re-digonalize")));
  }

  STSIntegrator* NormModeSmplMin::doMake(string& , const vector<Value>& values,ForceGroup* fg)const{
    return new NormModeSmplMin(values[0],values[1],values[2],values[3],values[4],values[5],fg);
  }

}

