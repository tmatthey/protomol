#include "NormModeMin.h"
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


namespace ProtoMol
{
	//__________________________________________________ NormModeMin

	const string NormModeMin::keyword("NormModeMin");

	NormModeMin::NormModeMin() : STSIntegrator()
	{
	}

	NormModeMin::NormModeMin(Real timestep, Real minimlim, bool fpec, bool masswt, bool rforce, bool rediag,
	                         ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces), minLim(minimlim), forcePEcheck(fpec), massWeight(masswt), randforce(rforce), reDiag(rediag)
	{
	}


	NormModeMin::~NormModeMin()
	{
	}

	void NormModeMin::initialize(GenericTopology* topo,
	                             Vector3DBlock* positions,
	                             Vector3DBlock* velocities,
	                             ScalarStructure* energies)
	{
		//cout << "VALUE CHECK MIN: " << mhQ[0] << endl;
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
		myPrevMTS = (NormModeInt*)myPreviousIntegrator;
		//calc sizes/initialize
		_N = (int)myPositions->size(); //get N
		pk.resize(_N);
		gk.resize(_N);
		oldGk.resize(_N);
		randStp = 0.0; //fix random step for first loop
		//***********************************************************
	}

	//Conjugate gradient minimizer for all modes outside subspace.
	//PK update respects 'c' subspace positions.
	int NormModeMin::updtConstrainedModes(Real peLim)
	{
		int in, itr;
		Real oldPot, lambda, lambda1, betaD, betaN, lambdaSlp, lambdaSlp1, newPE, lambdaRef, curLambdaRef;
		bool checkedPE, lambNeg;
		int rsCG;

		rsCG = 0;
		checkedPE = lambNeg = false;
		lambda = 0;
		totGamma = 0.0;
		numGam = 0;
		//Set first value of /lambda to be 1/eigval times the max value of the eigvector
		lambdaRef = 0.5 / myPrevMTS->maxEigval; //Factor of 2 for non-mass weight
		if (!forcePEcheck && !massWeight) lambdaRef *= 2.0; //compatibility with the old method
		if (massWeight) lambdaRef *= 2.0;
		curLambdaRef = lambdaRef;
		//
		int _3N = 3 * _N;
		for (int i = 0; i < _3N; i++) pk[i / 3][i % 3] = 0.0;
		//find new forces
		itr = 0;
		for (in = 0; in < 100; in++)
		{
			itr++;
			//
			if (!checkedPE && !rsCG)
			{
				calculateForces();
				forceCalc++;
			}
			//
			report.precision(10);
			report << debug(6) << "[NormModeMin::updtConstrainedModes] PE= " << myEnergies->potentialEnergy() << endr;
			//****find search direction vector pk
			betaD = betaN = 0.0;
			oldGk = gk;
			//find forces in compliment space.
			myPrevMTS->nonSubspaceForce(myForces, myForces);
			//sift so position move is in compliment space
			if (massWeight && !rsCG)
			{
				for (int i = 0; i < _N; i++) (*myForces)[i] /= myTopo->atoms[i].scaledMass;
			}
			else
			{
				myPrevMTS->nonSubspacePosition(myForces, myForces);
			}
			//set gk t0 -force
			gk.intoWeighted(-1.0, *myForces);
			//PR algorithm for finding /beta, length along search vector so orthogonal
			for (int k = 0; k < _N; k++)
			{
				betaD += oldGk[k].dot(oldGk[k]);//OR use HS: pk[k].dot(gk[k] - oldGk[k]);//
				betaN += gk[k].dot(gk[k] - oldGk[k]);
			}
			if (betaD != 0.0) betaN = betaN / betaD;
			else betaN = 0.0;
			//PK requres reset /beta if negative
			if (betaN < 0) betaN = 0;
			//steepest decent is mass weighted
			if (massWeight) betaN = 0.0;
			//set position change vector
			pk.intoWeighted(betaN, pk);
			pk.intoSubtract(gk);
			//
			//****line search along proposed vector pk
			//save PE at /lambda=0
			oldPot = myEnergies->potentialEnergy();
			//find slope of original PE with /lambda=0 here. If +ve then no solution so reset CG
			lambdaSlp1 = 0.0;
			for (int k = 0; k < _N; k++) lambdaSlp1 -= pk[k].dot((*myForces)[k]);
			if (lambdaSlp1 > 0)
			{
				for (int i = 0; i < _3N; i++) pk[i / 3][i % 3] = 0.0; //clear pk
				pk.intoSubtract(gk);
				report << hint << "[NormModeMin::updtConstrainedModes] Reset CG, nabla U^T pk fail." << endl;
				lambdaSlp1 = 0.0;
				for (int k = 0; k < _N; k++) lambdaSlp1 -= pk[k].dot((*myForces)[k]);
			}
			//Use eigenvalues as initial /lambda as solution in linear case
			//Set first value of /lambda to be 1/eigval 
			lambda = curLambdaRef;
			//
			report << debug(7) << "[NormModeMin::updtConstrainedModes] lambd= " << lambda << endl;
			//find force at new position pos+/lambda*pk
			(*myPositions).intoWeightedAdd(lambda, pk);
			calculateForces();
			forceCalc++;
			//no projection as guarenteed to be in comp space by method
			//find slope of PE with /lambda here
			lambdaSlp = 0.0;
			for (int k = 0; k < _N; k++) lambdaSlp -= pk[k].dot((*myForces)[k]);
			//solve for minimum for quadratic fit using two PE vales and the slope with /lambda
			Real a, b, oldLambda, a1;
			oldLambda = lambda;
			a = -((myEnergies->potentialEnergy() - oldPot) / lambda - lambdaSlp) / lambda;
			b = lambdaSlp - 2.0 * a * lambda;
			lambda = -b / (2 * a);
			//calc from first slope (overdefined)
			a1 = (myEnergies->potentialEnergy() - oldPot - lambdaSlp1 * oldLambda) / (oldLambda * oldLambda);
			lambda1 = -lambdaSlp1 / (2 * a1);
			newPE = a * lambda * lambda + b * lambda + oldPot;
			//Diagnostics
			report << debug(7) << "[NormModeMin::updtConstrainedModes] lambdaSlp= " << lambdaSlp << " lambdaSlp1= " << lambdaSlp1 <<
				" lambda= " << lambda << " lambda1= " << lambda1 << " PEdiff= " << myEnergies->potentialEnergy() - oldPot <<
				" newPE2= " << newPE << " newPE1= " << a1 * lambda1 * lambda1 + lambdaSlp1 * lambda1 + oldPot << " lambda Ratio= " << lambda * myPrevMTS->maxEigval << endl;
			if (lambda < 0)
			{
				lambda = lambda1;
				lambNeg = true;
			}
			if (rsCG && lambda > curLambdaRef) lambda = curLambdaRef;
			//Put solution into positions (but remove temporary solution for quadratic fit via oldLambda)
			(*myPositions).intoWeightedAdd(lambda - oldLambda, pk);
			//update total gamma
			totGamma += lambda;
			numGam++;
			//check PE?
			if (lambda1 == 0) lambda1 = 1e-6;
			Real tempRat = lambda / lambda1;
			if (lambNeg || tempRat > 2.0 || tempRat < 0.5 || forcePEcheck)
			{
				calculateForces();
				forceCalc++;
				newPE = myEnergies->potentialEnergy();
				checkedPE = true;
				lambNeg = false;
				if (!forcePEcheck) report << hint << "[NormModeMin::updtConstrainedModes] Re-calc PE." << endl;
			}
			else
			{
				checkedPE = false;
			}
			//test for end, CG reset first
			if ((oldPot - newPE) < 0)
			{
				if (rsCG > 4) report << error << "[NormModeMin::updtConstrainedModes] Minimization failed. Rediagonalization maybe required for the current conformation. Aborting." << endr;
				else
				{
					if (!reDiag)
					{ //allow minimization if mode at angle to sub-space
						(*myPositions).intoWeightedAdd(-lambda, pk); //reset positions
						//update total gamma
						totGamma -= lambda;
						numGam--;
						for (int i = 0; i < _3N; i++) pk[i / 3][i % 3] = 0.0; //and pk
						calculateForces();
						forceCalc++;
						rsCG++;
						curLambdaRef /= 2.0;
						report << hint << "[NormModeMin::updtConstrainedModes] Reset CG, PE fail. rsCG= " << rsCG << endl;
					}
					else
					{
						(*myPositions).intoWeightedAdd(-lambda, pk); //reset positions
						myPrevMTS->mhQ = NULL; //force rediagonalization
						return -1; //flag aborted
					}
				}
			}
			else
			{
				rsCG = 0;
				curLambdaRef = lambdaRef;
			}
			if ((oldPot - newPE) < peLim && !rsCG) break;
			//
		}
		return itr;
	}

	void NormModeMin::run(int numTimesteps)
	{
		if (numTimesteps != -1 && numTimesteps != -2)
		{
			preStepModify();
			postStepModify();
			return;
		}
		else
		{
			preStepModify();
			//restore previous pre-purturbed positions
			if (randforce) (*myPositions).intoWeightedAdd(-randStp, pk);
			//do minimization with local forces
			forceCalc = 0;
			myPrevMTS->itrs = updtConstrainedModes(minLim);
			myPrevMTS->minForceCalc = forceCalc;
			//
			if (randforce && myPrevMTS->itrs != -1)
			{ //add random force, but not if rediagonalizing
				//add random force
				if (totGamma > 0 && numGam > 0)
				{
					randStp = sqrt(totGamma / (Real)numGam);
					//randStp = sqrt((totGamma * Constant::TIMEFACTOR * Constant::TIMEFACTOR)  / (Real)numGam) * Constant::INV_TIMEFACTOR * Constant::INV_TIMEFACTOR;
					pk = myPrevMTS->gaussRandCoordM; //save random purtubation
					(*myPositions).intoWeightedAdd(randStp, myPrevMTS->gaussRandCoordM);
				}
			}
			else
			{
				randStp = 0.0;
			}
			//
			postStepModify();
		}
	}

	void NormModeMin::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("minimlim", Value(minLim, ConstraintValueType::NotNegative()), 0.1, Text("Minimizer target PE difference kcal mole^{-1}")));
		parameters.push_back(Parameter("forcePEcheck", Value(forcePEcheck, ConstraintValueType::NoConstraints()), true, Text("Force PE/calcForces check at end of loop.")));
		parameters.push_back(Parameter("massweight", Value(massWeight, ConstraintValueType::NoConstraints()), true, Text("Position move based on mass weighted forces.")));
		parameters.push_back(Parameter("randforce", Value(randforce, ConstraintValueType::NoConstraints()), true, Text("Add random force")));
		parameters.push_back(Parameter("rediag", Value(reDiag, ConstraintValueType::NoConstraints()), false, Text("Force re-digonalize")));
	}

	STSIntegrator* NormModeMin::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new NormModeMin(values[0], values[1], values[2], values[3], values[4], values[5], fg);
	}
}
