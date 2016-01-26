#include "NormModeVisual.h"
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

#include "ConfigurationReader.h"
#include "../frontend/InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"

#if defined(HAVE_LAPACK)
#include "LapackProtomol.h"
#else
#if defined(HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;

using std::string;
using std::vector;


namespace ProtoMol {
  //__________________________________________________ NormModeVisual

  const string NormModeVisual::keyword( "NormModeVisual" );

  NormModeVisual::NormModeVisual() : STSIntegrator() {ex0=NULL; eigvalp=NULL; tmpFX=NULL; cPos=NULL; secondPos=NULL;}

  NormModeVisual::NormModeVisual(Real timestep, int fixmodes, std::string  evals, int mstrt, int mend, int cyclest, Real tempr, bool scle, std::string scndp, std::string mfil,
                     ForceGroup *overloadedForces) 
    : STSIntegrator(timestep,overloadedForces), fixMod(fixmodes), evalfile(evals),
        modeStart(mstrt-1), modeEnd(mend-1), cycleSteps(cyclest), tempT(tempr), fScale(scle), secondPosFile(scndp), modeFile(mfil)
  {
      ex0=NULL; eigvalp=NULL; tmpFX=NULL; cPos=NULL; secondPos=NULL;
  }


  NormModeVisual::~NormModeVisual() 
  {  
      if(ex0!=NULL) delete ex0;
      if(eigvalp!=NULL) delete eigvalp;
      if(secondPos!=NULL) delete secondPos;
      if(tmpFX!=NULL) delete [] tmpFX;
      if(cPos!=NULL) delete [] cPos;

  }


  void NormModeVisual::initialize(GenericTopology *topo,
                      Vector3DBlock   *positions,
                      Vector3DBlock   *velocities,
                      ScalarStructure *energies){
  
    STSIntegrator::initialize(topo, positions, velocities, energies);
    initializeForces();
    //NM initialization
    //save initial positions
    ex0 = new Vector3DBlock;
    *ex0 = *myPositions;
    //set size-of paramiters
    _N = (int)myPositions->size(); //get N
    _3N = 3*_N;	//used everywhere, must be initialized
    _m = numEigvects; //get m
    if(fixMod > _m) fixMod = _3N-_m;
    //number of low frequency modes
    _rfM = _3N-fixMod;
    //Degrees of freedom
    myTopo->degreesOfFreedom = _rfM - 6;
    //modes
    cPos = new double[_3N];
    for(int i=0;i<_3N;i++) cPos[i] = 0.0;
    //define temporary position/force vector
    tmpFX = new double[_3N];
    //
    numSteps = 0;	//total steps
    //check modes to be simulated
    if(modeStart<0) modeStart = 0;
    if(modeStart>=(int)numEigvects) modeStart = numEigvects - 1;
    if(modeEnd>=(int)numEigvects) modeEnd = numEigvects - 1;
    if(modeEnd<modeStart) modeEnd = modeStart;
    currMode = modeStart;
    //
    InputPosVel reader;
    eigvalp = new Vector3DBlock;
    if(!reader.open(evalfile))
      report << error << "Can't open eigenvalue file \'"<<evalfile<<"\'."<<endr;
    else
      report << hint << "Eigenvalue file \'"<<evalfile<<"\'."<<endr;
    if(!(reader >> *eigvalp)) //test to see if eigenval.size() matches eigenvect.
      report << error << "Could not parse eigenvalue file \'"<<evalfile<<endr;
    else
      report << hint << "Parsed eigenvalue file \'"<<evalfile<<endr;
    if((*eigvalp).size() != (unsigned int)_N)
      report << error << "Eigenvalue file wrong size ("<<(*eigvalp).size()<<")."<<endr;
    else
      report << hint << "Eigenvalue file correct size ("<<(*eigvalp).size()<<")."<<endr;
    //
    if(secondPosFile != ""){
        // Positions
        if(!reader.open(secondPosFile))
            report << error << "Can't open position file \'"<<secondPosFile<<"\'."<<endr;
        secondPos = new Vector3DBlock;
        vector<PDB::Atom> pdbAtoms;
        if(reader.tryFormat(InputPosVelType::PDB)){
            PDB pdb;
            if(!(reader >> pdb))
                report << error << "Could not parse PDB position file \'"<<secondPosFile<<"\'."<<endr;
            swap(*secondPos,pdb.coords);
            swap(pdbAtoms,pdb.atoms);
        }else if(!(reader >> *secondPos)) report << error << "Could not parse position file \'"<<secondPosFile
                                        <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;
        report << plain << "Using "<<reader.getType()<<" second posfile \'"<<secondPosFile<<"\' ("<<(*secondPos).size()<<")." << endr;
        if((*secondPos).size() == (unsigned int)_N){ //same size?
            //put positions into myPositions
            *myPositions = *secondPos;
            //put the positions into cPos
            for( int i=0; i < _N; i++) (*secondPos)[i] -= (*ex0)[i]; //subtract ex0
            vector3DBlockTOvect(secondPos, tmpFX);	//put second positions-x_0 into linear array
            //calculate M^{-1/2}QQ^TM^{1/2}f using BLAS
            //v'=M^{1/2}*v
            for( int i=0; i < _3N; i++)
                    tmpFX[i] *= sqrt( myTopo->atoms[i/3].scaledMass );
            //c=Q^T*M^{-1/2}*v
            char *transA = "T";							//Transpose, LAPACK checks only first character N/V 
            int m = _3N; int n = _rfM; int incxy = 1;	//sizes
            double alpha = 1.0;	double beta = 0.0;
            //
#if defined(HAVE_LAPACK)
            dgemv_ (transA, &m, &n, &alpha, mhQ, &m, tmpFX, &incxy, &beta, cPos, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
            int len_transa = 1;	
            dgemv_ (*transA, m, n, alpha, mhQ, m, tmpFX, incxy, beta, cPos, incxy, len_transa);
#endif
#endif
            if(modeFile != ""){
                //Output modes for analysis
                ofstream myFile;
                myFile.open(modeFile.c_str(),ofstream::out);
                myFile.precision(10);
                for(int ii=0;ii<_3N;ii++) myFile << cPos[ii] << endl;
                //close file
                myFile.close();
            }
        }else{
            report << error << "Second position file different dimension ( "<<(*secondPos).size()<<" ) to system ( "<<_N<< " )."<<endr;
        }
    }
    //####hack to display current mode on web application
    numEigvects = modeStart+1;
    //####
    //
  }

  //Convert Vector3DBlock formal to linear array for BLAS
  double* NormModeVisual::vector3DBlockTOvect(Vector3DBlock* blkDat, double* vecDat){
      for( int i=0; i<3*_N; i++) vecDat[i] = (*blkDat)[i/3][i%3];
      return vecDat;
  }

  //Convert linear array from BLAS to Vector3DBlock 
  Vector3DBlock* NormModeVisual::vectTOvector3DBlock(double* vecDat, Vector3DBlock* blkDat){
      for( int i=0; i<3*_N; i++) (*blkDat)[i/3][i%3] = vecDat[i];
      return blkDat;
  }

  //Project from subspace to 3D space
  Vector3DBlock* NormModeVisual::subspaceProj(double *tmpC, Vector3DBlock * iPos){
    //
    char *transA = "N";							// Transpose Q, LAPACK checks only first character N/V
    int m = _3N; int n = _rfM; int incxy = 1;	//sizes
    double alpha = 1.0;	double beta = 0.0;		//multiplyers, see Blas docs.
    //
#if defined(HAVE_LAPACK)
    dgemv_ (transA, &m, &n, &alpha, mhQ, &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
    int len_transa = 1;							//length of transA
    dgemv_ (*transA, m, n, alpha, mhQ, m, tmpC, incxy, beta, tmpFX, incxy, len_transa);
#endif
#endif
    //
    for( int i=0; i < _3N; i++)
            tmpFX[i] /= sqrt( myTopo->atoms[i/3].scaledMass );
    //put back into vector3DBlocks
    vectTOvector3DBlock(tmpFX, iPos);
    //add ex0
    for( int i=0; i < _N; i++)
        (*iPos)[i] += (*ex0)[i];
    //delete temporary array
    return iPos;
  }

  void NormModeVisual::run(int numTimesteps) {
    //Real h = getTimestep();
    Real tempFrq;
    Real tempKt = sqrt(2.0 * tempT * Constant::BOLTZMANN);
    //Real totalT;

    if( numTimesteps < 1 )
      return;

    //main loop
    for( int i = 0; i < numTimesteps; i++ ) {
      numSteps++;

      //set mode
      if(!(numSteps%cycleSteps)){
          cPos[currMode++] = 0.0;
          *myPositions = *ex0;
          if(currMode>modeEnd) currMode = modeStart;
          numEigvects = currMode + 1;
      }
      //****Analytical mode integrator loop*****
      tempFrq = sqrt(fabs((*eigvalp)[currMode/3][currMode%3]));
      cPos[currMode] = tempKt * sin((double)(numSteps%cycleSteps)/(double)cycleSteps*2.0*3.14159265);
      if(fScale) cPos[currMode] /= tempFrq;
      else cPos[currMode] /= sqrt(tempFrq);
      //cPos[j] = tempKt * sin(totalT*tempFrq) / tempFrq;
      subspaceProj(cPos, myPositions);
    }
    //calculateForces();
  }  

  void NormModeVisual::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("fixmodes",Value(fixMod,ConstraintValueType::NotNegative()),0.0,Text("Number of high frequency modes constrained")));
    parameters.push_back(Parameter("eigvalFile",Value(evalfile,ConstraintValueType::NoConstraints()),std::string(""),Text("Eigenvalue filename")));
    parameters.push_back(Parameter("modeStart",Value(modeStart+1,ConstraintValueType::NotNegative()),1,Text("First mode")));
    parameters.push_back(Parameter("modeEnd",Value(modeEnd+1,ConstraintValueType::NotNegative()),1,Text("Last mode")));
    parameters.push_back(Parameter("cycleSteps",Value(cycleSteps,ConstraintValueType::NotNegative()),100,Text("Number of steps per mode cycle")));
    parameters.push_back(Parameter("temperature",Value(tempT,ConstraintValueType::NotNegative()),300.0,Text("Temperature")));
    parameters.push_back(Parameter("frequScale",Value(fScale,ConstraintValueType::NoConstraints()),true,Text("Scale for frequency")));
    parameters.push_back(Parameter("secondPosFile",Value(secondPosFile,ConstraintValueType::NoConstraints()),std::string(""),Text("Second position filename")));
    parameters.push_back(Parameter("modeFile",Value(modeFile,ConstraintValueType::NoConstraints()),std::string(""),Text("mode output filename")));
  }

  STSIntegrator* NormModeVisual::doMake(string& , const vector<Value>& values,ForceGroup* fg)const{
    return new NormModeVisual(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],fg);
  }

}

