#include "NormModeInt.h"
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

#if defined(HAVE_LAPACK)
#include "LapackProtomol.h"
#else
#if defined(HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

using namespace std;


using namespace ProtoMol::Report;

using std::string;
using std::vector;


#define FIXCOFM	0

namespace ProtoMol
{
	//__________________________________________________ NormModeInt

	const string NormModeInt::keyword("NormModeInt");

	NormModeInt::NormModeInt() : MTSIntegrator(), fixMod(-1), myGamma(-1), mySeed(-1), myTemp(-1), myNVE(0), myBerendsen(0), fDof(0)
	{
		ex0 = NULL;
		tmpFX = NULL;
		tmpC = NULL;
		invSqrtMass = NULL;
		sqrtMass = NULL;
	}

	NormModeInt::NormModeInt(int cycles, int fixmodes, Real gamma, int seed, Real temperature, int nve, Real bertau, int fdof,
	                         ForceGroup* overloadedForces, StandardIntegrator* nextIntegrator)
		: MTSIntegrator(cycles, overloadedForces, nextIntegrator), fixMod(fixmodes),
		  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)), mySeed(seed), myTemp(temperature), myNVE(nve), myBerendsen(bertau * Constant::INV_TIMEFACTOR), fDof(fdof)
	{
		ex0 = NULL;
		tmpFX = NULL;
		tmpC = NULL;
		invSqrtMass = NULL;
		sqrtMass = NULL;
		if (fDof > 6) fDof = 6;
	}


	NormModeInt::~NormModeInt()
	{
		report.precision(5);
		if (numSteps)
			report << hint << "[NormModeInt::~NormModeInt] Minimizer iterations = " << itrs << ", average = " << (float)avItrs / (float)numSteps << ", force calcs = " << minForceCalc << ", average = " << (float)avMinForceCalc / numSteps << endr;
		//
		if (ex0 != NULL) delete ex0;
		if (tmpFX != NULL) delete [] tmpFX;
		if (tmpC != NULL) delete [] tmpC;
		if (invSqrtMass != NULL) delete [] invSqrtMass;
		if (sqrtMass != NULL) delete [] sqrtMass;
	}

	void NormModeInt::initialize(GenericTopology* topo,
	                             Vector3DBlock* positions,
	                             Vector3DBlock* velocities,
	                             ScalarStructure* energies)
	{
		MTSIntegrator::initialize(topo, positions, velocities, energies);
		report << hint << "[NormModeInt::initialize] mpi = " << (myPreviousIntegrator == NULL) << endl;
		//check valid eigenvectors
		if (MTSIntegrator::mhQ == NULL && myPreviousIntegrator == NULL)
			report << error << "No Eigenvectors for NormMode integrator, use NormModeDiag or eigfile." << endr;
		//
		initializeForces();
		//NM initialization
#if FIXCOFM
	Vector3D cofm = centerOfMass(myPositions, myTopo);	
	for( unsigned int i = 0; i < myPositions->size(); i++ ) (*myPositions)[i] -= cofm;
#endif
		//save initial positions
		ex0 = new Vector3DBlock;
		*ex0 = *myPositions;

		//set size-of paramiters
		_N = (int)myPositions->size(); //get N
		_3N = 3 * _N; //used everywhere, must be initialized
		_m = numEigvects; //get m
		//check fixMod
		if (fixMod >= _3N) report << error << "fixedModes = " << fixMod << ", must be less than 3N = " << _3N << "." << endr;
		//if(fixMod > _m) fixMod = _3N-_m;
		if (fixMod < _3N - _m) report << error << "Insufficient eigenvectors: fixedModes = " << fixMod << ", must be greater than or equal to (3N - m) = " << _3N - _m << "." << endr;
		//number of low frequency modes
		_rfM = _3N - fixMod;
		//Degrees of freedom
		myTopo->degreesOfFreedom = _rfM - fDof;
		if (fDof == 0) myTopo->degreesOfFreedom -= 3;
		//temporary mode space variable for intermediate calculations
		tmpC = new double[_3N];
		//mhQ = new double[3*_N*_rfM];	//define matrix hQ
		//
		report << debug(5) << "[NormModeInt::initialize] calc hQ, Q." << endl;
		//put vectorBlock3D read from eigvector file into matrix
		//Since this creates two copies of the matrix in memory we should load it in as a linear array initially
		/*for(int i=0;i<_3N;i++){
			for(int j=0;j<_rfM;j++) mhQ[i + j*_3N] = (*eigvectp)[i/3 + j*_N][i%3];
			}*/
		//define temporary position/force vector
		tmpFX = new double[_3N];
		//setup sqrt masses and their inverse 
		invSqrtMass = new double[_N];
		sqrtMass = new double[_N];
		for (int i = 0; i < _N; i++)
		{
			sqrtMass[i] = sqrt(myTopo->atoms[i].scaledMass);
			if (sqrtMass[i]) invSqrtMass[i] = 1.0 / sqrtMass[i];
		}
		//
		report << debug(5) << "[NormModeInt::initialize] setup velocities." << endl;
		//take initial C velocites from system and remove non-subspace part
		if (MTSIntegrator::mhQ != NULL) subspaceVelocity(myVelocities, myVelocities);
		//****Other**************************************************
		//do first force calculation, and remove non sub-space part
		myEnergies->clear();
		calculateForces();
		if (MTSIntegrator::mhQ != NULL) subspaceForce(myForces, myForces);
		//
		avItrs = 0; //average number of minimizer iterations/force calcs
		avMinForceCalc = 0;
		numSteps = 0; //total steps
		posUpdt.resize(3 * _N);
		//***********************************************************
		//gauss randoms
		gaussRandCoord1.resize(myPositions->size());
		gaussRandCoord2.resize(myPositions->size());
		gaussRandCoordM.resize(myPositions->size());
		//Calculate initial total energy
		orgEnergy = kineticEnergy(myTopo, myVelocities) + myEnergies->potentialEnergy();
	}

	//****Routines to convert between Vector3DBlocks and linear arrays*********************

	//Convert Vector3DBlock formal to linear array for BLAS
	double* NormModeInt::vector3DBlockTOvect(Vector3DBlock* blkDat, double* vecDat)
	{
		for (int i = 0; i < 3 * _N; i++) vecDat[i] = (*blkDat)[i / 3][i % 3];
		return vecDat;
	}

	//Convert linear array from BLAS to Vector3DBlock 
	Vector3DBlock* NormModeInt::vectTOvector3DBlock(double* vecDat, Vector3DBlock* blkDat)
	{
		for (int i = 0; i < 3 * _N; i++) (*blkDat)[i / 3][i % 3] = vecDat[i];
		return blkDat;
	}

	//*************************************************************************************
	//****Projectors for complement sub space and sub space********************************
	//*************************************************************************************

	//Find forces acting outside subspace
	Vector3DBlock* NormModeInt::nonSubspaceForce(Vector3DBlock* force, Vector3DBlock* iPforce)
	{
		//
		vector3DBlockTOvect(force, tmpFX); //put positions-x_0 into linear array
		//calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
		//f'=M^{-1/2}*f
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= invSqrtMass[i / 3];
		//c=hQ^T*M^{-1/2}*f
		char* transA = "T"; // Transpose Q, LAPACK checks only first character N/V
		int m = _3N;
		int n = _rfM;
		int incxy = 1; //sizes
		double alpha = 1.0;
		double beta = 0.0; //multiplyers, see Blas docs.
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transA, &m, &n, &alpha, mhQ, &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	int len_transa = 1;							//length of transA
	dgemv_ (*transA, m, n, alpha, mhQ, m, tmpFX, incxy, beta, tmpC, incxy, len_transa);
#endif
#endif
		//
		//calculate f''=M^{-1/2}*f'-hQc using BLAS
		char* transB = "N";
		alpha = -1.0;
		beta = 1.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transB, &m, &n, &alpha, mhQ, &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	dgemv_ (*transB, m, n, alpha, mhQ, m, tmpC, incxy, beta, tmpFX, incxy, len_transa);
#endif
#endif
		//
		//f'''=M^{1/2}*f''
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= sqrtMass[i / 3];
		//put back into vector3DBlocks
		vectTOvector3DBlock(tmpFX, iPforce);
		//delete temporary array
		return iPforce;
	}

	//Find positions outside subspace
	Vector3DBlock* NormModeInt::nonSubspacePosition(Vector3DBlock* force, Vector3DBlock* iPforce)
	{
		//
		vector3DBlockTOvect(force, tmpFX); //put positions-x_0 into linear array
		//calculate M^{1/2}(I-M^{1/2}\hat{Q}\hat{Q}^TM^{-1/2})M^{-1/2}f using BLAS
		//f'=M^{-1/2}*f
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= sqrtMass[i / 3];
		//c=hQ^T*M^{-1/2}*f
		char* transA = "T"; // Transpose Q, LAPACK checks only first character N/V
		int m = _3N;
		int n = _rfM;
		int incxy = 1; //sizes
		double alpha = 1.0;
		double beta = 0.0; //multiplyers, see Blas docs.
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transA, &m, &n, &alpha, mhQ, &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	int len_transa = 1;							//length of transA
	dgemv_ (*transA, m, n, alpha, mhQ, m, tmpFX, incxy, beta, tmpC, incxy, len_transa);
#endif
#endif
		//
		//calculate f''=M^{-1/2}*f'-hQc using BLAS
		char* transB = "N";
		alpha = -1.0;
		beta = 1.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transB, &m, &n, &alpha, mhQ, &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	dgemv_ (*transB, m, n, alpha, mhQ, m, tmpC, incxy, beta, tmpFX, incxy, len_transa);
#endif
#endif
		//
		//f'''=M^{1/2}*f''
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= invSqrtMass[i / 3];
		//put back into vector3DBlocks
		vectTOvector3DBlock(tmpFX, iPforce);
		//delete temporary array
		return iPforce;
	}

	//Find forces acting inside subspace
	Vector3DBlock* NormModeInt::subspaceForce(Vector3DBlock* force, Vector3DBlock* iPforce)
	{
		//
		vector3DBlockTOvect(force, tmpFX); //put positions-x_0 into linear array
		//calculate M^{1/2}QQ^TM^{-1/2}f using BLAS
		//f'=M^{-1/2}*f
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= invSqrtMass[i / 3];
		//c=Q^T*M^{-1/2}*f
		char* transA = "T"; // Transpose, LAPACK checks only first character N/V
		int m = _3N;
		int n = _rfM - fDof;
		int incxy = 1; //sizes
		double alpha = 1.0;
		double beta = 0.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transA, &m, &n, &alpha, &mhQ[_3N*fDof], &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	int len_transa = 1;	
	dgemv_ (*transA, m, n, alpha, &mhQ[_3N*fDof], m, tmpFX, incxy, beta, tmpC, incxy, len_transa);
#endif
#endif
		//
		//f''=Qc
		char* transB = "N"; /* LAPACK checks only first character N/V */
		alpha = 1.0;
		beta = 0.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transB, &m, &n, &alpha, &mhQ[_3N*fDof], &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	dgemv_ (*transB, m, n, alpha, &mhQ[_3N*fDof], m, tmpC, incxy, beta, tmpFX, incxy, len_transa);
#endif
#endif
		//f'''=M^{1/2}*f''
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= sqrtMass[i / 3];
		//put back into vector3DBlocks
		vectTOvector3DBlock(tmpFX, iPforce);
		//delete temporary array
		return iPforce;
	}

	//Find velocities acting inside subspace
	Vector3DBlock* NormModeInt::subspaceVelocity(Vector3DBlock* force, Vector3DBlock* iPforce)
	{
		//
		vector3DBlockTOvect(force, tmpFX); //put positions-x_0 into linear array
		//calculate M^{-1/2}QQ^TM^{1/2}f using BLAS
		//v'=M^{1/2}*v
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= sqrtMass[i / 3];
		//c=Q^T*M^{-1/2}*v
		char* transA = "T"; //Transpose, LAPACK checks only first character N/V 
		int m = _3N;
		int n = _rfM - fDof;
		int incxy = 1; //sizes
		double alpha = 1.0;
		double beta = 0.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transA, &m, &n, &alpha, &mhQ[_3N*fDof], &m, tmpFX, &incxy, &beta, tmpC, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	int len_transa = 1;	
	dgemv_ (*transA, m, n, alpha, &mhQ[_3N*fDof], m, tmpFX, incxy, beta, tmpC, incxy, len_transa);
#endif
#endif
		//
		//v''=Qc
		char* transB = "N"; /* LAPACK checks only first character N/V */
		alpha = 1.0;
		beta = 0.0;
		//
#if defined(HAVE_LAPACK)
	dgemv_ (transB, &m, &n, &alpha, &mhQ[_3N*fDof], &m, tmpC, &incxy, &beta, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
	dgemv_ (*transB, m, n, alpha, &mhQ[_3N*fDof], m, tmpC, incxy, beta, tmpFX, incxy, len_transa);
#endif
#endif
		//v'''=M^{-1/2}*v''
		for (int i = 0; i < _3N; i++)
			tmpFX[i] *= invSqrtMass[i / 3];
		//put back into vector3DBlocks
		vectTOvector3DBlock(tmpFX, iPforce);
		//delete temporary array
		return iPforce;
	}

	//*************************************************************************************
	//****Normal run routine***************************************************************
	//*************************************************************************************

	void NormModeInt::run(int numTimesteps)
	{
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		Real actTime;

		if (numTimesteps < 1)
			return;

		const unsigned int count = myPositions->size();
		//time calculated in forces! so fix here
		actTime = myTopo->time + numTimesteps * getTimestep();
		//
		//main loop
		for (int i = 0; i < numTimesteps; i++)
		{
			numSteps++;

			//****Real force integrator loop********** 
			//****main loop*************************************
			preStepModify();
			//doHalfKick();
			for (unsigned int ii = 0; ii < count; ++ii)
				(*myVelocities)[ii] += (*myForces)[ii] * h * 0.5 / myTopo->atoms[ii].scaledMass;
			//
			if (!myNVE && !myBerendsen)
			{ //constant temp/Berendsen or langevin?
				drift();
			}
			else
			{ //or energy?
				subspaceVelocity(myVelocities, myVelocities);
				for (unsigned int ii = 0; ii < count; ++ii)
					(*myPositions)[ii] += (*myVelocities)[ii] * h;
			}
			//constraints?
			itrs = 0;
			myEnergies->clear();
			//run minimizer
			if (fixMod) myNextIntegrator->run(-1);
			if (mhQ != NULL)
			{ //not rediagonalize?
				//add number of iterations required
				avItrs += itrs;
				avMinForceCalc += minForceCalc;
				report << debug(5) << "[NormModeInt::run] iterations = " << itrs << " average = " << (float)avItrs / (float)numSteps << " force calcs = " << minForceCalc << " average = " << (float)avMinForceCalc / (float)numSteps << endl;
				//calculate sub space forces
				myEnergies->clear();
				calculateForces();
				subspaceForce(myForces, myForces);
				//
				for (unsigned int i = 0; i < count; ++i)
					(*myVelocities)[i] += (*myForces)[i] * h * 0.5 / myTopo->atoms[i].scaledMass;
				//
				postStepModify();
				//**************************************************
				//
				if (myBerendsen)
				{ //weak Berendsen thermostat?
					Real Tsystem = 2.0 * kineticEnergy(myTopo, myVelocities) / (Constant::BOLTZMANN * myTopo->degreesOfFreedom);
					//reset energy to constant value
					Real efact;
					efact = sqrt(1.0 + (h / myBerendsen) * (myTemp / Tsystem - 1.0));
					for (int k = 0; k < _N; k++) (*myVelocities)[k] *= efact;
				}
			}
			else
			{ //rediagonalize?
				report << debug(5) << "[NormModeInt::run] Re-diagonalizing" << endl;
				myTopo->time = actTime - (i - numTimesteps) * getTimestep();
				if (myPreviousIntegrator == NULL)
					report << error << "[NormModeInt::Run] Re-diagonalization forced with NormModeInt as outermost Integrator. Aborting." << endr;
				return;
			}
			//
		}
		//fix time
		myTopo->time = actTime;
		//
		//Real exhpe, exhcou;
		//expansionPe(&exhpe, &exhcou, 0);
		//report <<debug(4)<<"[NormModeInt::run] PE= "<<myEnergies->potentialEnergy()<<" Tay= "<<exhpe+0.5*exhcou+kineticEnergy(myTopo,myVelocities)<<" exPE= "<<exhpe<<" Coupling= "<<exhcou<<endl;
		//
	}

	//*************************************************************************************
	//****Output int paramiters************************************************************
	//*************************************************************************************

	void NormModeInt::getParameters(vector<Parameter>& parameters) const
	{
		MTSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("fixmodes", Value(fixMod, ConstraintValueType::NotNegative()), 0.0, Text("Number of high frequency modes constrained")));
		parameters.push_back(Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR), ConstraintValueType::NotNegative()), 80.0, Text("Langevin Gamma")));
		parameters.push_back(Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()), 1234, Text("Langevin random seed")));
		parameters.push_back(Parameter("temperature", Value(myTemp, ConstraintValueType::NotNegative()), 300.0, Text("Langevin temperature")));
		parameters.push_back(Parameter("nve", Value(myNVE, ConstraintValueType::NotNegative()), 0, Text("NVE simulation if not 0")));
		parameters.push_back(Parameter("Berendsen", Value(myBerendsen / Constant::INV_TIMEFACTOR, ConstraintValueType::NotNegative()), 0, Text("Berendsen tau in fs")));
		parameters.push_back(Parameter("fdof", Value(fDof, ConstraintValueType::NotNegative()), 6, Text("Fixed degrees of freedom")));
	}

	MTSIntegrator* NormModeInt::doMake(string&, const vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const
	{
		return new NormModeInt(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], fg, nextIntegrator);
	}

	//*************************************************************************************
	//****Langevin Thermostat**************************************************************
	//*************************************************************************************

	// fluctuation using Dr. Skeel's LI scheme which involves a semi-update
	// of velocities and a complete update of positions
	void NormModeInt::drift()
	{
		const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
		const Real tau1 = (1.0 - exp(-myGamma * dt)) / myGamma;
		const Real tau2 = (1.0 - exp(-2 * myGamma * dt)) / (2 * myGamma);
		const Real forceConstant = 2 * Constant::BOLTZMANN * myTemp * myGamma;
		const Real sqrtTau2 = sqrt(tau2);

		//  It is possible that roundoff and/or truncation can make this value < 0.
		//  It should be set to 0 if this happens.
		Real sqrtVal1 = (dt - tau1 * tau1 / tau2);

		if (sqrtVal1 < 0.)
			sqrtVal1 = 0;
		else
			sqrtVal1 = sqrt(sqrtVal1);

		Real sqrtFCoverM = sqrt(forceConstant);/// mass );
		Real langDriftVal = sqrtFCoverM / myGamma;
		Real langDriftZ1 = langDriftVal * (tau1 - tau2) / sqrtTau2;
		Real langDriftZ2 = langDriftVal * sqrtVal1;

		//generate 1st set of _m random force variables and project into sub space
		for (int i = 0; i < _3N; i++)
			gaussRandCoord1[i / 3][i % 3] = randomGaussianNumber(mySeed);//
		for (int i = 0; i < _N; i++)
			gaussRandCoord1[i] *= sqrtMass[i];
		gaussRandCoordM = gaussRandCoord1; //get randoms for minimizer
		subspaceForce(&gaussRandCoord1, &gaussRandCoord1);
		//inner random force
		Real randc = sqrt(2 * Constant::BOLTZMANN * myTemp);
		for (int i = 0; i < _N; i++)
		{
			gaussRandCoordM[i] -= gaussRandCoord1[i];
			gaussRandCoordM[i] *= randc / myTopo->atoms[i].scaledMass;
		}
		//sort mass weighting
		for (int i = 0; i < _N; i++)
			gaussRandCoord1[i] /= myTopo->atoms[i].scaledMass;

		//generate 2nd set of _m random force variables and project into sub space
		for (int i = 0; i < _3N; i++)
			gaussRandCoord2[i / 3][i % 3] = randomGaussianNumber(mySeed);//
		for (int i = 0; i < _N; i++)
			gaussRandCoord2[i] *= sqrtMass[i];
		subspaceForce(&gaussRandCoord2, &gaussRandCoord2);
		for (int i = 0; i < _N; i++)
			gaussRandCoord2[i] /= myTopo->atoms[i].scaledMass;
		//subspaceVelocity(myVelocities, myVelocities);	//this is now looked after by the thermostat
		for (unsigned int i = 0; i < myPositions->size(); i++)
		{
			posUpdt[i] = (gaussRandCoord1[i] * langDriftZ1 + gaussRandCoord2[i] * langDriftZ2 + (*myVelocities)[i]) * tau1;
			// semi-update velocities
			(*myVelocities)[i] = (*myVelocities)[i] * exp(-myGamma * dt) + gaussRandCoord1[i] * sqrtFCoverM * sqrtTau2;
		}
		//subspaceVelocity(myVelocities, myVelocities);	//this is now looked after by the thermostat
		//Check position change is in the sub space and add to positions
		subspaceVelocity(&posUpdt, &posUpdt);
		for (unsigned int i = 0; i < myPositions->size(); i++)
			(*myPositions)[i] += posUpdt[i];
		//fix COM/momentum (not conserved)
		buildMolecularCenterOfMass(myPositions, myTopo);
		buildMolecularMomentum(myVelocities, myTopo);
	}

	void NormModeInt::subSpaceSift()
	{
		//sift current data into subspace
		subspaceVelocity(myVelocities, myVelocities);
		subspaceForce(myForces, myForces);
	}

	Real NormModeInt::expansionPe(Real* hatpe, Real* coupl, int typ)
	{
		Vector3DBlock rlPos, exHatF;

		rlPos = *myPositions; //save real pos
		(*myPositions).intoSubtract(*ex0); //find ex hat pe
		subspaceVelocity(myPositions, myPositions);
		(*myPositions).intoAdd(*ex0);
		myEnergies->clear();
		calculateForces();
		//subspaceForce(myForces, myForces);
		*hatpe = myEnergies->potentialEnergy();
		*myPositions = rlPos; //find coupling term
		(*myPositions).intoSubtract(*ex0); //find ex bar
		nonSubspacePosition(myPositions, myPositions);
		*coupl = 0.0;
		for (unsigned int i = 0; i < myPositions->size() * 3; i++)
			*coupl -= (*myPositions)[i / 3][i % 3] * (*myForces)[i / 3][i % 3];
		*myPositions = rlPos; //set system back to normal
		myEnergies->clear();
		calculateForces();
		subspaceForce(myForces, myForces);
		return *hatpe;
	}

	void NormModeInt::saveState(CheckpointOutputStream& os)
	{
		os << avItrs;
		os << itrs;
		for (int i = 0; i < _3N; i++)
			os << tmpFX[i];
		for (int i = 0; i < _3N; i++)
			os << tmpC[i];
		unsigned int ex0size = ex0->size();
		os << ex0size;
		for (unsigned int i = 0; i < ex0->size(); i++)
			os << (*ex0)[i].x << (*ex0)[i].y << (*ex0)[i].z;
		os << fixMod;
		for (int i = 0; i < _N; i++)
			os << invSqrtMass[i];
		for (int i = 0; i < _N; i++)
			os << sqrtMass[i];
		os << numSteps;
		for (unsigned int i = 0; i < gaussRandCoord1.size(); i++)
			os << gaussRandCoord1[i].x << gaussRandCoord1[i].y << gaussRandCoord1[i].z;
		for (unsigned int i = 0; i < gaussRandCoord2.size(); i++)
			os << gaussRandCoord2[i].x << gaussRandCoord2[i].y << gaussRandCoord2[i].z;
		os << orgEnergy;
		myNextIntegrator->saveState(os);
	}

	void NormModeInt::restoreState(CheckpointInputStream& is)
	{
		is >> avItrs;
		is >> itrs;
		for (int i = 0; i < _3N; i++)
			is >> tmpFX[i];
		for (int i = 0; i < _3N; i++)
			is >> tmpC[i];
		unsigned int ex0size = ex0->size();
		is >> ex0size;
		for (unsigned int i = 0; i < ex0->size(); i++)
			is >> (*ex0)[i].x >> (*ex0)[i].y >> (*ex0)[i].z;
		is >> fixMod;
		for (int i = 0; i < _N; i++)
		{
			Real iSM = invSqrtMass[i];
			is >> iSM;
		}
		for (int i = 0; i < _N; i++)
		{
			Real sM = sqrtMass[i];
			is >> sM;
		}
		is >> numSteps;
		for (unsigned int i = 0; i < gaussRandCoord1.size(); i++)
		{
			is >> gaussRandCoord1[i].x;
			is >> gaussRandCoord1[i].y;
			is >> gaussRandCoord1[i].z;
		}
		for (unsigned int i = 0; i < gaussRandCoord2.size(); i++)
		{
			is >> gaussRandCoord2[i].x;
			is >> gaussRandCoord2[i].y;
			is >> gaussRandCoord2[i].z;
		}
		is >> orgEnergy;
		myNextIntegrator->restoreState(is);
	}
}
