#include "NormModeDiag.h"
#include "Report.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "NormModeInt.h"
#include "GenericTopology.h"

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

namespace ProtoMol {
  //__________________________________________________ NormModeDiag

  const string NormModeDiag::keyword( "NormModeDiag" );

  NormModeDiag::NormModeDiag() : MTSIntegrator(){
		eigVec=NULL;eigVal=NULL;eigIndx=NULL;mhQhQt=NULL;
  }

  NormModeDiag::NormModeDiag(int cycles, int avs,  Real avss, int redi, int rayf, int raya, std::string ray_s, bool raysf, 
	  int mins, Real minl, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), noAvStep(avs), 
		avStep(avss), rediagCount(redi), raylFrequ(rayf), raylAverage(raya), raylFile(ray_s), raylSift(raysf), minSteps(mins), minLim(minl) {
		eigVec=NULL;eigVal=NULL;eigIndx=NULL;mhQhQt=NULL;
		//
		hsn.findForces(overloadedForces);	//find forces and parameters
		if(raylFrequ && raylFrequ < raylAverage) raylFrequ = raylAverage;
  }


  NormModeDiag::~NormModeDiag() 
  {  
	if(eigVec!=NULL && eigAlloc) delete [] eigVec;
	if(eigVal!=NULL) delete [] eigVal;
	if(eigIndx!=NULL) delete [] eigIndx;
	if(mhQhQt!=NULL) delete [] mhQhQt;
  }

  void NormModeDiag::initialize(GenericTopology *topo,
				      Vector3DBlock   *positions,
				      Vector3DBlock   *velocities,
				      ScalarStructure *energies){
    //vars
	sz = 3*positions->size();
	//should check to see if mhQ is large enough and re-use, then fix de-allocation
	if(MTSIntegrator::mhQ != NULL && MTSIntegrator::numEigvects == sz){
		eigVec = mhQ;
		eigAlloc = false;
		validMaxEigv = true;
	}else{
		eigVec = new double[sz*sz];
		eigAlloc = true;
		validMaxEigv = false;
	}
	eigVal = new double[sz];
	eigIndx = new int[sz];
	//initialize Hessian array
	hsn.initialData(sz);
	if(eigVec==NULL || eigVal==NULL || eigIndx==NULL || hsn.hessM==NULL) report << error << "Eigenvector array allocation error."<<endr;
	//set next integrator pointers
	myNextIntegrator->mhQ = MTSIntegrator::mhQ;
	myNextIntegrator->maxEigval = MTSIntegrator::maxEigval;
	//set max eigvects to dof if not file, else file size
	if(mhQ != NULL) myNextIntegrator->numEigvects = MTSIntegrator::numEigvects;
	else myNextIntegrator->numEigvects = sz;
	//
    MTSIntegrator::initialize(topo, positions, velocities, energies);
    initializeForces();
	//setup rediag/rayl counter incase valid
	nextRediag = rediagCount;
	nextRayl = raylFrequ;
	//clear do raylegh quo
	raylDo = false;
	//save positions where diagonalized for checkpoint save (assume I.C. if file)
	diagAt = *myPositions;
	//remove old file and/or define QQ^T
	if(raylFrequ && eigAlloc){
		if(raylFile != ""){
			ofstream myFile;
			myFile.open(raylFile.c_str(),ofstream::out);
			myFile.close();
		}
		if(raylSift){
			mhQhQt = new double[sz*sz];	//define matrix \hat{Q}\hat{Q}^T
		}
	}
  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormModeDiag::run(int numTimesteps) {

    if( numTimesteps < 1 )
      return;
	if(rediagCount && raylDo == false && (int)(myTopo->time/getTimestep()) >= nextRediag){
		nextRediag += rediagCount;
		myNextIntegrator->mhQ = NULL;	//force rediag
	}
	//report <<hint<<"[HessianDiag::run] steps= "<<myTopo->time/getTimestep()+1<<endr;

	//check valid eigenvectors
	if(myNextIntegrator->mhQ == NULL){
		//Diagonalize if no input file
		report <<hint<<"[NormModeDiag::run] Finding diagonalized Hessian."<<endr;
		//save positions where diagonalized for checkpoint save
		diagAt = *myPositions;
		Vector3DBlock tmpVel = *myVelocities;
		//index for sort by absolute
		for(unsigned int i=0;i<sz;i++) eigIndx[i] = i;
		//mass re-weighted hessian
		hsn.clear();
		hsn.evaluate(myPositions, myTopo, true);	//true for mass re-weight;
		//minimized?
		if(minSteps > 0 && validMaxEigv){ //minimizer AND valid maximum eigenvalue 
			int its = SDminimize(minLim, minSteps);
		}
		//Avergaged?
		if(noAvStep > 1){
			Real h = avStep * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();
			calculateForces();
			for(int nos=1;nos<noAvStep;nos++){
				report <<debug(5)<<"[NormModeDiag::run] averaging step = "<<nos<<" step size = "<<avStep<<endr;
				for( unsigned int i = 0; i < count; ++i ) {
					(*myVelocities)[i] += (*myForces)[i]     * 0.5 * h / myTopo->atoms[i].scaledMass;
					(*myPositions)[i]  += (*myVelocities)[i] * h;
				} 
				hsn.evaluate(myPositions, myTopo, true);	//true for mass re-weight;
				calculateForces();
				for( unsigned int i = 0; i < count; ++i ) {
					(*myVelocities)[i] += (*myForces)[i]     * 0.5 * h / myTopo->atoms[i].scaledMass;
				} 
			}
			for(unsigned int i=0 ; i<sz*sz ; i++) hsn.hessM[i] /= (double)noAvStep;	//divide sum of Hessians
		}
		//average or minimized? then reset positions/velocities
		if(minSteps > 0 || noAvStep > 1){
			*myPositions = diagAt;	//back to original positions
			*myVelocities = tmpVel;
		}
		//diagonalize
		int info = diagHessian(eigVec, eigVal);
		if( info == 0 ){
			//find number of -ve eigs
			unsigned int ii;
			for(ii=0;ii<sz-3;ii++) if(eigVal[ii+3] > 0) break; 
			//
			report <<hint<<"[NormModeDiag::run] diagonalized. No. negative eigenvales = "<<ii<<endr;
			absSort();
		}else{
			report << error << "Diagonalization failed."<<endr;
		}
		numEigvects = sz;
		maxEigval = eigVal[sz-1];
		validMaxEigv = true;
		//set next integrator pointers
		//
		myNextIntegrator->mhQ = eigVec;
		myNextIntegrator->maxEigval = maxEigval;
		myNextIntegrator->numEigvects = numEigvects;
		//sift current velocities/forces
		((NormModeInt*)myNextIntegrator)->subSpaceSift();
		//If Rayleigh AND sifted, then calc QQ^T otherwise add eigvals to file/screen
		if(raylFrequ && eigAlloc){
			int rfm = ((NormModeInt*)myNextIntegrator)->_rfM;
			if(raylFile != ""){
				ofstream myFile;
				myFile.open(raylFile.c_str(),ofstream::app);
				myFile.precision(10);
				for(int ofs=1;ofs<5;ofs++)	myFile  << eigVal[rfm - ofs] << " , " << 0 << " , ";
				myFile << endl;
				myFile.close();
			}
			report.precision(10);
			for(int ofs=1;ofs<5;ofs++)	report <<debug(4)<<"[NormModeDiag::run] EigVal "<<ofs<<" ="<<eigVal[rfm - ofs]<<endl;
			if(raylSift){
				//calculate hQhQ^T using BLAS
				char *transA = "N"; char *transB = "T";
				int m = sz; int n = sz; int k = ((NormModeInt*)myNextIntegrator)->_rfM;
				double alpha = 1.0;	double beta = 0.0;
				//
#if defined(HAVE_LAPACK)
				dgemm_ (transA, transB, &m, &n, &k, &alpha, eigVec, &m, eigVec, &m, &beta, mhQhQt, &m);
#else
#if defined(HAVE_SIMTK_LAPACK)
				int len_transa = 1;	int len_transb = 1;
				dgemm_ (*transA, *transB, m, n, k, alpha, eigVec, m,
						eigVec, m, beta, mhQhQt, m, len_transa, len_transb);
#endif
#endif
				//report <<hint<<"[NormModeDiag::initialize] finding QQ^T e11="<<mhQhQt[0]<<" e21="<<mhQhQt[1]<<endr;
			}
		}
	}
	//do rayleigh quotient? Check allocated mhQ or too small for calc!
	if(raylFrequ && raylDo == false && eigAlloc==true && (int)(myTopo->time/getTimestep()) >= nextRayl){
		nextRayl += raylFrequ;
		raylDo = true;
		raylAvCount = 0;
		hsn.clear();	//clear Hessian
	}
	//main loop
	myEnergies->clear();
	if(raylDo){	//if calculating rayleigh do individula steps
		int ii;
		for(ii=0;ii<numTimesteps && raylAvCount < raylAverage;ii++,raylAvCount++){
			myNextIntegrator->run(1);
			hsn.evaluate(myPositions, myTopo, true);	//true for mass re-weight;
		}
		if(ii<numTimesteps) myNextIntegrator->run(numTimesteps-ii);
		if(raylAvCount >= raylAverage){	//do calculation here
			raylDo = false;
			calcRayleigh();
		}
	}else{
		myNextIntegrator->run(numTimesteps);
	}

  }  

  void NormModeDiag::calcRayleigh(){
    int _rfM;
	ofstream myFile;

	//Norm Mode diagnostics-Rayleigh Quotient
	_rfM = ((NormModeInt*)myNextIntegrator)->_rfM;	//get first eigenvector of reduced set
	//
	if(raylFile != ""){
		myFile.open(raylFile.c_str(),ofstream::app);
		myFile.precision(10);
	}
	//Get current Hessian
	for(unsigned int i=0 ; i<sz*sz ; i++) hsn.hessM[i] /= (double)raylAverage;	//divide sum of Hessians
	//set BLAS variables
	char *transA = "N";	
	int m = sz; int n = sz; 
	int incxy = 1;	//sizes
	double alpha = 1.0;	double beta = 0.0;
	if(raylSift){
		//Sift hessian through diagonal
		char *transB = "T";
		//
#if defined(HAVE_LAPACK)
		dgemm_ (transA, transB, &m, &n, &n, &alpha, hsn.hessM, &m,
						mhQhQt, &m, &beta, hsn.hessM, &m);
		dgemm_ (transA, transB, &m, &n, &n, &alpha, mhQhQt, &m,
						hsn.hessM, &m, &beta, hsn.hessM, &m);
#else
#if defined(HAVE_SIMTK_LAPACK)
		int len_transa = 1;	int len_transb = 1;
		dgemm_ (*transA, *transB, m, n, n, alpha, hsn.hessM, m,
						mhQhQt, m, beta, hsn.hessM, m, len_transa, len_transb);
		dgemm_ (*transA, *transB, m, n, n, alpha, mhQhQt, m,
						hsn.hessM, m, beta, hsn.hessM, m, len_transa, len_transb);
#endif
#endif
	}
	//calculate Rayleigh Quo/bound
	double rQ, boundRq;
	//define temporary position/force vector
	double *tmpFX = new double[sz];
	//loop
	for(int ofs = 1;ofs<5;ofs++){
		//
#if defined(HAVE_LAPACK)
		dgemv_ (transA, &m, &n, &alpha, hsn.hessM, &m, &eigVec[sz*(_rfM-ofs)], &incxy, &beta, tmpFX, &incxy);
		rQ = ddot_ (&n, &eigVec[sz*(_rfM-ofs)], &incxy, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
		int len_transa = 1;
		dgemv_ (*transA, m, n, alpha, hsn.hessM, m, &eigVec[sz*(_rfM-ofs)], incxy, beta, tmpFX, incxy, len_transa);
		rQ = ddot_ (n, &eigVec[sz*(_rfM-ofs)], incxy, tmpFX, incxy);
		//yTy = ddot_ (n, &mhQ[_3N*(_rfM-ofs)], incxy, &mhQ[_3N*(_rfM-ofs)], incxy); //normalized so dont need!
#endif
#endif
		//bound
		for(unsigned int i=0;i<sz;i++) hsn.hessM[i+sz*i] -= rQ;
#if defined(HAVE_LAPACK)
		dgemv_ (transA, &m, &n, &alpha, hsn.hessM, &m, &eigVec[sz*(_rfM-ofs)], &incxy, &beta, tmpFX, &incxy);
		boundRq = dnrm2_(&n, tmpFX, &incxy);
#else
#if defined(HAVE_SIMTK_LAPACK)
		dgemv_ (*transA, m, n, alpha, hsn.hessM, m, &eigVec[sz*(_rfM-ofs)], incxy, beta, tmpFX, incxy, len_transa);
		boundRq = dnrm2_(n, tmpFX, incxy);
		//normY = dnrm2_(n, &mhQ[_3N*(_rfM-ofs)], incxy); //normY = sqrt(yTy);//normalized so dont need!
#endif
#endif
		for(unsigned int i=0;i<sz;i++) hsn.hessM[i+sz*i] += rQ; //restore
		//
		report <<debug(4)<<"[NormModeDiag::calcRayleigh] rQ = "<<rQ<<" Bound="<<boundRq<<" _rfM="<<_rfM<<" _3N="<<sz<<endl;
		//
		if(myFile) myFile  << rQ << " , " << boundRq << " , ";// << endl;
	}
	if(myFile){
		myFile << endl;
		myFile.close();
    }
	delete [] tmpFX;
	//
  }

  int NormModeDiag::diagHessian(double *eigVecO, double *eigValO){ //
   double *wrkSp;
   int *isuppz, *iwork;

   wrkSp = new double[26*sz];
   isuppz = new int[2*sz];
   iwork = new int[10*sz];
   //Diagonalize
	char *jobz = "V"; char *range = "A"; char *uplo = "U"; /* LAPACK checks only first character N/V */
	int n = sz;             /* order of coefficient matrix a  */
	int lda = sz;           /* leading dimension of a array*/
	int info;				/* output 0=success */
	double vl = 1.0;
	double vu = 1.0; 
	int il = 1;
	int iu = 1;
	double abstol = 0;
	int m; int ldz = sz; int lwork = 26*sz; /* dimension of work array*/
	int liwork = 10*sz;						/* dimension of int work array*/
	//Recomended abstol for max precision
	char *cmach = "safe min";
	//call LAPACK 
	//	
#if defined( HAVE_LAPACK )
	abstol = dlamch_( cmach);	//find machine safe minimum  
	//
	dsyevr_( jobz, range, uplo, &n, hsn.hessM, &lda, &vl, &vu, &il, &iu, &abstol, &m, eigValO, eigVecO, &ldz, isuppz, 
				wrkSp, &lwork, iwork, &liwork, &info);
#else
#if defined( HAVE_SIMTK_LAPACK )
	int len_cmach = 8;
	int len_jobz = 1; int len_range = 1; int len_uplo = 1;
	abstol = dlamch_( *cmach, len_cmach);	//find machine safe minimum  
	//
	dsyevr_( *jobz, *range, *uplo, n, hsn.hessM, lda, &vl, &vu, &il, &iu, &abstol, m, eigValO, eigVecO, ldz, isuppz, 
				wrkSp, lwork, iwork, &liwork, info, len_jobz, len_range, len_uplo);
	//report <<hint<<"[HessianInt::diagHessian] m="<<m<<" wrkSp="<<wrkSp[0]<<" iwork="<<iwork[0]<<endr;
#endif
#endif

	//delete arrays
    delete [] iwork;
    delete [] isuppz;
    delete [] wrkSp;
	//return status
	return info;
  }

  void NormModeDiag::absSort(){	//sort for absolute magnitude
	unsigned int i;

	//find minimum abs value
	double minEv = fabs(eigVal[0]);
	for(i=1;i<sz;i++){
		if(minEv < fabs(eigVal[i])) break;
		else minEv = fabs(eigVal[i]);
	}
	i--;
	//sort around min
	if(i>0){
		int j = 0;
		eigIndx[j++] = i;
		int negp = i-1;
		unsigned int posp = i+1;
		while(negp >= 0 && posp < sz){
			if(fabs(eigVal[negp]) < fabs(eigVal[posp]))
				eigIndx[j++] = negp--;
			else eigIndx[j++] = posp++;
		}
		while(negp >= 0) eigIndx[j++] = negp--;
	}
	//Sort actual eigenvector array
	double *tmpVect = new double[sz];
	double tmpElt;
	unsigned int ii, k;
	for(i=0;i<sz;i++){
		if( eigIndx[i] != (int)i && eigIndx[i] != -1){		//need to swap?
			for(unsigned int j=0;j<sz;j++){
				tmpVect[j] = eigVec[i*sz+j];
				eigVec[i*sz+j] = eigVec[eigIndx[i]*sz+j];
			}
			eigIndx[i] = -1;								//flag swapped
			ii = i;
			do{
				for(k=0;k<sz && eigIndx[k]!=(int)ii;k++);	//find where tmpVect goes
				if(k==sz || k==ii) break;					//end chain where indeces are equal
				for(unsigned int j=0;j<sz;j++){				//put it there
					tmpElt = tmpVect[j];
					tmpVect[j] = eigVec[k*sz+j];
					eigVec[k*sz+j] = tmpElt;
				}
				eigIndx[k] = -1;							//flag swapped
				ii = k;
			}while(k<sz);
		}
	}
	delete [] tmpVect;
  }

  //Steepest descent minimizer for all modes outside subspace.
  int NormModeDiag::SDminimize(Real peLim, int numlp) {
	int in, itr;
	Real oldPot, lambda, lambdaSlp, lambdaSlp1; 
	Vector3DBlock pk;

	int _N = sz / 3;
	pk.resize(_N);
	//find new forces at start
	calculateForces();
	//
	itr = 0;
	for(in=0;in<numlp;in++){
		itr++;
		//
		report.precision(10);
		report <<debug(5)<<"[NormModeDiag::SDminimize] PE= "<<myEnergies->potentialEnergy()<<endr;
		//****find search direction vector pk
		//set pk to force
		pk.intoWeighted(1.0,*myForces);
		//save PE at /lambda=0
		oldPot = myEnergies->potentialEnergy();
		//find slope of original PE with /lambda=0 here.
		lambdaSlp1 = 0.0;
		for( int k = 0; k < _N; k++ ) lambdaSlp1 -= pk[k].dot((*myForces)[k]);
		//Set first value of /lambda to be 1/eigval times the max value of the eigvector
		lambda = 1.0 / maxEigval;	
		//
		report <<debug(5)<<"[NormModeDiag::SDminimize] lambd= "<<lambda<<endl;
		//find force at new position pos+/lambda*pk
		(*myPositions).intoWeightedAdd(lambda,pk);
		calculateForces();
		//no projection as guarenteed to be in comp space by method
		//find slope of PE with /lambda here
		lambdaSlp = 0.0;
		for( int k = 0; k < _N; k++ ) lambdaSlp -= pk[k].dot((*myForces)[k]);
		//solve for minimum for quadratic fit using two PE vales and the slope with /lambda
		Real a, b, oldLambda;
		oldLambda = lambda;
		a = -((myEnergies->potentialEnergy() - oldPot) / lambda - lambdaSlp) / lambda;
		b = lambdaSlp - 2.0 * a * lambda;
		lambda = -b / (2 * a);
		//Put solution into positions (but remove temporary solution for quadratic fit via oldLambda)
		(*myPositions).intoWeightedAdd(lambda-oldLambda,pk);
		//check
		calculateForces();
		//test for end
		if((oldPot - myEnergies->potentialEnergy()) < peLim) break;
		//
	}
	return itr;
  }

  void NormModeDiag::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("averageSteps", Value(noAvStep,ConstraintValueType::NotNegative()),1,Text("Hessian averaged over number of steps.")));
    parameters.push_back(Parameter("avStepSize",Value(avStep,ConstraintValueType::NotNegative()),1.0,Text("Step size for Hessian averaging.")));
    parameters.push_back(Parameter("reDiagFrequency", Value(rediagCount,ConstraintValueType::NotNegative()),0,Text("Frequency of re-diagonalization (steps).")));
    parameters.push_back(Parameter("raylFrequency", Value(raylFrequ,ConstraintValueType::NotNegative()),0,Text("Frequency of Rayleigh Quationt calculation (steps).")));
    parameters.push_back(Parameter("raylAverage", Value(raylAverage,ConstraintValueType::NotNegative()),1,Text("No. of steps to average Hessian for Rayleigh Quotient.")));
    parameters.push_back(Parameter("raylFile",Value(raylFile,ConstraintValueType::NoConstraints()),false,Text("Rayleigh Quotient output file")));
    parameters.push_back(Parameter("raylSift",Value(raylSift,ConstraintValueType::NoConstraints()),false,Text("Rayleigh, sift Hessian with QQ^T")));
    parameters.push_back(Parameter("minSteps", Value(minSteps,ConstraintValueType::NotNegative()),0,Text("Max. number of minimizer steps.")));
	parameters.push_back(Parameter("minLim",Value(minLim,ConstraintValueType::NotNegative()),1.0,Text("Minimization limit kcal mol^{-1}.")));
  }

  MTSIntegrator* NormModeDiag::doMake(string&, const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormModeDiag(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],fg,nextIntegrator);
  }

  void NormModeDiag::saveState(CheckpointOutputStream& os) {
    os << sz;
    unsigned int theSize = diagAt.size();
    os << theSize;
    for (unsigned int i = 0; i < theSize; i++)
      os << diagAt[i].x << diagAt[i].y << diagAt[i].z;
    for (unsigned int i = 0; i < sz*sz; i++)
      os << hsn.hessM[i];
    for (unsigned int i = 0; i < sz*sz; i++)
      os << eigVec[i];
    for (unsigned int i = 0; i < sz; i++)
      os << eigVal[i];
    for (unsigned int i = 0; i < sz; i++)
      os << eigIndx[i];
    os << eigAlloc;
    os << myNextIntegrator->maxEigval;
    myNextIntegrator->saveState(os);
  }

  void NormModeDiag::restoreState(CheckpointInputStream& is) {
    is >> sz;
    unsigned int theSize = diagAt.size();
    is >> theSize;
    for (unsigned int i = 0; i < theSize; i++) {
      is >> diagAt[i].x;
      is >> diagAt[i].y;
      is >> diagAt[i].z;
    }
    for (unsigned int i = 0; i < sz*sz; i++)
      is >> hsn.hessM[i];
    for (unsigned int i = 0; i < sz*sz; i++)
      is >> eigVec[i];
    for (unsigned int i = 0; i < sz; i++)
      is >> eigVal[i];
    for (unsigned int i = 0; i < sz; i++)
      is >> eigIndx[i];
    is >> eigAlloc;
    is >> myNextIntegrator->maxEigval;
    if (myNextIntegrator->mhQ == NULL)
      myNextIntegrator->mhQ = eigVec;
    myNextIntegrator->restoreState(is);
  }

};
