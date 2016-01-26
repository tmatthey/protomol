/* -*- c++ -*- */
#ifndef ISGNONBONDEDFULLEWALDSYSTEMFORCE_H
#define ISGNONBONDEDFULLEWALDSYSTEMFORCE_H

#include "SystemForce.h"
#include "Parallel.h"
#include "iSGNonbondedFullEwaldSystemForceBase.h"
#include "iSGCoulombForce.h"
#include "Topology.h"
#include "ScalarStructure.h"
#include "PeriodicBoundaryConditions.h"
#include "Timer.h"
#include "mathutilities.h"
#include "simpleTypes.h"
#include "pmconstants.h"
#include "CutoffSwitchingFunction.h"

using namespace ProtoMol::Report;

//#define DEBUG_EWALD_TIMING
//#define DEBUG_EWALD_ENERGIES
//#define USE_EWALD_EXACT_ERF
//#define USE_EWALD_NO_SINCOS_TABLE

namespace ProtoMol {

  //_________________________________________________________________ iSGNonbondedFullEwaldSystemForce
  
  template<class TBoundaryConditions, 
	   class TCellManager,
	   bool  real,
	   bool  reciprocal,
	   bool  correction,
	   class TSwitchingFunction>
  class iSGNonbondedFullEwaldSystemForce: public SystemForce, private iSGNonbondedFullEwaldSystemForceBase {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef Topology<TBoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef typename RealTopologyType::Enumerator::CellPair CellPairType;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGNonbondedFullEwaldSystemForce();
    iSGNonbondedFullEwaldSystemForce(Real alpha, Real accuracy, Real expansionFactor);

    virtual ~iSGNonbondedFullEwaldSystemForce();
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGNonbondedFullEwaldSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

  private:
    void initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions);

    void realTerm(const RealTopologyType* realTopo,
		  const Vector3DBlock* positions, 
		  Vector3DBlock* forces, 
		  ScalarStructure* energies,
		  Real& realEnergy,
		  Real& realDeltaMu,
		  unsigned int n);

    void reciprocalTerm(const RealTopologyType* realTopo,
			const Vector3DBlock* positions, 
			Vector3DBlock* forces, 
			ScalarStructure* energies,
			Real& reciprocalEnergy,
                        Real& reciprocalDeltaMu,
			unsigned int from, unsigned int to);

    void correctionTerm(const RealTopologyType* realTopo,
			const Vector3DBlock* positions, 
			Vector3DBlock* forces, 
			ScalarStructure* energies,
			Real& intraMolecularEnergy,
			Real& intraMolecularDeltaMu,
			unsigned int from, unsigned int to);

    void surfaceDipoleTerm(const RealTopologyType* realTopo,
			   const Vector3DBlock* positions, 
			   Vector3DBlock* forces, 
			   ScalarStructure* energies,
			   Real& surfaceDipoleEnergy);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual void evaluate(const GenericTopology*, 
			  const Vector3DBlock*, 
			  Vector3DBlock*, 
			  ScalarStructure*);

    virtual void parallelEvaluate(const GenericTopology*, 
				  const Vector3DBlock*, 
				  Vector3DBlock*, 
				  ScalarStructure*);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual unsigned int numberOfBlocks(const GenericTopology*,const Vector3DBlock*);
    virtual std::string getKeyword() const{return keyword;}
    virtual void uncache(){myCached=false;};
  private:
    virtual Force* doMake(std::string& errMsg, std::vector<Value>) const;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return (TBoundaryConditions::PERIODIC ? 2:3);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
    bool myCached;

    std::vector<Vector3D> myK;  // Reciprocal lattice vector (h*2PI/Lx,k*2PI/Ly,l*2PI/Lz)
    std::vector<Real> myKSquared;  // Squared norm of the reciprocal lattice vector
#ifndef USE_EWALD_NO_SINCOS_TABLE
    std::vector<TripleInt> myKInt;      // Reciprocal lattice vector (h,k,l)
#endif

    Real myExpansionFactor;
    Real myTRate;             //
    Real myAccuracy;

    Real myAlpha;             // 
    bool myAlphaDefault;

    Real myAlphaSquared;
    Real myAlphaSquaredr;
    Real my2AlphaPI;
    Real myFac;

    Real myRc;                // Cutoff real term
    Real myRcSquared;         // Cutoff squared real term
    Real myKc;                // Cutoff reciprocal term

    int myHmax;
    int myKmax;
    int myLmax;
    int myHKLmax;

    Real* mySinCosA;          // Look up tables
    Real* myLastSinCos;
    Vector3D* mySinCos;

    Real* myDMU_SinCosA;
    Real* myDMU_LastSinCos;
    Vector3D* myDMU_SinCos;

    Real myLX, myLY, myLZ;
    Real myLXr,myLYr,myLZr;
    Real myV,myVr;
    Vector3D myOrigin;

    Real myPointSelfEnergy;   // Precomputed energy terms
    Real myPointSelfDeltaMu;
    Real myChargedSystemEnergy;
    Real myChargedSystemDeltaMu;
                                         
    PeriodicBoundaryConditions boundaryConditions;

    TSwitchingFunction switchingFunction;
    EnumeratorType enumerator;
    std::vector<Vector3D> myLattice;
    unsigned int myOldAtomCount; // Keep track of old values and hope
    int myOldHKLmax;             // we do not need to reallocate memory ...
#if defined(DEBUG_EWALD_TIMING)
    Timer myReal;
    Timer myReciprocal;
    Timer myIntra;
    Timer mySurface;
#endif
  };

  //______________________________________________________________________ INLINES
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
				   TCellManager,
				   real,
				   reciprocal,
				   correction,
				   TSwitchingFunction>::iSGNonbondedFullEwaldSystemForce()
				     : SystemForce(),
				       myCached(false),
				       myExpansionFactor(0.0),
				       myTRate(0.0),
				       myAccuracy(0.0),
				       myAlpha(-1.0),
				       myAlphaDefault(true),
				       myRc(0.0),
				       myKc(0.0),
				       mySinCosA(NULL),
				       myLastSinCos(NULL),
				       mySinCos(NULL),
				       myDMU_SinCosA(NULL),
				       myDMU_LastSinCos(NULL),
				       myDMU_SinCos(NULL),
				       myV(-1.0),
				       myOldAtomCount(0),
				       myOldHKLmax(0){
#if defined(DEBUG_EWALD_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }

  //________________________________________________________________constructor
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
				   TCellManager,
				   real,reciprocal,
				   correction,
				   TSwitchingFunction>::iSGNonbondedFullEwaldSystemForce(Real alpha, Real accuracy, Real expansionFactor)
				     : SystemForce(),
				       myCached(false),
				       myExpansionFactor(expansionFactor),
				       myTRate(0.0),
				       myAccuracy(accuracy),
				       myAlpha(alpha),
				       myAlphaDefault(alpha<= 0.0),
				       myRc(0.0),
				       myKc(0.0),
				       mySinCosA(NULL),
				       myLastSinCos(NULL),
				       mySinCos(NULL),
				       myDMU_SinCosA(NULL),
				       myDMU_LastSinCos(NULL),
				       myDMU_SinCos(NULL),
				       myV(-1.0),
				       myOldAtomCount(0),
				       myOldHKLmax(0) {
#if defined(DEBUG_EWALD_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }
  
  //________________________________________________________________destructor
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
				   TCellManager,
				   real,
				   reciprocal,
				   correction,
				   TSwitchingFunction>::~iSGNonbondedFullEwaldSystemForce(){
    if(reciprocal){
      delete [] mySinCosA;
      delete [] mySinCos;
      delete [] myLastSinCos;
      delete [] myDMU_SinCosA;
      delete [] myDMU_SinCos;
      delete [] myDMU_LastSinCos;
    }
#if defined(DEBUG_EWALD_TIMING)
    if(boundaryConditions.getVolume() > EPSILON){
      report.setf(std::ios::showpoint|std::ios::fixed);
      report << allnodes << plain <<"Timing ("
	     <<Parallel::getId()<<") Ewald:"
	     <<": real:" <<myReal.getTime().getProcessTime()<<"[s]"
	     <<", reciprocal:"<<myReciprocal.getTime().getProcessTime()<<"[s]"
	     <<", intra:"<<myIntra.getTime().getProcessTime()<<"[s]"
	     <<", dipole:"<<mySurface.getTime().getProcessTime()<<"[s]"
	     <<", Rc="<<myRc
	     <<", Kc="<<myKc
	     <<", Kn="<<myK.size()
	     <<", n="<<myLattice.size()
	     <<", alpha="<<myAlpha
	     <<", accuracy="<<myAccuracy<<"."<<endr;
    }
#endif
  }
  
  //________________________________________________________________evaluate function
  template <class TBoundaryConditions, 
	    class TCellManager, 
	    bool  real, 
	    bool  reciprocal, 
	    bool  correction, 
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::evaluate(const GenericTopology* topo, 
								      const Vector3DBlock* positions, 
								      Vector3DBlock* forces, 
								      ScalarStructure* energies) {

    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);  
    
    // Initialize data members and precompute tables & short cuts
    if(!myCached)
      initialize(realTopo,positions);
    
    // Intra-molecular and surface diplol term
    Real intraMolecularEnergy = 0.0;
    Real surfaceDipoleEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
    Real chargedSystemEnergy = 0.0;
    Real intraMolecularDeltaMu = 0.0;
    Real pointSelfDeltaMu = 0.0;
    Real chargedSystemDeltaMu = 0.0;

    if(correction){
      correctionTerm(realTopo,positions,forces,energies,
		     intraMolecularEnergy,intraMolecularDeltaMu,
		     0,realTopo->exclusions.getTable().size());
      pointSelfEnergy      = myPointSelfEnergy;
      chargedSystemEnergy  = myChargedSystemEnergy;
      pointSelfDeltaMu     = myPointSelfDeltaMu;
      chargedSystemDeltaMu = myChargedSystemDeltaMu;
      (*energies)[ScalarStructure::VIRIALXX] += myChargedSystemEnergy;
      (*energies)[ScalarStructure::VIRIALYY] += myChargedSystemEnergy;
      (*energies)[ScalarStructure::VIRIALZZ] += myChargedSystemEnergy;
      if(false)
	surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
    }

    // Real-space term
    Real realEnergy = 0.0;
    Real realDeltaMu = 0.0;
    if(real) {
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myRc);

      realTerm(realTopo,positions,forces,energies,realEnergy,realDeltaMu,realTopo->cellLists.size());
    }  

    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    Real reciprocalDeltaMu = 0.0;
    if(reciprocal)
      reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy,reciprocalDeltaMu,
		     0,myK.size());

    // Sum of all energy terms
    // Sum of all energy terms
    Real e = 
      realEnergy+
      reciprocalEnergy+
      intraMolecularEnergy+
      surfaceDipoleEnergy+
      pointSelfEnergy+
      chargedSystemEnergy;
    Real d_mu = 
      realDeltaMu+
      reciprocalDeltaMu+
      intraMolecularDeltaMu+
      pointSelfDeltaMu+
      chargedSystemDeltaMu;

    (*energies)[ScalarStructure::COULOMB] += e;
    (*energies)[ScalarStructure::COULOMB_DELTAMU] += d_mu;

#if defined(DEBUG_EWALD_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << plain <<"Ewald: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   << endr;
    report << plain <<"Ewald: point="<<myPointSelfDeltaMu
	   <<", charged="<<myChargedSystemDeltaMu
	   <<", real="<<realDeltaMu
	   <<", reciprocal="<<reciprocalDeltaMu
	   <<", intra="<<intraMolecularDeltaMu
	   <<", total="<<d_mu
	   << endr;
#endif
  }

  //________________________________________________________________parallelEvaluateFunction
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::parallelEvaluate(const GenericTopology* topo, 
									      const Vector3DBlock* positions, 
									      Vector3DBlock* forces, 
									      ScalarStructure* energies) {

    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);  
    
    // Initialize data members and precompute tables & short cuts
    if(!myCached)
      initialize(realTopo,positions);
    
    // Intra-molecular and surface diplol term
    Real intraMolecularEnergy = 0.0;
    Real surfaceDipoleEnergy = 0.0;
    Real pointSelfEnergy = 0.0;
    Real chargedSystemEnergy = 0.0;
    Real intraMolecularDeltaMu = 0.0;
    Real pointSelfDeltaMu = 0.0;
    Real chargedSystemDeltaMu = 0.0;

    if(correction) {
      unsigned int n = realTopo->exclusions.getTable().size();
      if(n > 0) {
	unsigned int count = std::min(n,static_cast<unsigned int>(Parallel::getAvailableNum()));    
	for(unsigned int i = 0;i<count;i++) {
	  if(Parallel::next()) {
	    unsigned int to = (n*(i+1))/count;
	    if(to > n)
	      to = n;
	    unsigned int from = (n*i)/count;
	    correctionTerm(realTopo,positions,forces,energies,
			   intraMolecularEnergy,intraMolecularDeltaMu,from,to);
	  }
	}
      }
      if(Parallel::getAvailableId() == 0){
	pointSelfEnergy      = myPointSelfEnergy;
	chargedSystemEnergy  = myChargedSystemEnergy;
	pointSelfDeltaMu     = myPointSelfDeltaMu;
	chargedSystemDeltaMu = myChargedSystemDeltaMu;
	(*energies)[ScalarStructure::VIRIALXX] += myChargedSystemEnergy;
	(*energies)[ScalarStructure::VIRIALYY] += myChargedSystemEnergy;
	(*energies)[ScalarStructure::VIRIALZZ] += myChargedSystemEnergy;
	if(false)
	  surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
      }
    }
    
    // Real-space term
    Real realEnergy = 0.0;
    Real realDeltaMu = 0.0;
    if(real){
      realTopo->updateCellLists(positions);
      enumerator.initialize(realTopo, myRc);
      unsigned int n = realTopo->cellLists.size();
      unsigned int count = Parallel::getNumberOfPackages(n);

      for(unsigned int i = 0;i<count;i++){
	unsigned int l = (n*(i+1))/count - (n*i)/count;
	if(Parallel::next())
	  realTerm(realTopo,positions,forces,energies,realEnergy,realDeltaMu,l); 
	else
	  enumerator.nextNewPair(l);	
      }
    }  
    
    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    Real reciprocalDeltaMu = 0.0;
    if(reciprocal){   
      unsigned int count = std::min(static_cast<unsigned int>(myK.size()),static_cast<unsigned int>(Parallel::getAvailableNum()));
      
      for(unsigned int i = 0;i<count;i++)
	if(Parallel::next())	
	  reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy,reciprocalDeltaMu,
			 (myK.size()*i)/count,(myK.size()*(i+1))/count);   
      
    }
    
    // Sum of all energy terms
    Real e = 
      realEnergy+
      reciprocalEnergy+
      intraMolecularEnergy+
      surfaceDipoleEnergy+
      pointSelfEnergy+
      chargedSystemEnergy;

    Real d_mu = 
      realDeltaMu+
      reciprocalDeltaMu+
      intraMolecularDeltaMu+
      pointSelfDeltaMu+
      chargedSystemDeltaMu;

    (*energies)[ScalarStructure::COULOMB] += e;
    (*energies)[ScalarStructure::COULOMB_DELTAMU] += d_mu;

#if defined(DEBUG_EWALD_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << allnodes << plain <<"Ewald: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   << endr;
    report << allnodes << plain <<"Ewald: point="<<myPointSelfDeltaMu
	   <<", charged="<<myChargedSystemDeltaMu
	   <<", real="<<realDeltaMu
	   <<", reciprocal="<<reciprocalDeltaMu
	   <<", intra="<<intraMolecularDeltaMu
	   <<", total="<<d_mu
	   << endr;
#endif
  }

  //________________________________________________________________numberOfBlocks function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  unsigned int iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
						TCellManager,
						real,
						reciprocal,
						correction,
						TSwitchingFunction>::numberOfBlocks(const GenericTopology* topo,
										    const Vector3DBlock* positions){
    unsigned int n = 0;
    if(correction)
      n +=std::min(static_cast<int>(topo->exclusions.getTable().size()),static_cast<int>(Parallel::getAvailableNum()));
    
    if(reciprocal){
      if(!myCached)
	initialize(dynamic_cast<const RealTopologyType*>(topo),positions);
      
      n += std::min(static_cast<unsigned int>(myK.size()),static_cast<unsigned int>(Parallel::getAvailableNum()));
    }
    
    if(real){
      const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
      realTopo->updateCellLists(positions);    
      n += Parallel::getNumberOfPackages(realTopo->cellLists.size());
    }

    return n;
  }
  
  //________________________________________________________________initialize function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions) {
    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint)
      Report::report << Report::donthint;
    
    if(TBoundaryConditions::VACUUM) {
      // We are have a non-periodic case. Do an expansion of the simulation 
      // box by the expansion factor, if need
      realTopo->updateCellLists(positions);
      Vector3D d2 = (boundaryConditions.getMax()-boundaryConditions.getMin());
      Vector3D d0 = d2/myExpansionFactor;
      Vector3D d1 = realTopo->max-realTopo->min;
      if (fabs(1-d0.x/d1.x) > 0.1 || fabs(1-d0.y/d1.y) > 0.1 || fabs(1-d0.y/d1.y) > 0.1 || 
	  d2.x <= d1.x || d2.y <= d1.y || d2.z <= d1.z ){
	// The boundaries are to big or small, adjust
	Vector3D e = (realTopo->max-realTopo->min)*myExpansionFactor;
	Vector3D origin = realTopo->min*myExpansionFactor + e*0.5;
	boundaryConditions.set(Vector3D(e.x,0,0), Vector3D(0,e.y,0), Vector3D(0,0,e.z), origin);
	report << hint << "The boundaries for Ewald force evaluation were re-sized to "<<boundaryConditions.getMin()<<"-"<<boundaryConditions.getMax()<<"."<<endr;
      }
      else {
	// The particles are just inside +-10% our previous simulation box. No update needed
	myV = boundaryConditions.getVolume();
	return;      
      }
    }
    else {
      boundaryConditions.set(realTopo->boundaryConditions.e1(),
			     realTopo->boundaryConditions.e2(),
			     realTopo->boundaryConditions.e3(),
			     realTopo->boundaryConditions.origin());
    }
    
    if(!boundaryConditions.isOrthogonal())
      report << error << "[NonbondedFullEwaldSystemForce::initialize] Not orthogonal, aborting."<<endr;

    const unsigned int atomCount = realTopo->atoms.size();
    
    
    myTRate     = 5.5;     // From Moldy, rate between real and reciprocal
    //myAccuracy  = 0.00001; // From Moldy, accuracy
    //myAccuracy  = 1.e-6; // From NAMD, accuracy
  
    // Dimension of the simulation box
    myLX   = boundaryConditions.e1().x;
    myLY   = boundaryConditions.e2().y;
    myLZ   = boundaryConditions.e3().z;
    myV    = boundaryConditions.getVolume();
    myLXr  = boundaryConditions.e1r().x;
    myLYr  = boundaryConditions.e2r().y;
    myLZr  = boundaryConditions.e3r().z;
    myVr   = 1.0/myV;
    myOrigin = boundaryConditions.origin();

    // Short cuts
    if(myAlphaDefault)
      myAlpha         = sqrt(M_PI)*pow(myTRate*atomCount/(power<2>(myV)),1.0/6.0);
    myAlphaSquared  = myAlpha*myAlpha;
    myAlphaSquaredr = 1.0/myAlphaSquared;
    my2AlphaPI      = 2.0*myAlpha/sqrt(M_PI);
    Real p          = -log(myAccuracy);
    myRc            = sqrt(p)/myAlpha;
    myRcSquared     = myRc*myRc;
    myKc            = 2.0*myAlpha*sqrt(p);
    myFac           =  1.0/(4.0*myAlpha*myAlpha);

    switchingFunction=TSwitchingFunction(myRc);

    // Reciprocal part
    // Build the lattice k-vectors 2*PI(n0/Lx,n1/Ly,n2/Lz)
    // Maximum values of h, k, l  s.t. |k| < myKc
    myHmax = (int)floor(myKc/(2.0*M_PI)*myLX);
    myKmax = (int)floor(myKc/(2.0*M_PI)*myLY);
    myLmax = (int)floor(myKc/(2.0*M_PI)*myLZ);
    myHKLmax = std::max(2,std::max(myHmax,std::max(myKmax,myLmax))+1);
    myK.clear();
    myKSquared.clear();
#ifndef USE_EWALD_NO_SINCOS_TABLE
    myKInt.clear();
#endif
    int lastH = Constant::MAX_INT;
    int lastK = Constant::MAX_INT;
    int misses = 0;
    Real kcSquared = myKc*myKc;
    if(reciprocal){
      for(int h = 0; h <= myHmax; h++){
	Real kx = 2.0*M_PI*h*myLXr;
	for(int k = (h==0 ? 0 : -myKmax); k <= myKmax; k++){
	  Real ky = 2.0*M_PI*k*myLYr;
	  for(int l = (h==0 && k==0 ? 1 : -myLmax); l <= myLmax; l++) {
	    Real kz = 2.0*M_PI*l*myLZr;
	    if(kx*kx + ky*ky + kz*kz < kcSquared) {
	      myK.push_back(Vector3D(kx,ky,kz));
	      myKSquared.push_back(kx*kx + ky*ky + kz*kz);
#ifndef USE_EWALD_NO_SINCOS_TABLE
	      TripleInt tmp(h,k,l);
	      myKInt.push_back(tmp);
	      if(lastH != h || lastK != k)
		misses++;
	      lastH = h;
	      lastK = k;
#endif
	    }
          }
	}
      }
      
      // Allocate memory for our tables
      if(mySinCosA != NULL && (atomCount != myOldAtomCount)){
	delete [] mySinCosA;
	mySinCosA = NULL;
      }
      
      if(myDMU_SinCosA != NULL && (atomCount != myOldAtomCount)){
	delete [] myDMU_SinCosA;
	myDMU_SinCosA = NULL;
      }
      
      if(mySinCos != NULL && (atomCount != myOldAtomCount || myHKLmax != myOldHKLmax)){
	delete [] mySinCos;
	mySinCos = NULL;
      }
      
      if(myDMU_SinCos != NULL && (atomCount != myOldAtomCount || myHKLmax != myOldHKLmax)){
	delete [] myDMU_SinCos;
	myDMU_SinCos = NULL;
      }
      
      if(myLastSinCos != NULL && (atomCount != myOldAtomCount)){
	delete [] myLastSinCos;
	myLastSinCos = NULL;
      }
      
      if(myDMU_LastSinCos != NULL && (atomCount != myOldAtomCount)){
	delete [] myDMU_LastSinCos;
	myDMU_LastSinCos = NULL;
      }
      
      if(mySinCosA == NULL)
        mySinCosA    = new Real[2*atomCount];
      if(mySinCosA == NULL)
        report << error << "[ISGNonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "<<2*atomCount*sizeof(Real)<<" bytes."<<endr;

      if(myDMU_SinCosA == NULL)
        myDMU_SinCosA    = new Real[2*atomCount];
      if(myDMU_SinCosA == NULL)
        report << error << "[ISGNonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "<<2*atomCount*sizeof(Real)<<" bytes."<<endr;

#ifndef USE_EWALD_NO_SINCOS_TABLE
      if(mySinCos == NULL)
        mySinCos     = new Vector3D[2*atomCount*myHKLmax];
      if(myLastSinCos == NULL)
        myLastSinCos = new Real[2*atomCount];
      if(mySinCos == NULL || myLastSinCos == NULL)
        report << error << "[ISGNonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "
	       <<2*atomCount*myHKLmax *sizeof(Real)*3 + 2*atomCount *sizeof(Real)<<" bytes."<<endr;
      
      if(myDMU_SinCos == NULL)
        myDMU_SinCos     = new Vector3D[2*atomCount*myHKLmax];
      if(myDMU_LastSinCos == NULL)
        myDMU_LastSinCos = new Real[2*atomCount];
      if(myDMU_SinCos == NULL || myDMU_LastSinCos == NULL)
        report << error << "[ISGNonbondedFullEwaldSystemForce::evaluate] Not enough memory, requesting "
	       <<2*atomCount*myHKLmax *sizeof(Real)*3 + 2*atomCount *sizeof(Real)<<" bytes."<<endr;
      
      myOldHKLmax = myHKLmax;
#endif
      myOldAtomCount = atomCount;
    }
    
    //
    // Point self-energy and chemical potential difference
    //
    myPointSelfEnergy = 0.0;
    myPointSelfDeltaMu = 0.0;
    if(correction) {
      Real q = 0.0;
      Real DeltaMu_q = 0.0;
      for(unsigned int i=0;i<atomCount;i++){
	// point self interactions
	q += realTopo->atoms[i].scaledCharge*realTopo->atoms[i].scaledCharge;
	DeltaMu_q += realTopo->atoms[i].scaledCharge * realTopo->atoms[i].deltaQ;
      }  // end for loop
      
      myPointSelfEnergy = -q*myAlpha/sqrt(M_PI);
      myPointSelfDeltaMu = -2 * DeltaMu_q*myAlpha/sqrt(M_PI);
    }
  
    //
    // Charged system energy
    //
    myChargedSystemEnergy   = 0.0;
    myChargedSystemDeltaMu  = 0.0;
    if(correction) {
      Real q = 0.0;
      Real DeltaMu_q = 0.0;
      for(unsigned int i=0;i<atomCount;i++) {
	q += realTopo->atoms[i].scaledCharge;
	DeltaMu_q += realTopo->atoms[i].deltaQ;
      }
      if(fabs(q * 0.00268283) > 1.0e-5) {
	myChargedSystemEnergy = -M_PI/(2.0*myV*myAlphaSquared)*q*q;
	myChargedSystemDeltaMu = -M_PI/(myV*myAlphaSquared)*q*DeltaMu_q;
      }
    }

    myLattice = boundaryConditions.buildLatticeVectors(myRc);
    myLattice.insert(myLattice.begin(),Vector3D(0,0,0));

    report << hint <<"Ewald";
#if defined(USE_EWALD_EXACT_ERF) && defined(USE_EWALD_NO_SINCOS_TABLE)
    report << "(Exact)";
#else
#ifdef USE_EWALD_EXACT_ERF
    report << "(Exact erf)";
#endif
#ifdef USE_EWALD_NO_SINCOS_TABLE
    report << "(No sincos table)";
#endif
#endif  
    report <<": alpha="<<myAlpha<<", V="<<myV<<", Rc="
	   <<myRc<<", Kc ("<<myK.size()<<")="<<myKc<<", n="<<myLattice.size()<<", accuracy="
	   <<myAccuracy<<", misses="
	   <<(Real)100*misses/(Real)(myK.size()> 0 ? myK.size():1) <<"%."<<endr;

    // Now we have all pre-computed stuff, non-periodic case may need an update when 
    // the simulation box expands or shrinks to much (+-10%)
    if(TBoundaryConditions::PERIODIC)
      myCached = true;

    myV = boundaryConditions.getVolume();

    if(dontHint)
      Report::report << Report::dohint;
  }

  //________________________________________________________________realTerm function  
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::realTerm(const RealTopologyType* realTopo,
								      const Vector3DBlock* positions, 
								      Vector3DBlock* forces, 
								      ScalarStructure* energies,
								      Real& realEnergy,
								      Real& realDeltaMu,
								      unsigned int n) {
#if defined(DEBUG_EWALD_TIMING)
    myReal.start();
#endif

    // Real-space term
    CellPairType thisPair;
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();
    unsigned int count = 0;
    for (; !enumerator.done(); enumerator.next()) {
      enumerator.get(thisPair);
      bool notSameCell = enumerator.notSameCell();
    
      if(!notSameCell){
	count++;
	if(count > n)
	  break;
      }
      
      for(int i=thisPair.first; i!=-1; i=realTopo->atoms[i].cellListNext){
	Real qi  = realTopo->atoms[i].scaledCharge;
	Vector3D ri((*positions)[i]), fi;
        int mi = realTopo->atoms[i].molecule;
	for(int j=(notSameCell ? thisPair.second:i); j!=-1; j=realTopo->atoms[j].cellListNext){

          Vector3D rijMinimal(boundaryConditions.minimalDifference(ri,(*positions)[j]));
          int mj = realTopo->atoms[j].molecule;	
          bool same = (mi==mj);  
          ExclusionClass excl = (same?realTopo->exclusions.check(i,j):EXCLUSION_NONE);
	  if(i == j)
	    excl = EXCLUSION_FULL;
	  for(unsigned int k=0;k<myLattice.size();k++){
	    // Check for an exclusion.
	    if (excl == EXCLUSION_FULL){ 
	      excl = EXCLUSION_NONE;
	      continue;
	    }

	    Vector3D rij(rijMinimal+myLattice[k]);
	    Real rSquared = rij.normSquared();
	    
	    // Do switching function rough test.
	    if (!switchingFunction.roughTest(rSquared))
	      continue;
	    
	    // we must determine which type of pair interaction this is:
	    // type 0 = untransformed-untransformed
	    // type 1 = untransformed-transformed
	    // type 2 = transformed-transformed (intramolecular)
	    bool atom1_scaled, atom2_scaled;
	    atom1_scaled = atom2_scaled = false;
	    int myStage, /* OldStage,*/ atom_stage;
            myStage /*= OldStage */ = atom_stage = 0;
                      
	    // Use the molecule types to determine if the atoms are being transformed or not
	    if (realTopo->molecules[mi].lambda != 0.0) {atom1_scaled = true;}
	    if (realTopo->molecules[mj].lambda != 0.0) {atom2_scaled = true;}
	    
	    // integer used to select which interaction type to compute   
	    int Choice;
	    
	    // determine the value of Choice
            // neither atom is being transformed
            if ( !(atom1_scaled) && !(atom2_scaled) ) {Choice = 0;}
            // intramolecular interaction on a transforming molecule
            else if ( (atom1_scaled) && (atom2_scaled) ) {             

              // determine the current transformation stage of the molecule
              myStage = static_cast<int>(floor(realTopo->molecules[mi].lambda)) + 1;

              // get the transformation stage #'s for these atoms
              int atom1_stage = realTopo->atoms[i].stageNumber;
              int atom2_stage = realTopo->atoms[j].stageNumber;

              // if the current stage is not the same as either atom1_stage or atom2_stage
              // the neither atom is currently being transformed
              if (myStage != atom1_stage && myStage != atom2_stage) Choice = 0;
              else Choice = 2;
            }
            // one of the atoms is being transformed
            else {
              // determine the current transformation stage of the molecule
              // and get the transformation stage # for this atoms
              if (atom1_scaled) {
                myStage = static_cast<int>(floor(realTopo->molecules[mi].lambda)) + 1;
                atom_stage = realTopo->atoms[i].stageNumber;}
              else {
                myStage = static_cast<int>(floor(realTopo->molecules[mj].lambda)) + 1;
                atom_stage = realTopo->atoms[j].stageNumber;}
        
              // if the current stage is not the same as atom_stage
              // then neither atom is currently being transformed
              if (myStage != atom_stage) Choice = 0;
              else Choice = 1;
            }
	    
	    // Compute the appropriate electrostatic interaction
	    Real deltaMu = 0.0;
	    Real energy = 0.0;
	    Real force = 0.0;
	    Real qq, rr, ar, e, q_Dq, dmu, myLambda;
	    Real qq_correction_add, qq_correction_minus, q_Dq_correction_add, q_Dq_correction_minus;
	    Real e_correction, myCorrection, myDMUcorrection;
	    qq = rr = ar = e = q_Dq = dmu = myLambda = 0.0;
	    qq_correction_add = qq_correction_minus = q_Dq_correction_add = q_Dq_correction_minus = 0.0;
	    e_correction = myCorrection = myDMUcorrection = 0.0;
#ifdef USE_EWALD_EXACT_ERF
	    Real a = 0.0;
#endif

	    // get the charge on atom j
	    Real qj = realTopo->atoms[j].scaledCharge;
	    
	    // compute the separation distance
	    Real r = sqrt(rSquared);
		  
	    switch (Choice) {
	      //----------------------------------------------------------
	    case 0:
	      // interaction between two non-transforming atoms
	      // compute the raw interaction energy and force

	      // skip if either of the charges is zero
	      if (qi == 0.0 || qj == 0.0) break;

	      qq = qi*qj;
	      
	      if (excl == EXCLUSION_MODIFIED) 
		qq *= realTopo->coulombScalingFactor;
	      
	      // Approximation Abramowitz & Stegun p299.
	      // Energy
#ifndef USE_EWALD_EXACT_ERF
	      rr = 1.0/r;
	      ar = myAlpha*r;
	      e = qq*exp(-ar*ar);
	      energy = poly5(ar)*e*rr;
	      force = ((energy+my2AlphaPI*e)*rr*rr);
#else	  
	      a = erfc(myAlpha*r)/r;
	      energy = qq*a;
	      force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	      break;
	      //----------------------------------------------------------
	    case 1:  
	      // interaction between non-transforming atom and a transforming atom
	      // compute the raw interaction energy, force, and chemical potential difference

	      // skip if either of the charges is zero
	      if (qi == 0.0 || qj == 0.0) break;

	      qq = qi*qj;
	      q_Dq = (qi * realTopo->atoms[j].deltaQ
		      + qj * realTopo->atoms[i].deltaQ);
	      
	      // Approximation Abramowitz & Stegun p299.
	      // Energy
#ifndef USE_EWALD_EXACT_ERF
	      rr = 1.0/r;
	      ar = myAlpha*r;
	      e = qq*exp(-ar*ar);
	      dmu = q_Dq*exp(-ar*ar);
	      energy = poly5(ar)*e*rr;
	      deltaMu = poly5(ar)*dmu*rr;
	      force = ((energy+my2AlphaPI*e)*rr*rr);
#else	  
	      a = erfc(myAlpha*r)/r;
	      energy = qq*a;
	      deltaMu = q_Dq*a;
	      force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	      break;
	      //----------------------------------------------------------
	    case 2:  
	      // intramolecular interaction on a transforming molecule
	      // get the value of lambda for this molecule
	      myLambda = realTopo->molecules[mi].lambda;

              // determine the current transformation stage of the molecule
              int OldStage = myStage - 1;
              //if (myLambda > myStage) myLambda = myStage;
              //else if (myLambda < OldStage) myLambda = OldStage;
        	      
	      // Energy and chemical potential difference
	      // The interaction energy is the sum of the interaction in the
	      // old state scaled by (1 - lambda) plus the interaction in the new state
	      // scaled by lambda 
	      r = sqrt(rSquared);
	      qj = realTopo->atoms[j].scaledCharge;
	      qq = (myStage - myLambda) * realTopo->atoms[i].Qold * realTopo->atoms[j].Qold
		+ (myLambda - OldStage) * realTopo->atoms[i].Qnew * realTopo->atoms[j].Qnew;
	      
	      // The chemical potential difference is the derivative of the 
	      // energy above (qq) with respect to lambda
	      q_Dq = realTopo->atoms[i].Qnew * realTopo->atoms[j].Qnew
		- realTopo->atoms[i].Qold * realTopo->atoms[j].Qold;

	      // correction terms for intramolecular self energy and chemical potential difference
	      // these correction terms are needed because the reciprocal space term contains an
	      // intramolecular self term that is proportional to atoms[i].scaledCharge * atoms[j].scaledCharge,
	      // which is NOT the correct intramolecular self interaction (see qq and q_Dq terms above).
	      // Correcting the reciprocal intramolecular interaction amounts to the addition of a factor
	      // of qq to the energy and q_Dq to the chemical potential difference, and subtraction of a
	      // factor of atoms[i].scaledCharge * atoms[j].scaledCharge from the energy and
	      // subtraction of (atoms[i].scaledCharge * atoms[j].deltaQ + atoms[j].scaledCharge * atoms[i].deltaQ)
	      // from the chemical potential difference
	      qq_correction_add = qq;
	      q_Dq_correction_add = q_Dq;
	      qq_correction_minus = 0.0;
	      q_Dq_correction_minus = 0.0;
	      
	      if (excl == EXCLUSION_MODIFIED) {
		qq *= realTopo->coulombScalingFactor;
		q_Dq *= realTopo->coulombScalingFactor;
		qq_correction_add *= realTopo->coulombScalingFactor;
		q_Dq_correction_add *= realTopo->coulombScalingFactor;
	      }
	      else {
		// this is an unmodified intramolecular interaction, so subtract the reciprocal space
		// energy and chemical potential difference
		qq_correction_minus = qi * realTopo->atoms[j].scaledCharge;
		q_Dq_correction_minus = qi * realTopo->atoms[j].deltaQ + qj * realTopo->atoms[i].deltaQ;
	      }
      
	      // Approximation Abramowitz & Stegun p299.
	      // Energy
#ifndef USE_EWALD_EXACT_ERF
	      rr = 1.0/r;
	      ar = myAlpha*r;
	      e = qq*exp(-ar*ar);
	      e_correction = erf(myAlpha*r)*rr;
	      dmu = q_Dq*exp(-ar*ar);
	      energy = poly5(ar)*e*rr;
	      deltaMu = poly5(ar)*dmu*rr;
	      force = ((energy+my2AlphaPI*e)*rr*rr);
#else	  
	      a = erfc(myAlpha*r)/r;
	      e_correction = erf(myAlpha*r)/r;
	      energy = qq*a;
	      deltaMu = q_Dq*a;
	      force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	      
	      // add in the reciprocal space correction terms
	      myCorrection = e_correction * (qq_correction_add - qq_correction_minus);
	      myDMUcorrection = e_correction * (q_Dq_correction_add - q_Dq_correction_minus);
	      energy += myCorrection;
	      deltaMu += myDMUcorrection;
	      break;
	      //----------------------------------------------------------
	    } // end switch structure
	    
	    
	    // Calculate the switched force, energy, and delta_mu.
	    Real switchingValue, switchingDeriv;
	    switchingFunction(switchingValue, switchingDeriv, rSquared);
	    energy = energy * switchingValue;
	    deltaMu = deltaMu * switchingValue;
	    // This has a - sign because the force is the negative of the 
	    // derivative of the energy (divided by the distance between the atoms).
	    force = force * switchingValue - energy * switchingDeriv;
	    // Force, F_ij
	    realEnergy += energy; 
	    realDeltaMu += deltaMu;
	    Vector3D fij(rij*force);
	    fi -= fij;
	    (*forces)[j] += fij;
	    
	    // compute the vector between molecular centers of mass
	    if(!same && doMolVirial){
	      // Add to the atomic and molecular virials
	      energies->addVirial(fij,rij,realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position, realTopo->molecules[mj].position));
	    }
	    else if(doVirial) {
	      energies->addVirial(fij,rij);
	    }
	    excl = EXCLUSION_NONE;
	  }
	}
        (*forces)[i] += fi;
      }
    }
#if defined(DEBUG_EWALD_TIMING)
    myReal.stop();
#endif
  }

  //________________________________________________________________reciprocalTerm function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::reciprocalTerm(const RealTopologyType* realTopo,
									    const Vector3DBlock* positions, 
									    Vector3DBlock* forces, 
									    ScalarStructure* energies,
									    Real& reciprocalEnergy,
									    Real& reciprocalDeltaMu,
									    unsigned int from, unsigned int to) {
  
#if defined(DEBUG_EWALD_TIMING)
    myReciprocal.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    
    Real energy = 0.0;
    Real deltaMu = 0.0;
#ifndef USE_EWALD_NO_SINCOS_TABLE
    // Precompute/ cache cos/ sin (r*N*2*PI/L) for the lattice vectors in each dimension 
    for(unsigned int j=0;j<atomCount;j++){
      int l = 2*j*myHKLmax;
      Vector3D r = boundaryConditions.minimalPosition((*positions)[j]);
      // Multiply charge only with x-coord of each particle
      // since we use the add theorem 

      Real qi = realTopo->atoms[j].scaledCharge;
      Real Dqi = realTopo->atoms[j].deltaQ;

      Real x = r.x*2.0*M_PI*myLXr;
      Real y = r.y*2.0*M_PI*myLYr;
      Real z = r.z*2.0*M_PI*myLZr;
      Real xsin = sin(x);
      Real ysin = sin(y);
      Real zsin = sin(z);
      Real xcos = cos(x);
      Real ycos = cos(y);
      Real zcos = cos(z);
      
      // The first two cos/ sin values
      // sin(r*0*2*PI/L)
      mySinCos[l  ].x = 0.0;
      mySinCos[l  ].y = 0.0;
      mySinCos[l  ].z = 0.0;
      myDMU_SinCos[l  ].x = 0.0;
      myDMU_SinCos[l  ].y = 0.0;
      myDMU_SinCos[l  ].z = 0.0;
            
      // cos(r*0*2*PI/L)
      mySinCos[l+1].x = qi*1.0;
      mySinCos[l+1].y = 1.0;
      mySinCos[l+1].z = 1.0;
      myDMU_SinCos[l+1].x = Dqi*1.0;
      myDMU_SinCos[l+1].y = 1.0;
      myDMU_SinCos[l+1].z = 1.0;
      
      // sin(r*1*2*PI/L)
      mySinCos[l+2].x = qi*xsin;
      mySinCos[l+2].y = ysin;
      mySinCos[l+2].z = zsin;
      myDMU_SinCos[l+2].x = Dqi*xsin;
      myDMU_SinCos[l+2].y = ysin;
      myDMU_SinCos[l+2].z = zsin;
      
      // cos(r*1*2*PI/L)
      mySinCos[l+3].x = qi*xcos;
      mySinCos[l+3].y = ycos;
      mySinCos[l+3].z = zcos;
      myDMU_SinCos[l+3].x = Dqi*xcos;
      myDMU_SinCos[l+3].y = ycos;
      myDMU_SinCos[l+3].z = zcos;
      
      // Using add theorem to compute sin(r*2*2*PI/L to r*myHKLmax*2*PI/L) 
      // and cos(r*2*2*PI/L to r*myHKLmax*2*PI/L)
      for(int i=4;i<2*myHKLmax;i+=2){
	mySinCos[l+i  ].x = xsin*mySinCos[l+i-1].x+xcos*mySinCos[l+i-2].x;
	mySinCos[l+i  ].y = ysin*mySinCos[l+i-1].y+ycos*mySinCos[l+i-2].y;
	mySinCos[l+i  ].z = zsin*mySinCos[l+i-1].z+zcos*mySinCos[l+i-2].z;
	mySinCos[l+i+1].x = xcos*mySinCos[l+i-1].x-xsin*mySinCos[l+i-2].x;
	mySinCos[l+i+1].y = ycos*mySinCos[l+i-1].y-ysin*mySinCos[l+i-2].y;
	mySinCos[l+i+1].z = zcos*mySinCos[l+i-1].z-zsin*mySinCos[l+i-2].z;
	
	myDMU_SinCos[l+i  ].x = xsin*myDMU_SinCos[l+i-1].x+xcos*myDMU_SinCos[l+i-2].x;
	myDMU_SinCos[l+i  ].y = ysin*myDMU_SinCos[l+i-1].y+ycos*myDMU_SinCos[l+i-2].y;
	myDMU_SinCos[l+i  ].z = zsin*myDMU_SinCos[l+i-1].z+zcos*myDMU_SinCos[l+i-2].z;
	myDMU_SinCos[l+i+1].x = xcos*myDMU_SinCos[l+i-1].x-xsin*myDMU_SinCos[l+i-2].x;
	myDMU_SinCos[l+i+1].y = ycos*myDMU_SinCos[l+i-1].y-ysin*myDMU_SinCos[l+i-2].y;
	myDMU_SinCos[l+i+1].z = zcos*myDMU_SinCos[l+i-1].z-zsin*myDMU_SinCos[l+i-2].z;
      }
    }
    
    //int nk = myK.size();
    int lastH = Constant::MAX_INT;
    int lastK = Constant::MAX_INT;
#endif
    // atomic virial
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;

    // molecular virial
    bool doMolVirial = energies->molecularVirial();
    bool doVirial = energies->virial();

    for(unsigned int l=from;l<to;l++){
      // Energy
      Vector3D k = myK[l];
      Real sumSin = 0.0;
      Real sumCos = 0.0;
      Real DMU_sumSin = 0.0;
      Real DMU_sumCos = 0.0;
      Real kSquared = myKSquared[l];

#ifndef USE_EWALD_NO_SINCOS_TABLE
      int indexH = myKInt[l].h;
      int indexK = myKInt[l].k;
      int indexKabs = abs(indexK);
      int indexKsign = 1;
      if(indexK < 0)
        indexKsign = -1;
      int indexL = myKInt[l].l;
      int indexLabs = abs(indexL);
      int indexLsign = 1;
      if(indexL < 0)
        indexLsign = -1;
    
      // Precompute and cache sin/ cos for h and k
      // using the precompute table of sin/ cos
      // Hit/miss rate vary from 6:1 to 12:1
      if(indexH != lastH || indexK != lastK){
	for(unsigned int i=0;i<atomCount;i++){
	  Real xsin =            mySinCos[i*myHKLmax*2+indexH*2  ].x;
	  Real xcos =            mySinCos[i*myHKLmax*2+indexH*2+1].x;
	  Real ysin = indexKsign*mySinCos[i*myHKLmax*2+indexKabs*2  ].y;
	  Real ycos =            mySinCos[i*myHKLmax*2+indexKabs*2+1].y;
	  
	  myLastSinCos[i*2  ] = xsin*ycos + xcos*ysin;
	  myLastSinCos[i*2+1] = xcos*ycos - xsin*ysin;
	  
	  Real DMU_xsin =            myDMU_SinCos[i*myHKLmax*2+indexH*2  ].x;
	  Real DMU_xcos =            myDMU_SinCos[i*myHKLmax*2+indexH*2+1].x;
	  Real DMU_ysin = indexKsign*myDMU_SinCos[i*myHKLmax*2+indexKabs*2  ].y;
	  Real DMU_ycos =            myDMU_SinCos[i*myHKLmax*2+indexKabs*2+1].y;
	  
	  myDMU_LastSinCos[i*2  ] = DMU_xsin*DMU_ycos + DMU_xcos*DMU_ysin;
	  myDMU_LastSinCos[i*2+1] = DMU_xcos*DMU_ycos - DMU_xsin*DMU_ysin;
	}
	lastH = indexH;
	lastK = indexK;
      }
#endif    
      for(unsigned int i=0;i<atomCount;i++){
#ifdef USE_EWALD_NO_SINCOS_TABLE
	Real qi = realTopo->atoms[i].scaledCharge;
	Real Dqi = realTopo->atoms[i].deltaQ;
	// It does not matter if coordinates are not in the minimal image since
	// they are multiplied by 2PI/l, which is a shift of 2PI of a. 
	Real a = k.dot(boundaryConditions.minimalPosition((*positions)[i]));
	//Real a = k.dot((*positions)[i]);
	Real sinA = qi*sin(a);
	Real cosA = qi*cos(a);
	Real DMU_sinA = Dqi*sin(a);
	Real DMU_cosA = Dqi*cos(a);
#else
	
	Real xysin = myLastSinCos[i*2  ];
	Real xycos = myLastSinCos[i*2+1];
	
	Real zsin = indexLsign*mySinCos[i*myHKLmax*2+indexLabs*2  ].z;
	Real zcos = mySinCos[i*myHKLmax*2+indexLabs*2+1].z;
	
	Real sinA = xysin*zcos + xycos*zsin;
	Real cosA = xycos*zcos - xysin*zsin;
	
	Real DMU_xysin = myDMU_LastSinCos[i*2  ];
	Real DMU_xycos = myDMU_LastSinCos[i*2+1];
	
	Real DMU_zsin = indexLsign*myDMU_SinCos[i*myHKLmax*2+indexLabs*2  ].z;
	Real DMU_zcos = myDMU_SinCos[i*myHKLmax*2+indexLabs*2+1].z;
	
	Real DMU_sinA = DMU_xysin*DMU_zcos + DMU_xycos*DMU_zsin;
	Real DMU_cosA = DMU_xycos*DMU_zcos - DMU_xysin*DMU_zsin;
#endif
	mySinCosA[i*2  ] = sinA;
	mySinCosA[i*2+1] = cosA;
	sumSin += sinA;
	sumCos += cosA;
	
	myDMU_SinCosA[i*2  ] = DMU_sinA;
	myDMU_SinCosA[i*2+1] = DMU_cosA;
	DMU_sumSin += DMU_sinA;
	DMU_sumCos += DMU_cosA;
      }

      // Energy
      Real b = 1.0/kSquared*exp(-kSquared*myAlphaSquaredr/4.0);
      Real e = b*(sumSin*sumSin+sumCos*sumCos);
      energy += e;
      
      // Chemical potential difference
      Real dMu = 2*b*( DMU_sumSin*sumSin + DMU_sumCos*sumCos );
      deltaMu += dMu;
      
      // Virial
      if(doVirial){
	Real c = 2.0*(1.0/kSquared+myFac);
	virialxx += e * (1.0-c*k.x * k.x);
	virialxy -= e * c * k.x * k.y;
	virialxz -= e * c * k.x * k.z;
	virialyy += e * (1.0-c*k.y * k.y);
	virialyz -= e * c * k.y * k.z;
	virialzz += e * (1.0-c*k.z * k.z);
      }
      
      // Force, F_i 
      Real a = 8.0*M_PI*myVr*b;
      for(unsigned int i=0;i<atomCount;i++){

	// compute the force on atom i from the reciprocal space part
	Vector3D fi(k*(a*(mySinCosA[i*2  ]*sumCos - mySinCosA[i*2+1]*sumSin)));
	(*forces)[i] += fi;     
	
	// compute the reciprocal space contribution to the molecular virial
	// this expression is taken from Alejandre, Tildesley, and Chapela, J. Chem. Phys. 102 (11), 4574.
	if(doMolVirial){
	  // get the ID# of the molecule to which this atom belongs
	  int Mi = realTopo->atoms[i].molecule;
	  
	  // compute the vector from atom i to the center of mass of the molecule
	  Vector3D ria(boundaryConditions.minimalPosition((*positions)[i]));
	  Vector3D mri(realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[Mi].position));

	  energies->addMolVirial(fi,mri);
	}
      }
    }
    
    Real c = 4.0*M_PI*myVr;
    reciprocalEnergy += c*energy;
    reciprocalDeltaMu += c*deltaMu;
    
    // atomic virial
    (*energies)[ScalarStructure::VIRIALXX] += c*virialxx;
    (*energies)[ScalarStructure::VIRIALXY] += c*virialxy;
    (*energies)[ScalarStructure::VIRIALXZ] += c*virialxz;
    (*energies)[ScalarStructure::VIRIALYX] += c*virialxy;
    (*energies)[ScalarStructure::VIRIALYY] += c*virialyy;
    (*energies)[ScalarStructure::VIRIALYZ] += c*virialyz;
    (*energies)[ScalarStructure::VIRIALZX] += c*virialxz;
    (*energies)[ScalarStructure::VIRIALZY] += c*virialyz;
    (*energies)[ScalarStructure::VIRIALZZ] += c*virialzz;
  
    // molecular virial
    (*energies)[ScalarStructure::MOLVIRIALXX] += c*virialxx;
    (*energies)[ScalarStructure::MOLVIRIALXY] += c*virialxy;
    (*energies)[ScalarStructure::MOLVIRIALXZ] += c*virialxz;
    (*energies)[ScalarStructure::MOLVIRIALYX] += c*virialxy;
    (*energies)[ScalarStructure::MOLVIRIALYY] += c*virialyy;
    (*energies)[ScalarStructure::MOLVIRIALYZ] += c*virialyz;
    (*energies)[ScalarStructure::MOLVIRIALZX] += c*virialxz;
    (*energies)[ScalarStructure::MOLVIRIALZY] += c*virialyz;
    (*energies)[ScalarStructure::MOLVIRIALZZ] += c*virialzz;

#if defined(DEBUG_EWALD_TIMING)
    myReciprocal.stop();
#endif
  }

  //________________________________________________________________correctionTerm function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::correctionTerm(const RealTopologyType* realTopo,
									    const Vector3DBlock* positions, 
									    Vector3DBlock* forces, 
									    ScalarStructure* energies,
									    Real& intraMolecularEnergy,
									    Real& intraMolecularDeltaMu,
									    unsigned int from, unsigned int to){
#if defined(DEBUG_EWALD_TIMING)
    myIntra.start();
#endif
    // Intra-molecular term  
    bool doVirial = energies->virial();  
    const std::vector<ExclusionPair>& exclusions = realTopo->exclusions.getTable();
    for(unsigned int i=from;i<to;i++){
      ExclusionPair excl = exclusions[i];
      Real rSquared;
      Vector3D rij(realTopo->boundaryConditions.minimalDifference((*positions)[excl.a1],(*positions)[excl.a2],rSquared));

      // get the value of lambda for this molecule
      int M1 = realTopo->atoms[excl.a1].molecule;
      Real myLambda = realTopo->molecules[M1].lambda;

      // intramolecular self energy and chemical potential difference contributions
      // If this is a transforming molecule (lambda /= 0) then we must completely subtract the
      // excluded interaction.  If this is not a transforming molecule, then we apply the 
      // 1-4 scaling factor (assuming one has been specified) as normal.
      Real qq = realTopo->atoms[excl.a1].scaledCharge*realTopo->atoms[excl.a2].scaledCharge;
      Real q_Dq = realTopo->atoms[excl.a1].deltaQ * realTopo->atoms[excl.a2].scaledCharge
        + realTopo->atoms[excl.a1].scaledCharge * realTopo->atoms[excl.a2].deltaQ;

      // check to see if this is not a transforming molecule
      // determine the current transformation stage of the molecule
      int myStage = static_cast<int>(floor(myLambda)) + 1;

      // get the transformation stage #'s for these atoms
      int atom1_stage = realTopo->atoms[excl.a1].stageNumber;
      int atom2_stage = realTopo->atoms[excl.a2].stageNumber;

      // if the current stage is not the same as either atom1_stage or atom2_stage, or lambda = 0,
      // then neither atom is currently being transformed       
      if ( (myLambda == 0.0 || (myStage != atom1_stage && myStage != atom2_stage))
           && excl.excl == EXCLUSION_MODIFIED) {
        qq *= 1-realTopo->coulombScalingFactor;
        q_Dq *= 1-realTopo->coulombScalingFactor;
      }
      
      Real r  = sqrt(rSquared);
      Real rr = 1/r;
      Real e = erf(myAlpha*r)*rr;
      // Intra-molecular selv energy
      intraMolecularEnergy -= qq*e;
      intraMolecularDeltaMu -= q_Dq*e;
      // Intra-molecular selv force
      Vector3D fij(rij*(qq*(my2AlphaPI*exp(-myAlphaSquared*rSquared)-e)*rr*rr));
      (*forces)[excl.a1] -= fij;
      (*forces)[excl.a2] += fij;
      if(doVirial)	
        energies->addVirial(fij,rij);
    }
#if defined(DEBUG_EWALD_TIMING)
    myIntra.stop();
#endif
  }

  // Surface diplole term
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::surfaceDipoleTerm(const RealTopologyType* realTopo,
									       const Vector3DBlock* positions, 
									       Vector3DBlock* forces, 
									       ScalarStructure* energies,
									       Real& surfaceDipoleEnergy){
#if defined(DEBUG_EWALD_TIMING)
    mySurface.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    Vector3D sum(0,0,0);
    for(unsigned int i=0;i<atomCount;i++)
      sum += boundaryConditions.minimalPosition((*positions)[i])*realTopo->atoms[i].scaledCharge;
    // Energy
    surfaceDipoleEnergy = 2.0/3.0*M_PI*myVr*sum.normSquared();
  
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;
    // Force, F_i and virial_i (not confirmed)
    sum *= 2.0/3.0*M_PI*myVr;
    for(unsigned int i=0;i<atomCount;i++){
      Vector3D force(sum*realTopo->atoms[i].scaledCharge);
      Vector3D ri(boundaryConditions.minimalPosition((*positions)[i]));
      (*forces)[i] += force;
      virialxx += force.x*ri.x;
      virialxy += force.x*ri.y;
      virialxz += force.x*ri.z;
      virialyy += force.y*ri.y;
      virialyz += force.y*ri.z;
      virialzz += force.z*ri.z;
    }
    (*energies)[ScalarStructure::VIRIALXX] += virialxx;
    (*energies)[ScalarStructure::VIRIALXY] += virialxy;
    (*energies)[ScalarStructure::VIRIALXZ] += virialxz;
    (*energies)[ScalarStructure::VIRIALYX] += virialxy;
    (*energies)[ScalarStructure::VIRIALYY] += virialyy;
    (*energies)[ScalarStructure::VIRIALYZ] += virialyz;
    (*energies)[ScalarStructure::VIRIALZX] += virialxz;
    (*energies)[ScalarStructure::VIRIALZY] += virialyz;
    (*energies)[ScalarStructure::VIRIALZZ] += virialzz;
  
#if defined(DEBUG_EWALD_TIMING)
    mySurface.stop();
#endif
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  std::string iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					       TCellManager,
					       real,
					       reciprocal,
					       correction,
					       TSwitchingFunction>::getIdNoAlias() const{
    return (iSGCoulombForce::keyword + " -algorithm "+ keyword + 
	    std::string((real)       ? std::string(" -real")       : std::string("")) + 
	    std::string((reciprocal) ? std::string(" -reciprocal") : std::string("")) +
	    std::string((correction) ? std::string(" -correction") : std::string("")) +	  
	    std::string((TSwitchingFunction::getId() != CutoffSwitchingFunction::getId()) ? 
			std::string(std::string(" -switchingFunction " + TSwitchingFunction::getId())) : std::string("")));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  void iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					TCellManager,
					real,
					reciprocal,
					correction,
					TSwitchingFunction>::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-alpha",Value(myAlpha),-1.0));    
    parameters.push_back(Parameter("-accuracy",Value(myAccuracy,ConstraintValueType::Positive()),0.00001));
    if(TBoundaryConditions::VACUUM)
      parameters.push_back(Parameter("-j",Value(myExpansionFactor,ConstraintValueType::Positive()),3.0));
  }

  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TSwitchingFunction>
  Force*  iSGNonbondedFullEwaldSystemForce<TBoundaryConditions,
					   TCellManager,
					   real,
					   reciprocal,
					   correction,
					   TSwitchingFunction>::doMake(std::string& errMsg, std::vector<Value> values) const{
    Real alpha           = values[0];
    Real accuracy        = values[1];
    Real expansionFactor = (TBoundaryConditions::VACUUM?(Real)values[2]:3.0);
    std::string err      = "";

    if(!values[0].valid())
      err +=" alpha \'"+values[0].getString()+"\' not valid.";

    if(!values[1].valid())
      err +=" accuracy \'"+values[1].getString()+"\' not valid.";

    if(TBoundaryConditions::VACUUM && !values[2].valid())
      err +=" expansionFactor \'"+values[2].getString()+"\' not valid.";
    else if(expansionFactor <= 1.0)
      err += keyword + " simulation box expansion factor (="+toString(expansionFactor)+") > 1.0.";

    if(!err.empty()){
      errMsg += " force "+keyword+" :"+err;
      return NULL;
    }

    return (new iSGNonbondedFullEwaldSystemForce(alpha, accuracy,expansionFactor));
  }
}
#endif /* ISGNONBONDEDFULLEWALDSYSTEMFORCE_H */
