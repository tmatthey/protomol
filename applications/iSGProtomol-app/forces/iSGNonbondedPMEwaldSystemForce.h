/* -*- c++ -*- */
#ifndef ISGNONBONDEDPMEWALDSYSTEMFORCE_H
#define ISGNONBONDEDPMEWALDSYSTEMFORCE_H

#include "SystemForce.h"
#include "iSGCoulombForce.h"
#include "CutoffSwitchingFunction.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "PeriodicBoundaryConditions.h"
#include "iSGCoulombForce.h"
#include "mathutilities.h"
#include "iSGNonbondedPMEwaldSystemForceBase.h"
#include "iSGGrid.h"
#include "Timer.h"
#include "CutoffSwitchingFunction.h"

using namespace ProtoMol::Report;

//#define DEBUG_PME_TIMING
//#define DEBUG_PME_ENERGIES
//#define USE_PME_EXACT_ERF
namespace ProtoMol {
  //_________________________________________________________________ ISGNonbondedPMEwaldSystemForce
  
  template<class TBoundaryConditions, 
	   class TCellManager,
	   bool  real,
	   bool  reciprocal,
	   bool  correction,
	   class TInterpolation,
	   class TSwitchingFunction>
  class iSGNonbondedPMEwaldSystemForce: public SystemForce, private iSGNonbondedPMEwaldSystemForceBase {
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
    iSGNonbondedPMEwaldSystemForce();
    iSGNonbondedPMEwaldSystemForce(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int order, Real cutoff, Real accuracy, Real alpha, Real expansionFactor);

    virtual ~iSGNonbondedPMEwaldSystemForce();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ISGNonbondedPMEwaldSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

  private:
    void initialize(const RealTopologyType* realTopo,
		    const Vector3DBlock* positions);
    
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
			Real& reciprocalDeltaMu);
    
    void reciprocalTermParallel(const RealTopologyType* realTopo,
				const Vector3DBlock* positions, 
				Vector3DBlock* forces, 
				ScalarStructure* energies,
				Real& reciprocalEnergy,
				Real& reciprocalDeltaMu);
    
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
    virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual std::string getIdNoAlias() const;
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return (TBoundaryConditions::PERIODIC ? 7:8);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
    bool myCached;
    Real myExpansionFactor;
    Real myTRate;             //
    Real myAccuracy;

    Real myAlpha;             // 
    bool myAlphaDefault;

    Real myAlphaSquared;
    Real myAlphaSquaredr;
    Real my2AlphaPI;

    Real myRc;                // Cutoff real term
    Real myRcSquared;         // Cutoff squared real term
    Real myKc;                // Cutoff reciprocal term

    Real myLX, myLY, myLZ;    // cell basis vectors
    Real myLXr,myLYr,myLZr;   // reciprocal cell basis vectors
    Real myV,myVr;            // Volume and inverse volume
    Vector3D myOrigin;        // xyz coordinates of the origin
    unsigned int myNX;        // # of grid points in x-direction
    unsigned int myNY;        // # of grid points in y-direction
    unsigned int myNZ;        // # of grid points in z-direction

    iSGGrid<TInterpolation> myISGGrid;   // object which computes the FFT's

    unsigned int myInterOrder; // order of the B-spline interpolation

    Real myPointSelfEnergy;    // Precomputed energy terms
    Real myPointSelfDeltaMu;
    Real myChargedSystemEnergy;
    Real myChargedSystemDeltaMu;

    PeriodicBoundaryConditions boundaryConditions;

    TSwitchingFunction switchingFunction;
    EnumeratorType enumerator;
    std::vector<Vector3D> myLattice;
#if defined(DEBUG_PME_TIMING)
    Timer myReal;
    Timer myReciprocal;
    Timer myIntra;
    Timer mySurface;
#endif
  };

  //______________________________________________________________________ INLINES
  
  //_______________________________________________________________________ empty constructor
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,
				 TCellManager,
				 real,
				 reciprocal,
				 correction,
				 TInterpolation,
				 TSwitchingFunction>::iSGNonbondedPMEwaldSystemForce() : SystemForce(),myCached(false),myExpansionFactor(3.0),myTRate(0.0),myAccuracy(0.0),myAlpha(-1.0),myAlphaDefault(true),myRc(0.0),myKc(0.0),myV(-1.0),myNX(0),myNY(0),myNZ(0),myInterOrder(0) {
#if defined(DEBUG_PME_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }

  
  //_______________________________________________________________________ constructor
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,
				 TCellManager,
				 real,
				 reciprocal,
				 correction,
				 TInterpolation,
				 TSwitchingFunction>::iSGNonbondedPMEwaldSystemForce(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int order, Real cutoff, Real accuracy, Real alpha, Real expansionFactor) : SystemForce(),myCached(false),myExpansionFactor(expansionFactor),myTRate(0.0),myAccuracy(accuracy),myAlpha(alpha),myAlphaDefault(alpha<= 0.0),myRc(cutoff),myKc(0.0),myV(-1.0),myNX(nx),myNY(ny),myNZ(nz),myInterOrder(order){
    switchingFunction=TSwitchingFunction(myRc);
#if defined(DEBUG_PME_TIMING)
    myReal.reset();
    myReciprocal.reset();
    myIntra.reset();
    mySurface.reset();
#endif
  }


  //_______________________________________________________________________ destructor
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,
				 TCellManager,
				 real,
				 reciprocal,
				 correction,
				 TInterpolation,
				 TSwitchingFunction>::~iSGNonbondedPMEwaldSystemForce(){
#if defined(DEBUG_PME_TIMING)
    if(boundaryConditions.getVolume() > EPSILON){
      report.setf(std::ios::showpoint|std::ios::fixed);
      report << allnodes << plain <<"Timing ("
	     <<Parallel::getId()<<") PME:"
	     <<" real:" <<myReal.getTime().getProcessTime()<<"[s]"
	     <<", reciprocal:"<<myReciprocal.getTime().getProcessTime()<<"[s]"
	     <<", intra:"<<myIntra.getTime().getProcessTime()<<"[s]"
	     <<", dipole:"<<mySurface.getTime().getProcessTime()<<"[s]"
	     <<", Rc="<<myRc
	     <<", Kc="<<myKc
	     <<", n="<<myLattice.size()
	     <<", alpha="<<myAlpha
	     <<", accuracy="<<myAccuracy<<"."<<endr;
    }
#endif
  }
  
	    
  //_______________________________________________________________________ evaluate function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,
				      TCellManager,
				      real,
				      reciprocal,
				      correction,
				      TInterpolation,
				      TSwitchingFunction>::evaluate(const GenericTopology* topo, 
								    const Vector3DBlock* positions, 
								    Vector3DBlock* forces, 
								    ScalarStructure* energies) {
    
    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

    // Initialize data members and precompute tables & short cuts
    if(!myCached)
      initialize(realTopo,positions);
    
    // Intra-molecular and surface dipole term
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
      pointSelfEnergy     = myPointSelfEnergy;
      chargedSystemEnergy = myChargedSystemEnergy;
      pointSelfDeltaMu     = myPointSelfDeltaMu;
      chargedSystemDeltaMu = myChargedSystemDeltaMu;
      (*energies)[ScalarStructure::VIRIALXX] += myChargedSystemEnergy;
      (*energies)[ScalarStructure::VIRIALYY] += myChargedSystemEnergy;
      (*energies)[ScalarStructure::VIRIALZZ] += myChargedSystemEnergy;
      
      if(false) surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
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
      reciprocalTerm(realTopo,positions,forces,energies,reciprocalEnergy,reciprocalDeltaMu);

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

#if defined(DEBUG_PME_ENERGIES)
    report.setf(std::ios::showpoint|std::ios::fixed);
    report << plain <<"PME: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   <<endr;
    report << plain <<"PME: point="<<myPointSelfDeltaMu
	   <<", charged="<<myChargedSystemDeltaMu
	   <<", real="<<realDeltaMu
	   <<", reciprocal="<<reciprocalDeltaMu
	   <<", intra="<<intraMolecularDeltaMu
	   <<", total="<<d_mu
	   << endr;
#endif
  }


  //_______________________________________________________________________ parallelEvaluate function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,TSwitchingFunction>::parallelEvaluate(const GenericTopology* topo, 
														      const Vector3DBlock* positions, 
														      Vector3DBlock* forces, 
														      ScalarStructure* energies) {
    
    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint)
      Report::report << Report::donthint;

    const RealTopologyType* realTopo = (const RealTopologyType *)topo;  

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
	  if(Parallel::next()){
	    unsigned int to = (n*(i+1))/count;
	    if(to > n) to = n;
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
	if(false) surfaceDipoleTerm(realTopo,positions,forces,energies,surfaceDipoleEnergy);
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
	if(Parallel::next()) realTerm(realTopo,positions,forces,energies,realEnergy,realDeltaMu,l);
	else enumerator.nextNewPair(l);	
      }
    }  
    
    // Reciprocal-space term
    Real reciprocalEnergy = 0.0;
    Real reciprocalDeltaMu = 0.0;
    if(reciprocal){
      reciprocalTermParallel(realTopo,positions,forces,energies,reciprocalEnergy,reciprocalDeltaMu);      
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

#if defined(DEBUG_PME_ENERGIES)
    report << allnodes << plain <<"PME: point="<<myPointSelfEnergy
	   <<", charged="<<myChargedSystemEnergy
	   <<", real="<<realEnergy
	   <<", reciprocal="<<reciprocalEnergy
	   <<", intra="<<intraMolecularEnergy
	   <<", surface="<<surfaceDipoleEnergy
	   <<", total="<<e 
	   <<endr;
    report << allnodes << plain <<"PME: point="<<myPointSelfDeltaMu
	   <<", charged="<<myChargedSystemDeltaMu
	   <<", real="<<realDeltaMu
	   <<", reciprocal="<<reciprocalDeltaMu
	   <<", intra="<<intraMolecularDeltaMu
	   <<", total="<<d_mu
	   << endr;
#endif

    if(dontHint)
      Report::report << Report::dohint;

  }


  //_______________________________________________________________________ numberOfBlocks function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  unsigned int iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
					      real,reciprocal,correction,TInterpolation,TSwitchingFunction>::numberOfBlocks(const GenericTopology* topo,
															    const Vector3DBlock* positions){
		  
    unsigned int n = 0;
    if(correction)
      n +=std::min(static_cast<int>(topo->exclusions.getTable().size()),static_cast<int>(Parallel::getAvailableNum()));
    
    if(reciprocal)
      ;
    
    if(real){
      const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
      realTopo->updateCellLists(positions);
      n += Parallel::getNumberOfPackages(realTopo->cellLists.size());
    }
    
    return n;
  }


  //_______________________________________________________________________ initialize function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::initialize(const RealTopologyType* realTopo, const Vector3DBlock* positions) {
    
    bool dontHint = (myV >= 0.0 && TBoundaryConditions::PERIODIC);
    if(dontHint) Report::report << Report::donthint;
    
    if(TBoundaryConditions::VACUUM){
      // We are have a non-periodic case. Do an expansion of the simulation 
      // box by the expansion factor, if need
      realTopo->updateCellLists(positions);
      Vector3D d2 = (boundaryConditions.getMax()-boundaryConditions.getMin());
      Vector3D d0 = d2/myExpansionFactor;
      Vector3D d1 = realTopo->max-realTopo->min;
      if (fabs(1-d0.x/d1.x) > 0.1 || fabs(1-d0.y/d1.y) > 0.1 || fabs(1-d0.y/d1.y) > 0.1 || 
          d2.x <= d1.x || d2.y <= d1.y || d2.z <= d1.z){
	// The boundaries are to big or small, adjust
	Vector3D e = (realTopo->max-realTopo->min)*myExpansionFactor;
	Vector3D origin = realTopo->min*myExpansionFactor + e*0.5;
	boundaryConditions.set(Vector3D(e.x,0,0), Vector3D(0,e.y,0), Vector3D(0,0,e.z), origin);
	report << hint << "The boundaries for PME force evaluation were re-sized to "<<boundaryConditions.getMin()<<"-"<<boundaryConditions.getMax()<<"."<<endr;
      }
      else {
	// The particles are just inside +-10% our previous simulation box. No update needed
	myV = boundaryConditions.getVolume();
	return;      
      }
    }  // end if VACUUM boundary conditions case
    else {
      boundaryConditions.set(realTopo->boundaryConditions.e1(),
			     realTopo->boundaryConditions.e2(),
			     realTopo->boundaryConditions.e3(),
			     realTopo->boundaryConditions.origin());
    } // end if PERIODIC boundary conditions case
    
    if(!boundaryConditions.isOrthogonal())
      report << error << "[NonbondedFullEwaldSystemForce::initialize] Not orthogonal, aborting."<<endr;

    // get the number of atoms in the system
    const unsigned int atomCount = realTopo->atoms.size();
    
    
    myTRate     = 5.5;     // From Moldy, rate between real and reciprocal
    //myAccuracy  = 0.00001; // From Moldy, accuracy
    //myAccuracy  = 1.e-6; // From NAMD, accuracy
  
    // Dimension of the simulation box and box volume
    myLX   = boundaryConditions.e1().x;
    myLY   = boundaryConditions.e2().y;
    myLZ   = boundaryConditions.e3().z;
    myV    = boundaryConditions.getVolume();

    // reciprocal cell basis vectors and reciprocal box volume
    myLXr  = boundaryConditions.e1r().x;
    myLYr  = boundaryConditions.e2r().y;
    myLZr  = boundaryConditions.e3r().z;
    myVr   = 1.0/myV;

    // origin for the FFT grid
    myOrigin = boundaryConditions.origin();

    // Short cuts
    if(myAlphaDefault){
      myAlpha = 1.0;
      //myRc = 6.5;
      while ( erfc(myAlpha*myRc)/myRc >= myAccuracy ) 
        myAlpha *= 2.0;

      Real low = 0.;
      Real high = myAlpha;
      for(unsigned int i=0; i<100; ++i){
	myAlpha = 0.5 * (low + high);
	if ( erfc(myAlpha*myRc)/myRc >= myAccuracy ) {low = myAlpha;}
	else {high = myAlpha;}
      }
      //myAlpha         = sqrt(M_PI)*pow(myTRate*atomCount/(power<2>(myV)),1.0/6.0);
    }
    Real p          = -log(myAccuracy);
    myKc            = 2.0*myAlpha*sqrt(p);
    //myRc            = sqrt(p)/myAlpha;
    myRcSquared     = myRc*myRc;
    
    // compute powers of alpha
    myAlphaSquared  = myAlpha*myAlpha;
    myAlphaSquaredr = 1.0/myAlphaSquared;
    my2AlphaPI      = 2.0*myAlpha/sqrt(M_PI);
    
    // Initialize the charge distribution grids
    if(reciprocal) {
      myISGGrid.initialize(myLX,myLY,myLZ,myAlpha,myNX,myNY,myNZ,myInterOrder,atomCount);
    }

    // ___________________________________________________
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
    
    // ___________________________________________________
    // Charged system energy
    //
    myChargedSystemEnergy  = 0.0;
    myChargedSystemDeltaMu = 0.0;
    if(correction){
      Real q = 0.0;
      Real DeltaMu_q = 0.0;
      for(unsigned int i=0;i<atomCount;i++){

	// compute the net charge in the system
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
    
    report << hint <<"PME";
#ifdef USE_PME_EXACT_ERF
    report << "(Exact)";
#endif  
    report <<": alpha="<<myAlpha<<", V="<<myV<<", Rc="
	   <<myRc<<", Kc ="<<myKc<<", n="<<myLattice.size()<<", accuracy="
	   <<myAccuracy<<", interpolation="<<TInterpolation::keyword
	   <<", order="<<myInterOrder<<"."<<endr;

    // Now we have all pre-computed stuff
    if(TBoundaryConditions::PERIODIC)
      myCached = true;

    if(dontHint)
      Report::report << Report::dohint;
  }

  //_______________________________________________________________________ realTerm function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::realTerm(const RealTopologyType* realTopo,
								    const Vector3DBlock* positions, 
								    Vector3DBlock* forces, 
								    ScalarStructure* energies,
								    Real& realEnergy,
								    Real& realDeltaMu,
								    unsigned int n) {
#if defined(DEBUG_PME_TIMING)
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
	if(count > n) break;	
      }
      
      // loop over all pairs in the cell list
      // outer loop over atom i
      for(int i=thisPair.first; i!=-1; i=realTopo->atoms[i].cellListNext){
	
	// the charge on atom i
	Real qi  = realTopo->atoms[i].scaledCharge;
	
	// the xyz coordinates of atom i, force on atom i, and its molecule ID #
	Vector3D ri((*positions)[i]), fi;
        int mi = realTopo->atoms[i].molecule;
	
	// inner loop over atom j				
	for(int j=(notSameCell ? thisPair.second:i); j!=-1; j=realTopo->atoms[j].cellListNext){
  
          // molecule ID# for atom j
          int mj = realTopo->atoms[j].molecule;

	  // check to see if this pair should be excluded
          bool same = (mi==mj);
	  ExclusionClass excl = (same?realTopo->exclusions.check(i,j):EXCLUSION_NONE);
	  
	  // exclude any self interactions (i.e. atom interacting with it's own periodic image)
	  if(i == j) excl = EXCLUSION_FULL;
	  
	  for(unsigned int k=0;k<myLattice.size();k++){
	    // Check for an exclusion.
	    if (excl == EXCLUSION_FULL){ 
	      excl = EXCLUSION_NONE;
	      continue;	
	    }
	    
	    // minimal image vector from atom j to atom i (rj - ri)
	    Vector3D rijMinimal(boundaryConditions.minimalDifference(ri,(*positions)[j]));
	    Vector3D rij(rijMinimal+myLattice[k]);
	    
	    // separation distance squared
	    Real rSquared = rij.normSquared();
	    
	    // Do switching function rough test.
	    if (rSquared>myRcSquared) continue;
	    
	    // we must determine which type of pair interaction this is:
	    // type 0 = untransformed-untransformed
	    // type 1 = untransformed-transformed
	    // type 2 = transformed-transformed (intramolecular)
	    bool atom1_scaled, atom2_scaled;
	    atom1_scaled = atom2_scaled = false;
	    int myStage, OldStage /*, atom_stage*/;
            myStage = OldStage = /*atom_stage =*/ 0;
                        
	    // Use the molecule types to determine if the atoms are being transformed or not	    
	    // if lambda for a molecule is exactly zero then it is not a scaled/transforming molecule
	    if (realTopo->molecules[mi].lambda != 0.0) {atom1_scaled = true;}
	    if (realTopo->molecules[mj].lambda != 0.0) {atom2_scaled = true;}
	    
	    // integer used to select which interaction type to compute   
	    int Choice;
	    
	    // determine the value of Choice
            // neither atom is being transformed
            if ( !(atom1_scaled) && !(atom2_scaled) ) {Choice = 0;}
            // intramolecular interaction on a transforming molecule
            else if ( (atom1_scaled) && (atom2_scaled) ) Choice = 2;
            // one of the atoms is being transformed
            else Choice = 1;
      
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
#ifdef USE_PME_EXACT_ERF
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

	      // product of charges
	      qq = qi*qj;	
	      
	      // scale any 1-4 interactions if desired
	      if (excl == EXCLUSION_MODIFIED) qq *= realTopo->coulombScalingFactor;
	      
	      // multiply qq by the complementary error function
	      // Approximation Abramowitz & Stegun p299.
	      // Energy and force
#ifndef USE_PME_EXACT_ERF
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

	      // product of charges
	      qq = qi*qj;

	      // derivative of qq with respect to lambda	
	      q_Dq = (qi * realTopo->atoms[j].deltaQ
		      + qj * realTopo->atoms[i].deltaQ);
	      
	      // multiply qq and q_Dq by the complementary error function
	      // Approximation Abramowitz & Stegun p299.
	      // Energy
#ifndef USE_PME_EXACT_ERF
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
              OldStage = myStage - 1;
        
	      // Energy
	      // The interaction energy is the sum of the interaction in the
	      // old state scaled by (1 - lambda) plus the interaction in the new state
	      // scaled by lambda
	      
	      // product of charges for each identity
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
	      
	      // scale any 1-4 interactions if desired					
              if (excl == EXCLUSION_MODIFIED) {
		qq *= realTopo->coulombScalingFactor;
		q_Dq *= realTopo->coulombScalingFactor;
		qq_correction_add *= realTopo->coulombScalingFactor;
		q_Dq_correction_add *= realTopo->coulombScalingFactor;
              }
              else {
		// this is an unmodified intramolecular interaction, so subtract the total
		// reciprocal space energy and chemical potential difference
		qq_correction_minus = qi * realTopo->atoms[j].scaledCharge;
		q_Dq_correction_minus = qi * realTopo->atoms[j].deltaQ + qj * realTopo->atoms[i].deltaQ;
              }

              // multiply qq and q_Dq by the complementary error function
	      // Approximation Abramowitz & Stegun p299.
	      // Energy
#ifndef USE_PME_EXACT_ERF
	      rr = 1.0/r;
	      ar = myAlpha*r;
	      e = qq*exp(-ar*ar);
	      dmu = q_Dq*exp(-ar*ar);	
	      energy = poly5(ar)*e*rr;
	      deltaMu = poly5(ar)*dmu*rr;
	      force = ((energy+my2AlphaPI*e)*rr*rr);
	      // Reciprocal space correction term gets
	      // multiplied by the error function
	      e_correction = erf(myAlpha*r)*rr;
#else	  
	      a = erfc(myAlpha*r)/r;
	      energy = qq*a;
	      deltaMu = q_Dq*a;
	      force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
	      // Reciprocal space correction term gets
	      // multiplied by the error function
	      e_correction = erf(myAlpha*r)/r;
#endif

	      // add in the reciprocal space correction terms
	      myCorrection = e_correction * (qq_correction_add - qq_correction_minus);
	      myDMUcorrection = e_correction * (q_Dq_correction_add - q_Dq_correction_minus);
	      energy += myCorrection;
	      deltaMu += myDMUcorrection;
	      break;
	      //----------------------------------------------------------
            } // end switch structure
	    
	    
	    // Calculate the switched force, energy, and deltaMu.
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
	      energies->addVirial(fij,rij,realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
									                 realTopo->molecules[mj].position));
	    }
	    else if(doVirial) 
	      energies->addVirial(fij,rij);
	 
	    excl = EXCLUSION_NONE;
          }
        }
        (*forces)[i] += fi;
      }
    }
#if defined(DEBUG_PME_TIMING)
    myReal.stop();
#endif
  }


  //_______________________________________________________________________ reciprocalTermParallel function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::reciprocalTermParallel(const RealTopologyType* realTopo,
										  const Vector3DBlock* positions, 
										  Vector3DBlock* forces, 
										  ScalarStructure* energies,
										  Real& reciprocalEnergy,
										  Real& reciprocalDeltaMu) {
    
#if defined(DEBUG_PME_TIMING)
    myReciprocal.start();
#endif
    const unsigned int atomCount = realTopo->atoms.size();
    const unsigned int nBlocks = Parallel::getAvailableNum();
    const unsigned int block = Parallel::getAvailableId();
    const unsigned int i0 = (atomCount*block)/nBlocks;
    const unsigned int i1 = (atomCount*(block+1))/nBlocks;
    Real* begin = NULL;
    Real* end   = NULL;
    
    // clear the charge distribution grids
    myISGGrid.clear();

    // loop over all atoms
    for(unsigned int i=i0;i<i1;i++) {
      
      // add the atomic charge and deltaQ to the grid
      myISGGrid.anterpolateCharge(realTopo->atoms[i].scaledCharge,realTopo->atoms[i].deltaQ,(*positions)[i],i); 
    }

    
    myISGGrid.getQ(begin,end);
    Parallel::reduceSlaves(begin,end);

    // perform the inverse FFT operation on the two grids
    if(FFTComplex::isParallel()) {myISGGrid.fftBack();}
    else {
      if(block == 0) {myISGGrid.fftBack();}
      myISGGrid.getQ(begin,end);
      Parallel::bcastSlaves(begin,end);
    }
    
    // sum over the transformed charge grids to compute the energy
    // and chemical potential difference
    Real energy, deltaMu;
    myISGGrid.scalarSum(energies,energy,deltaMu,block,nBlocks);
    
    myISGGrid.getQ(begin,end);
    Parallel::reduceSlaves(begin,end);

    // take the product of the transformed Q grid with the arrays B and C (see Essmann et. al.)
    // and perform the forward FFT operation to get the convolution (Theta * Q)[m1,m2,m3] 
    if(FFTComplex::isParallel()) {myISGGrid.fftForward();}
    else{
      if(block == 0) {myISGGrid.fftForward();}
      myISGGrid.getQ(begin,end);
      Parallel::bcastSlaves(begin,end);
    }

    // loop over all atoms
    bool doMolVirial = energies->molecularVirial();
    for(unsigned int i=i0;i<i1;i++) {

      // compute the force on atom i from the reciprocal space part
      Vector3D fi;
      myISGGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,fi);
      (*forces)[i] += fi;

      if(doMolVirial){
        // compute the vector from atom i to the center of mass of the molecule
        Vector3D ria = boundaryConditions.minimalPosition((*positions)[i]);
        Vector3D mri = realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[realTopo->atoms[i].molecule].position);

        // compute the reciprocal space contribution to the molecular virial
        // this expression is taken from Darden, et al. J. Chem. Phys. 103 (19), 8577.
        energies->addMolVirial(fi,mri);
      }
    } // end loop over atoms (i)

   
    reciprocalEnergy += energy;
    reciprocalDeltaMu += deltaMu;
  
#if defined(DEBUG_PME_TIMING)
    myReciprocal.stop();
#endif
  }


  //_______________________________________________________________________ reciprocalTerm function  
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::reciprocalTerm(const RealTopologyType* realTopo,
									  const Vector3DBlock* positions, 
									  Vector3DBlock* forces, 
									  ScalarStructure* energies,
									  Real& reciprocalEnergy,
									  Real& reciprocalDeltaMu) {
    
#if defined(DEBUG_PME_TIMING)
    myReciprocal.start();
#endif
    
    // get the number of atoms in the system
    const unsigned int atomCount = realTopo->atoms.size();
    
    // clear the charge distribution grids   
    myISGGrid.clear();

    // loop over all atoms
    for(unsigned int i=0;i<atomCount;i++) {

      // add the atomic charge and deltaQ to the grid
      myISGGrid.anterpolateCharge(realTopo->atoms[i].scaledCharge,
				  realTopo->atoms[i].deltaQ,(*positions)[i],i);
    }

    // perform the inverse FFT operation on the charge distribution grids
    myISGGrid.fftBack();

    // sum over the transformed charge grids to compute the energy
    // and chemical potential difference				
    Real energy, deltaMu;
    myISGGrid.scalarSum(energies,energy,deltaMu);

    // take the product of the transformed Q grid with the arrays B and C (see Essmann et. al.)
    // and perform the forward FFT operation to get the convolution (Theta * Q)[m1,m2,m3]   
    myISGGrid.fftForward();

    // loop over all atoms
    bool doMolVirial = energies->molecularVirial();
    for(unsigned int i=0;i<atomCount;i++) {

      // compute the force on atom i from the reciprocal space part
      Vector3D fi;
      myISGGrid.interpolateForce(realTopo->atoms[i].scaledCharge,i,fi);
      (*forces)[i] += fi;

      if(doMolVirial){
        // get the ID# of the molecule to which this atom belongs
	int Mi = realTopo->atoms[i].molecule;

	// compute the vector from atom i to the center of mass of the molecule
	Vector3D ria(boundaryConditions.minimalPosition((*positions)[i]));
	Vector3D mri(realTopo->boundaryConditions.minimalDifference(ria,realTopo->molecules[Mi].position));

	// compute the reciprocal space contribution to the molecular virial
	// this expression is taken from Darden, et al. J. Chem. Phys. 103 (19), 8577.
	energies->addMolVirial(fi,mri);
      }
    } // end loop over atoms

    reciprocalEnergy += energy;
    reciprocalDeltaMu += deltaMu;
  
#if defined(DEBUG_PME_TIMING)
    myReciprocal.stop();
#endif
  }


  //_______________________________________________________________________ correctionTerm function  
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::correctionTerm(const RealTopologyType* realTopo,
									  const Vector3DBlock* positions, 
									  Vector3DBlock* forces, 
									  ScalarStructure* energies,
									  Real& intraMolecularEnergy,
									  Real& intraMolecularDeltaMu,
									  unsigned int from, unsigned int to){
    
#if defined(DEBUG_PME_TIMING)
    myIntra.start();
#endif

    // Intra-molecular term
    bool doVirial = energies->virial();
    const std::vector<ExclusionPair>& exclusions = realTopo->exclusions.getTable();

    // loop over all excluded pairs
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
      
      // separation distances
      Real r  = sqrt(rSquared);
      Real rr = 1/r;

      // error function
      Real e = erf(myAlpha*r)*rr;
      
      // Intra-molecular self energy and chemical potential difference
      intraMolecularEnergy -= qq*e;
      intraMolecularDeltaMu -= q_Dq*e;
      
      // Intra-molecular self force
      Vector3D fij(rij*(qq*(my2AlphaPI*exp(-myAlphaSquared*rSquared)-e)*rr*rr));
      (*forces)[excl.a1] -= fij;
      (*forces)[excl.a2] += fij;
      if (doVirial)	
        energies->addVirial(fij,rij);
    }
#if defined(DEBUG_PME_TIMING)
    myIntra.stop();
#endif
  }


  //_______________________________________________________________________ surfaceDipoleTerm function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,
				      TSwitchingFunction>::surfaceDipoleTerm(const RealTopologyType* realTopo,
									     const Vector3DBlock* positions, 
									     Vector3DBlock* forces, 
									     ScalarStructure* energies,
									     Real& surfaceDipoleEnergy){
				
#if defined(DEBUG_PME_TIMING)
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
  
#if defined(DEBUG_PME_TIMING)
    mySurface.stop();
#endif
  }

  //_______________________________________________________________________ getIdNoAlias function
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  std::string iSGNonbondedPMEwaldSystemForce<TBoundaryConditions, TCellManager,
					     real,reciprocal,correction,TInterpolation,TSwitchingFunction>::getIdNoAlias() const{
		    
    return (iSGCoulombForce::keyword + " -algorithm "+ keyword + 
	    std::string((real)       ? std::string(" -real")       : std::string("")) + 
	    std::string((reciprocal) ? std::string(" -reciprocal") : std::string("")) +
	    std::string((correction) ? std::string(" -correction") : std::string("")) +
	    " -interpolation "+TInterpolation::getKeyword()+
	    std::string((TSwitchingFunction::getId() != CutoffSwitchingFunction::getId()) ? 
			std::string(std::string(" -switchingFunction " + TSwitchingFunction::getId())) : std::string("")));
  }


  //_______________________________________________________________________ getParameters function 
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  void iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
				      real,reciprocal,correction,TInterpolation,TSwitchingFunction>::getParameters(std::vector<Parameter>& parameters) const{
		
    parameters.push_back(Parameter("-gridsize",Value(myNX,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("",Value(myNY,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("",Value(myNZ,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("-cutoff",Value(myRc,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("-order",Value(myInterOrder,ConstraintValueType::Positive()),4));
    parameters.push_back(Parameter("-accuracy",Value(myAccuracy,ConstraintValueType::Positive()),1.e-6));
    parameters.push_back(Parameter("-alpha",Value(myAlpha),-1.0));    
    if(TBoundaryConditions::VACUUM) parameters.push_back(Parameter("-j",Value(myExpansionFactor),3.0));
  }


  //_______________________________________________________________________ doMake function  
  template <class TBoundaryConditions, 
	    class TCellManager,
	    bool  real,
	    bool  reciprocal,
	    bool  correction,
	    class TInterpolation,
	    class TSwitchingFunction>
  Force* iSGNonbondedPMEwaldSystemForce<TBoundaryConditions,TCellManager,
					real,reciprocal,correction,TInterpolation,TSwitchingFunction>::doMake(std::string& errMsg,std::vector<Value> values) const{

    int  nx              = values[0];
    int  ny              = values[1];
    int  nz              = values[2];
    Real cutoff          = values[3];
    int  order           = values[4];
    Real accuracy        = values[5];
    Real alpha           = values[6];
    Real expansionFactor = (TBoundaryConditions::VACUUM?(Real)values[7]:3.0);
    std::string err = "";

    if(TBoundaryConditions::VACUUM && !values[7].valid())
      err +=" expansionFactor \'"+values[2].getString()+"\' for "+getId()+" not valid.";
    if(expansionFactor <= 1.0)
      err += keyword + " simulation box expansion factor (="+toString(expansionFactor)+") > 1.0.";

    if(order < 2 || (order % 2) != 0 || !values[4].valid())
      err += keyword + " order (="+values[4].getString()+") > 1 and even.";

    if(!values[0].valid() || !values[1].valid() || !values[2].valid() || nx < order || ny < order || nz < order)
      err += keyword + " force: "+values[4].getString()+" <= nx (="+values[0].getString()+"), "+values[4].getString()+" <= ny (="+values[1].getString()+"), "+values[4].getString()+" <= nz (="+values[2].getString()+").";

    if(cutoff <= 0.0 || !values[3].valid())
      err += keyword + " cutoff (="+values[3].getString()+") > 0.";

    if(accuracy <= 0.0 || !values[5].valid())
      err += keyword + " accuracy (="+values[5].getString()+") > 0.";

    if(!values[6].valid())
      err +=" alpha \'"+values[6].getString()+"\' not valid.";

    if(!err.empty()){
      errMsg += " force "+keyword+" :"+err;
      return NULL;
    }

    return (new iSGNonbondedPMEwaldSystemForce((unsigned int)nx,(unsigned int)ny,(unsigned int)nz,(unsigned int)order,cutoff,accuracy,alpha,expansionFactor));
  }
}
#endif /* ISGNONBONDEDPMEWALDSYSTEMFORCE_H */
