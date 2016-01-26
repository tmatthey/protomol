/*  -*- c++ -*-  */
#ifndef OUTPUTCACHE_H
#define OUTPUTCACHE_H

#include "Real.h"
#include "PDB.h"
#include "PAR.h"
#include "iSGPAR.h"
#include "PSF.h"
#include "Vector3D.h"

namespace ProtoMol {
  class Output;
  class Configuration;
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class OutputFactoryDetails;
  class Integrator;
  struct PDB;
  struct Atom;
  /**
     OutputCache caches all kind of values, which may be needed
     by Output objects and simplifies the access to values of interest.
     Add new cached values, if needed ..
     There are some (feature) values, which will only change when the 
     Topology changes
  */
  //________________________________________________________ OutputCache
  class OutputCache  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputCache();
    ~OutputCache();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputCache
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
    		    const Vector3DBlock* pos, const Vector3DBlock*  vel, const ScalarStructure* energies);

    // Methods to add additional data for output objects
    void add(const std::vector<PDB::Atom>& pdbAtoms){myAtoms = pdbAtoms;}
    void add(const PSF& psf){myPSF = psf;}
    void add(const PAR& par){myPAR = par;}
    void add(const iSGPAR& par){myiSGPAR = par;}
    void add(const std::vector<Real> &REMExchangeRate) {myREMRates = REMExchangeRate;}
    void add(const Real *replicaHistory) {myReplicaHistory = replicaHistory;}

    Real     totalEnergy() const;
    Real     potentialEnergy() const;
    Real     kineticEnergy() const;
    Real     temperature() const;
    Real     temperatureForWater() const;
    Real     temperatureForNonWater() const;
    Real     volume() const;
    Real     time() const;
    Real     pressure() const;
    Real     molecularPressure() const;
    Real     molecularTemperature() const;
    Real     molecularKineticEnergy() const;
    Vector3D linearMomentum() const;
    Vector3D angularMomentum() const;
    Vector3D centerOfMass() const;
    Real     diffusion() const;
    Real     density() const;
    Real     mass() const;
    Real     dihedralPhi(int index) const;
    Real     brent(Real ax, Real bx, Real cx, Real tol, Real &xmin, int dihindex, bool max) const;
    std::vector<Real>                 dihedralPhis(std::vector<int>) const;
    std::vector<std::vector< Real > > brentMaxima(std::vector<int>, bool) const;
    const Vector3DBlock*             minimalPositions() const;
    const Vector3DBlock*             PositionsNoWater() const;

    const std::vector<PDB::Atom>& pdb() const {return myAtoms;}
    const std::vector<Real> &REMRates() const {return myREMRates;}
    const Real        *replicaHistory() const {return myReplicaHistory;}
    const PSF&                       psf() const {return myPSF;}
    const PAR&                       par() const {return myPAR;}
    const iSGPAR&                    iSGpar() const {return myiSGPAR;}
    void uncache() const;
    ///< To be called before every run() or finialize()
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const Vector3DBlock* getPositions() const {return myPositions;}
    const Vector3DBlock* getVelocities() const  {return myVelocities;}
    const Integrator* getIntegrator() const {return myIntegrator;}
    const ScalarStructure* getEnergies() const {return myEnergies;}

    void setRestore() {myRestore = true;}
    void clearRestore() {myRestore = false;}
    bool restore() const {return myRestore;}

  private:
    const Configuration*   myConfig; 
    const GenericTopology* myTopology; 
    const Integrator*      myIntegrator;
    const ScalarStructure* myEnergies;
    const Vector3DBlock*   myPositions;
    const Vector3DBlock*   myVelocities; 
    Vector3DBlock*         myInitialPositions;
    mutable Vector3DBlock*   myMinimalPositions;
    mutable Vector3DBlock*   myPositionsNoWater;    

    // Additional data
    std::vector<PDB::Atom> myAtoms;
    std::vector<Real> myREMRates;
    const Real *myReplicaHistory;
    PSF myPSF;
    PAR myPAR;
    iSGPAR myiSGPAR;

    mutable bool myCachedKE;
    mutable Real myKE;
    mutable Real myT;

    mutable bool myCachedPE;
    mutable Real myPE;

    mutable bool myCachedV;
    mutable Real myV;

    mutable bool myCachedP;
    mutable Real myP;
    mutable bool myCachedMolP;
    mutable Real myMolP;

    mutable bool myCachedLinearMomentum;
    mutable Vector3D myLinearMomentum;

    mutable bool myCachedAngularMomentum;
    mutable Vector3D myAngularMomentum;

    mutable bool myCachedCenterOfMass;
    mutable Vector3D myCenterOfMass;

    mutable bool myCachedDiffusion;
    mutable Real myDiffusion;

    mutable bool myCachedDensity;
    mutable Real myDensity;

    mutable bool myCachedMass;
    mutable Real myMass;

    mutable int myCachedDihedralPhi;
    mutable Real myDihedralPhi;

    mutable bool myCachedDihedralPhis;
    mutable std::vector<Real>* myDihedralPhis;

    mutable bool myCachedBrentMaxima;
    mutable std::vector< std::vector< Real > >* myBrentMaxima;

    mutable bool myCachedMolT;
    mutable Real myMolT;

    mutable bool myCachedMolKE;
    mutable Real myMolKE;

    mutable bool myCachedWaterT;
    mutable Real myWaterT;

    mutable bool myCachedNonWaterT;
    mutable Real myNonWaterT;

    mutable bool myCachedMinimalPositions;

    mutable bool myCachedPositionsNoWater;

    bool myRestore;
  };
}
#endif
