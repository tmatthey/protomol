/*  -*- c++ -*-  */
#ifndef OSGINTEGRATOR_H
#define OSGINTEGRATOR_H

#include "STSIntegrator.h"
#include "XSC.h"

#include "PAR.h"
#include "PSF.h"
#include "STAGE.h"
#include "Array.h"

namespace ProtoMol {

  //________________________________________________________ oSGIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;
  class ModifierISGforces;
  class Modifier;

  template<class TIntegrator>
  class ModifierPreForceThermostat;
  template<class TIntegrator>
  class ModifierPostForceThermostat;
  template<class TIntegrator>
  class ModifierPreForceBarostat;
  template<class TIntegrator>
  class ModifierPostForceBarostat;

  class oSGModifierPreForceChemostat;
  class oSGModifierPostForceChemostat;

  class oSGIntegrator: public STSIntegrator {

    ///friend class ModifierOSG;

    template<class TIntegrator>
    friend class ModifierPreForceThermostat;
    template<class TIntegrator>
    friend class ModifierPostForceThermostat;
    template<class TIntegrator>
    friend class ModifierPreForceBarostat;
    template<class TIntegrator>
    friend class ModifierPostForceBarostat;

    friend class oSGModifierPreForceChemostat;
    friend class oSGModifierPostForceChemostat;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    oSGIntegrator();
    oSGIntegrator(Real timestep,
                  unsigned int numComp,
		  Real temperature,
                  Real pressure,
                  Real tauT,
		  Real tauV,
                  Real tauP,
		  Real tauD,
                  Real tauL,
                  Real MuTemp,
                  ForceGroup *overloadedForces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // new methods of class oSGIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Real getEpsilonVel() const {return myEpsilonVel;}
    Real getEtaVel() const {return myEtaVel;}
    Real getNumAtoms() const {return NumAtoms;}
    Real getAveDeltaMu(int NumSteps) {return (AveDeltaMu / NumSteps);}
    Real getAveLambdaT(int NumSteps) {return (LambdaT / static_cast<Real>(NumSteps));}
    void setTrajectoryFreq(int freq) {myTrajectoryFreq = freq;}

    // store the information for the transformMaps, stock molecules, and topology
    void indexTypes(GenericTopology *topo, const STAGE &tstage, const PAR& par, bool dihedralMultPSF, int theSeed);

    //  If an initial XSC file is specified, read in the values
    //  of the extended system coordinates
    virtual void readXSCs(const std::string myFile, GenericTopology *topo);
    
    // check to see if a transformation stage has been completed
    void checkForStageCompletion();
 
    // for successful transformations, update the topology
    void setForcesAfterTransformation();
    
    // this function decides iiif we will insert or delete a molecule
    void pickNewMolecule();

  protected:
    virtual void do2ndHalfKick()=0;
    void PreForceThermostat();
    void PostForceThermostat();
    void PreForceBarostat();
    void PostForceBarostat();
    virtual void PreForceChemostat();
    virtual void PostForceChemostat();
    virtual void deletionDeltaMu(const Real myFugacity, const int Osm)=0;
    virtual void insertionDeltaMu(const Real myFugacity, const int Osm)=0;
    
  private:
    void modifyForces();
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getScope() const{return scope;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
                            Vector3DBlock   *positions,
                            Vector3DBlock   *velocities,
                            ScalarStructure *energies);
    virtual void run(int numTimesteps);

  protected:
    virtual void addModifierBeforeInitialize();
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doDrift()=0;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doHalfKick()=0;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string scope;

  protected: 
    const unsigned int myNumComp; //  Number of components in the mixture
    unsigned int myNumStages;     //  Number of stages that the transformation will be split into
    const Real myTargetTemp;      //  Target temperature.  Units: (K)
    const Real myTargetPres;      //  Target pressure.  Units: (bar)
    const Real myMuTemp;          //  Target chemostat temperature.  Units: (K)
    const Real myTauT;            //  thermostat oscillation period.  Units: (fs)
    const Real myTauV;            //  volume thermostat oscillation period.  Units: (fs)
    const Real myTauP;            //  barostat oscillation period.  Units: (fs)
    const Real myTauD;            //  chemostat oscillation period.  Units: (fs)
    const Real myTauL;            //  chemostat thermostat oscillation period.  Units: (fs)
    const Real kbT;               //  Target temperature * Boltzmann's constant.  Units: (kcal/mol)
    Real myLambda;                //  the chemostat variable.  Units:  (dimensionless)
    unsigned int NumAtoms;        //  Total # of atoms in the system.
    unsigned int NumMols;         //  Total # of molecules in the system
    unsigned int myNumFree;       //  Total # of degrees of freedom = (3*Natoms - 3) - NumConstraints
    int T;                        //  The ID# of the molecule currently being transformed
    int thisStage;                //  The current transformation stage
    unsigned int NumTransSteps;   //  The # of timesteps it took to complete the transformation
    unsigned int OldType;         //  The type of molecule being transformed (species ID#)
    bool Insert;                  //  Indicates if the transforming molecule is being inserted or deleted
    bool Transformed;             //  flag so we know when the molecule has been completely transformed
    bool Succeeded;               //  flag so we know if the transformation attempt succeeded
    Real myTargetMu;              //  residual chemical potential.  Units (kcal/mol)
    Real myVolume;                //  Current cubic volume.  Units (AA^3)
    Real myEpsilonVel;            //  Barostat strain rate velocity.  Units: (fs)^-1
    Real Qo;                      //  Particle thermostat mass.  Units: (kcal fs^2 / mol)
    Real Qv;                      //  Volume thermostat mass.  Units: (kcal fs^2 / mol)
    Real W;                       //  Barostat mass.  Units: (kcal fs^2 / mol)
    Real Qd;                      //  Chemostat mass.  Units: (kcal fs^2 / mol)
    Real Ql;                      //  Chemostat thermostat mass.  Units: (kcal fs^2 / mol)
    Real myEta;                   //  Nose-Hoover particle thermostat variable.  Units: (dimensionless)
    Real myEtaV;                  //  Nose-Hoover volume thermostat variable.  Units: (dimensionless)
    Real myEtaLambda;             //  Nose-Hoover chemostat thermostat variable.  Units: (dimensionless)
    Real myEtaVel;                //  Velocity of the thermostat variable.  Units: (fs)^-1
    Real myEtaVolVel;             //  Velocity of the volume thermostat variable.  Units: (fs)^-1
    Real myLambdaVel;             //  Velocity of the chemostat variable.  Units: (fs)^-1
    Real myEtaLambdaVel;           //  Velocity of the chemostat thermostat variable.  Units: (fs)^-1
    std::vector<int> N;           //  The number of molecules of each mixture component
    // ~~~ For tracking of the conserved quantity (CQ) ~~~
    Real AveCQ, AveCQSq;          //  The instantaneous CQ, average CQ, and CQ-squared
    Real AveDeltaMu;              //  Quantity needed for TI calculations
    //~~~ For computing the average chemostat temperature
    Real LambdaT;
    unsigned int TotSteps;

    int myTrajectoryFreq;         //  Frequency of transformation attempts to output the system trajectory
    int numAttempts;              //  Total # of transformation attempts since the last trajectory output
      
    //______________________________________________________________________osmoticType
    struct osmoticType{
      // This structure holds the data for the atom types belonging to a particular osmotic
      // molecule.  It consists of AtomTypes, where AtomType is defined above, a molecule type
      // name, the number of stages to use for insertion/deletion, and the target fugacity.

      // constructors
      osmoticType(){}
      osmoticType(Real F, unsigned int S): fugacity(F), NumberOfStages(S) {}

      /// the target fugacity for this molecule type
      Real fugacity;

      /// the number of stages to use for insertion/deletion
      unsigned int NumberOfStages;

      /// the lists of atoms, atomTypes, and bonding information for this molecule
      std::vector<Atom>            atoms;
      std::vector<STAGE::AtomType> atomTypes;
      std::vector<Bond>            bonds;
      std::vector<Angle>           angles;
      std::vector<Torsion>         dihedrals;
      std::vector<Torsion>         impropers;

      /// the information for one molecule of this component
      Molecule myMolecule;

      /// stock xyz coordinates for a molecule of this component
      Vector3DBlock myCoordinates;
    };

    /// osmoticType container
    std::vector<osmoticType> transformMaps;

    /// the number of osmotic molecule types in the system
    unsigned int numOsmComp;

    /// constant conversion factor needed to that the product
    /// f*V/(NkT) is dimensionless
    const Real fugacityFactor;
    
    /// seed for random number generator
    int seed;
      
  };
}
#endif

