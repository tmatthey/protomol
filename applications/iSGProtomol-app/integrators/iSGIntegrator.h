/*  -*- c++ -*-  */
#ifndef ISGINTEGRATOR_H
#define ISGINTEGRATOR_H

#include "iSGPAR.h"
#include "STSIntegrator.h"
#include "XSC.h"
#include "Array.h"
#include "TRANS.h"
#include "PSF.h"

namespace ProtoMol {

  //________________________________________________________ iSGIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;
  class Modifier;

  class iSGModifyForces;

  class iSGIntegrator: public STSIntegrator {

    friend class iSGModifyForces;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGIntegrator();  
    iSGIntegrator(Real timestep,
                  unsigned int numComp,
		              Real temperature,
                  Real pressure,
                  const std::vector<Real>& fugacityFrac,
                  Real tauT,    
		              Real tauV,    
                  Real tauP, 
		              Real tauD,
                  Real tauL,
		              Real MuTemp,
                  ForceGroup *overloadedForces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // new methods of class iSGIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //  Store the values of the extended system coordinates at
    //  the end of a simulation into an XSC object
    XSC & getXSC() const;
    Real getEpsilonVel() const {return myEpsilonVel;}
    Real getEtaVel() const {return myEtaVel;}
    Real getNumAtoms() const {return NumAtoms;}
    Real getAveDeltaMu(int NumSteps) {return (AveDeltaMu / static_cast<Real>(NumSteps));}
    Real getAveLambdaT(int NumSteps) {return (LambdaT / static_cast<Real>(NumSteps));}

    //  store the information for the transforming molecule topology
    void indexTopology(GenericTopology *topo, const TRANS &trans, const PSF& psf,
                       const iSGPAR& par, bool dihedralMultPSF, int theSeed);

    //  If an initial XSC file is specified, read in the values
    //  of the extended system coordinates
    virtual void readXSCs(const std::string myFile, GenericTopology *topo);

  private:
    void do2ndHalfKick();
    void PreForceThermostat(); 
    void PostForceThermostat(); 
    void PreForceBarostat(); 
    void PostForceBarostat();
    void PreForceChemostat();
    void PostForceChemostat();

    void modifyForces();
    void pickNewMolecule();
    void checkForCompletion();
    void updateTopology();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 11;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
                            Vector3DBlock   *positions, 
                            Vector3DBlock   *velocities, 
                            ScalarStructure *energies);
    virtual void run(int numTimesteps);

    /// Create a Rattle modifier
    virtual Modifier* createRattleModifier(Real eps, int maxIter);
    /// Create a Shake modifier 
    virtual Modifier* createShakeModifier(Real eps, int maxIter);

  protected:
    virtual void addModifierBeforeInitialize();
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator* doMake(std::string&, const std::vector<Value>& values, ForceGroup* fg)const;
  protected:
    virtual void doDrift();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doHalfKick();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:  
    const unsigned int myNumComp; //  Number of components in the mixture
    unsigned int myNumStages;     //  Number of stages that the transformation will be split into
    const Real myTargetTemp;      //  Target temperature.  Units: (K)
    const Real myTargetPres;      //  Target pressure.  Units: (bar)
    const Real myMuTemp;          //  Target chemostat temperature.  Units: (K)
    std::vector<Real> myTargetMu; //  Logarithms of the target fugacity fractions.  Units: (kcal/mol)
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
    unsigned int OldType;         //  The type of molecule being transformed
    unsigned int NewType;         //  The new identity we wish to give to the transformin molecule
    unsigned int FinalType;       //  The identity of the molecule at the end of the transformation attempt
    bool Transformed;             //  flag so we know when the molecule has been completely transformed
    Array<Real,3> myDeltaMuIG;    //  matrix of ideal gas chemical potential differences.  Units: (kcal/mol)
    std::vector<Real> myFugacityFrac;  //  Target fugacity fraction or chemical potential.  Units: (kcal/mol)
    Real myTargetDeltaMu;         //  residual chemical potential difference.  Units (kcal/mol)
    Real mylnMassRatio;           //  Ideal gas delta mu that accounts for the transformation of atomic mass
    Real myDMuIG;                 //  Ideal gas delta mu for the current stage.  Units (kcal/mol)
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
    
    int seed;                     //  seed for the random number generator

    //_________________________________________________________________ ChargeMap
    struct ChargeMap {
      // This class contains information common to one type of atom.  
  
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ChargeMap(){}
      // Default constructor for one atom type's set of charges for each transformation stage.
  
      // constructor for iSGMD simulations that creates memory for all identities
      ChargeMap(int i, int s) {
        old_charge.resize(ArraySizes(i)(i)(s));
        new_charge.resize(ArraySizes(i)(i)(s));
        alphaLJ.resize(ArraySizes(i)(i));
      }

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Array<Real,3> old_charge;
      Array<Real,3> new_charge;
      // This atomtype's old_type and new_type charges for each transformation stage.
      Array<Real,2> alphaLJ;
      // This atomtype's alphaLJ parameter for a switch between two identities
    };

    //_________________________________________________________________ local topology storage elements
    std::vector<iSGPAR::Bond> bonds;
    std::vector<iSGPAR::Angle> angles;
    std::vector<iSGPAR::Dihedral> dihedrals;
    std::vector<iSGPAR::Improper> impropers;
    std::vector<iSGPAR::AtomType> atomTypes;
    std::vector<ChargeMap> chargeMaps;
 
  };
}
#endif

