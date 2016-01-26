/*  -*- c++ -*-  */
#ifndef NFPTINTEGRATOR_H
#define NFPTINTEGRATOR_H

#include "oSGIntegrator.h"
#include "XSC.h"

namespace ProtoMol {

  //________________________________________________________ NfPTIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;
  class Modifier;

  class oSGModifierPreForceChemostat;
  class oSGModifierPostForceChemostat;

  class NfPTIntegrator: public oSGIntegrator {

    friend class oSGModifierPreForceChemostat;
    friend class oSGModifierPostForceChemostat;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NfPTIntegrator();
    NfPTIntegrator(Real timestep,
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
    // new methods of class NfPTIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //  Store the values of the extended system coordinates at
    //  the end of a simulation into an XSC object
    XSC & getXSC() const;
        
  private:
    virtual void do2ndHalfKick();
    virtual void PreForceChemostat();
    virtual void PostForceChemostat();
    virtual void deletionDeltaMu(const Real myFugacity, const int Osm) {
      myTargetMu = -myFugacity + kbT * log (static_cast<Real>(N[Osm]));
      myTargetMu /= static_cast<Real>(myNumStages);
    }
    virtual void insertionDeltaMu(const Real myFugacity, const int Osm) { 
      myTargetMu = myFugacity - kbT * log (static_cast<Real>(N[Osm]+1));
      myTargetMu /= static_cast<Real>(myNumStages);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 10;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
                            Vector3DBlock   *positions,
                            Vector3DBlock   *velocities,
                            ScalarStructure *energies);

    /// Create a Rattle modifier
    virtual Modifier* createRattleModifier(Real eps, int maxIter);
    /// Create a Shake modifier
    virtual Modifier* createShakeModifier(Real eps, int maxIter);


  protected:
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

  };
}
#endif

