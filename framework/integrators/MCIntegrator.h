/*  -*- c++ -*-  */
#ifndef MCINTEGRATOR_H
#define MCINTEGRATOR_H

#include "MTSIntegrator.h"

namespace ProtoMol {
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;

  //____________________________________________________________ MCIntegrator
  /**

  Base class for Monte Carlo MTS integrators.

  */
  class MCIntegrator: public MTSIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MCIntegrator();  
    MCIntegrator(int cycles,
		 Real initialTemperature,
		 ForceGroup *overloadedForces, 
		 StandardIntegrator *nextIntegrator);

    virtual ~MCIntegrator();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void getParameters(std::vector<Parameter>& parameters) const;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
			    Vector3DBlock *positions,
			    Vector3DBlock *velocities,
			    ScalarStructure *energies);

    virtual void run(int numTimesteps);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods for class MCIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Real getInitialTemperature() const { return myInitialTemperature; }



  protected:

    /// Implements the MC algorithm
    virtual void walk(int steps);

    /// Metropolis test interface which calls metropolisTest with arguments.
    virtual bool metropolisTest();

    /// Implementation of the metropolis test
    bool metropolisTest( Real newEnergy, Real oldEnergy );

    /// Implementation of the metropolis test that returns the probability.
    bool metropolisTest( Real newEnergy, Real oldEnergy, Real & prob );

    ///  Each class must implement its own system perturbation method.
    virtual void perturbSystem() = 0;

    ///  Save/restore state implementations (default pos, vel & energies)
    virtual void saveValues();
    virtual void restoreValues();


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

    Real myInitialTemperature;

  protected:

    Real myOldKineticEnergy;

    Vector3DBlock*   myOldPositions;
    Vector3DBlock*   myOldVelocities;

    ScalarStructure* myOldEnergies;

  };

}
#endif
