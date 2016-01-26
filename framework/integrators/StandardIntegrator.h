/*  -*- c++ -*-  */
#ifndef STANDARDINTEGRATOR_H
#define STANDARDINTEGRATOR_H

#include "Integrator.h"

namespace ProtoMol {
  //_________________________________________________________________ StandardIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;

  class StandardIntegrator: public Integrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    friend class MTSIntegrator;
    friend class S2HMCIntegrator;
  public:
    StandardIntegrator();
    StandardIntegrator(ForceGroup* forceGroup);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doHalfKick();
    virtual void doKick();
    virtual void initializeForces();
    virtual void doDriftOrNextIntegrator() = 0;
    virtual void calculateForces();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void run (int numTimesteps);
    virtual void initialize(GenericTopology *topo, 
			    Vector3DBlock   *positions, 
			    Vector3DBlock   *velocities, 
			    ScalarStructure *energies);

    virtual Integrator* previous();
    virtual const Integrator* previous() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    StandardIntegrator *myPreviousIntegrator;

  };
  //______________________________________________________________________ INLINES

}
#endif
