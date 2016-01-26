/* -*- c++ -*- */
#ifndef EXTENDEDTIMEFORCE_H
#define EXTENDEDTIMEFORCE_H

#include "TimeForce.h"
#include "ExtendedForce.h"

namespace ProtoMol {

  //_________________________________________________________________ ExtendedTimeForce

  class ExtendedTimeForce : public TimeForce, public ExtendedForce {
    // This class contains the definition of one force
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ExtendedTimeForce(Force* actualForce);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ExtendedForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void evaluate(const GenericTopology* topo, 
			  const Vector3DBlock* positions,
			  const Vector3DBlock* velocities,
			  Vector3DBlock* forces, 
			  ScalarStructure* energies);

    virtual void parallelEvaluate(const GenericTopology* topo, 
				  const Vector3DBlock* positions, 
				  const Vector3DBlock* velocities,
				  Vector3DBlock* forces, 
				  ScalarStructure* energies);

  };

  //______________________________________________________________________ INLINES
}
#endif /* EXTENDEDTIMEFORCE_H */
