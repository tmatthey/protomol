/* -*- c++ -*- */
#ifndef SYSTEMCOMPAREFORCE_H
#define SYSTEMCOMPAREFORCE_H

#include "CompareForce.h"
#include "SystemForce.h"
#include "ReducedHessAngle.h"

namespace ProtoMol {

  //_________________________________________________________________ SystemCompareForce

  class SystemCompareForce : public CompareForce, public SystemForce {
    // This class contains the definition of one force
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SystemCompareForce(Force* actualForce, CompareForce* compareForce);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void evaluate(const GenericTopology* topo, 
			  const Vector3DBlock* positions,
			  Vector3DBlock* forces, 
			  ScalarStructure* energies);

    virtual void parallelEvaluate(const GenericTopology* topo, 
				  const Vector3DBlock* positions, 
				  Vector3DBlock* forces, 
				  ScalarStructure* energies);

  };

  //______________________________________________________________________ INLINES
}
#endif /* SYSTEMCOMPAREFORCE_H */
