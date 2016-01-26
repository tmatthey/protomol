#include "ExtendedTimeForce.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"

namespace ProtoMol {

  //_________________________________________________________________ ExtendedTimeForce

  ExtendedTimeForce::ExtendedTimeForce(Force* actualForce):TimeForce(actualForce){}

  void ExtendedTimeForce::evaluate(const GenericTopology* topo, 
				   const Vector3DBlock* positions,
				   const Vector3DBlock* velocities,
				   Vector3DBlock* forces, 
				   ScalarStructure* energies){
    preprocess(positions->size());
    (dynamic_cast<ExtendedForce*>(myActualForce))->evaluate(topo,
							    positions,
							    velocities,
							    forces,
							    energies);
    postprocess(topo,forces,energies);
  }

  void ExtendedTimeForce::parallelEvaluate(const GenericTopology* topo, 
					   const Vector3DBlock* positions, 
					   const Vector3DBlock* velocities,
					   Vector3DBlock* forces, 
					   ScalarStructure* energies){
    preprocess(positions->size());
    (dynamic_cast<ExtendedForce*>(myActualForce))->parallelEvaluate(topo,
								    positions,
								    velocities,
								    forces,
								    energies);
    postprocess(topo,forces,energies);
  }

}
