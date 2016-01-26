#include "CheckpointOutputStream.h"
#include "GenericTopology.h"
#include "SemiGenericTopology.h"
#include "PeriodicBoundaryConditions.h"
#include "VacuumBoundaryConditions.h"
#include "Integrator.h"

namespace ProtoMol {

void CheckpointOutputStream::writeArchive(GenericTopology* topo,
					  Vector3DBlock positions,
					  Vector3DBlock velocities,
					  Integrator* integrator,
					  ScalarStructure scalar)
{
      unsigned int syssize = positions.size();
      (*this) << syssize;
      unsigned int levels = integrator->level(); // THIS IS THE HIGHEST INTEGRATOR
      (*this) << levels;
      Real ts = topo->time;
      //Real ts = integrator->getTimestep();
      (*this) << ts;
      unsigned short* cRN;
      unsigned short junk[3] = {1, 1, 1};
      cRN = seed48(junk); // Get the last random number
      unsigned short cRN2[3] = {cRN[0], cRN[1], cRN[2]};
      (*this) << cRN2[0] << cRN2[1] << cRN2[2];
      seed48(cRN2); // Restore random number generator state
      if (((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.getKeyword() == "Periodic") {
	(*this) << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().z;
	(*this) << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().z;
	(*this) << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().z;
	(*this) << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().z;
      }
      else {
	Real zero = 0.0;
	(*this) << zero << zero << zero;
	(*this) << zero << zero << zero;
	(*this) << zero << zero << zero;
	(*this) << zero << zero << zero;
      }
      for (unsigned int i = 0; i < positions.size(); i++) {
	(*this) << positions[i].x << positions[i].y << positions[i].z;
	(*this) << velocities[i].x << velocities[i].y << velocities[i].z;
	Integrator* tmp = integrator;
	for (unsigned int j = 0; j <= levels; j++) {
	  (*this) << (*(tmp->getForces()))[i].x << (*(tmp->getForces()))[i].y << (*(tmp->getForces()))[i].z;
	  if (j != levels)
	    tmp = tmp->next();
	}
      }
      unsigned int first = scalar.FIRST;
      unsigned int last = scalar.LAST;
      for (unsigned int i = first; i < last; i++) {
	Real val = scalar[(ScalarStructure::Index)i];
	(*this) << val;
      }

}

};

