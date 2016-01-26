#include "registerTopologyExemplars.h"
#include "Topology.h"
#include "TopologyFactory.h"
#include "VacuumBoundaryConditions.h"
#include "CubicCellManager.h"
#include "PeriodicBoundaryConditions.h"
#include "Vector.h"

namespace ProtoMol {

  void registerTopologyExemplars(){
    // vacuum or normal boundary conditions
    TopologyFactory::registerExemplar(new Topology<VacuumBoundaryConditions,CubicCellManager>(),Vector<std::string>("NormalCubic"));
    // periodic boundary conditions
    TopologyFactory::registerExemplar(new Topology<PeriodicBoundaryConditions,CubicCellManager>());
  }
}
