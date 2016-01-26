#include "registerTopologyExemplars.h"
#include "Topology.h"
#include "TopologyFactory.h"
#include "VacuumBoundaryConditions.h"
#include "CubicCellManager.h"
#include "Vector.h"

namespace ProtoMol {

  void registerTopologyExemplars(){
    TopologyFactory::registerExemplar(new Topology<VacuumBoundaryConditions,CubicCellManager>(),Vector<std::string>("NormalCubic"));
  }
}
