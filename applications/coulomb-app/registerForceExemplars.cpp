#include "registerForceExemplars.h"
#include "Parallel.h"
#include "Topology.h"
#include "GenericTopology.h"
#include "ForceFactory.h"

#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"

#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "CubicCellManager.h"
#include "CoulombForce.h"
#include "OneAtomPair.h"
#include "UniversalSwitchingFunction.h"

#include "NonbondedMultiGridSystemForce.h"
#include "PaulTrapExtendedForce.h"
#include "NonbondedSimpleFullSystemForce.h"



#include "BSpline.h"
#include "Lagrange.h"
#include "Hermite.h"

#include "Vector.h"

namespace ProtoMol {

  void registerForceExemplars(const GenericTopology* topo){
    if(dynamic_cast<const SemiGenericTopology<VacuumBoundaryConditions>* >(topo) != NULL){
      ForceFactory::registerExemplar(new PaulTrapExtendedForce<VacuumBoundaryConditions>());
      ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForce> >());
      ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true,true>());
    }
  }
}
