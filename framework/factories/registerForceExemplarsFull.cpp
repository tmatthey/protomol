#include "registerForceExemplarsFull.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CnSwitchingFunction.h"
#include "CmpCnCnSwitchingFunction.h"
#include "CubicCellManager.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "CoulombForce.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "LennardJonesForce.h"
#include "NonbondedFullSystemForce.h"
#include "OneAtomPairFull.h"
#include "OneAtomPairTwoFull.h"
#include "Parameter.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "UniversalSwitchingFunction.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void registerForceExemplarsFull(const PeriodicBoundaryConditions*){
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,LennardJonesForce> >());


    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairTwoFull<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<OneAtomPairTwoFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
#ifdef USE_FULL_FORCE_SET
#endif
  }
  void registerForceExemplarsFull(const VacuumBoundaryConditions*){
    // Huh, are you really sure you'll like to add such forces for Vacumm?
  }
}
