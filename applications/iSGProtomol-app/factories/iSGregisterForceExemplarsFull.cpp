#include "iSGregisterForceExemplarsFull.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CubicCellManager.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "iSGCoulombForce.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "iSGLennardJonesForce.h"
#include "NonbondedFullSystemForce.h"
#include "iSGOneAtomPairFull.h"
#include "iSGOneAtomPairTwoFull.h"
#include "Parameter.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "UniversalSwitchingFunction.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void iSGregisterForceExemplarsFull(const PeriodicBoundaryConditions*){

    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGLennardJonesForce> >());


    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairTwoFull<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,ShiftSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<iSGOneAtomPairTwoFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombForce> >());

  }
   
  void iSGregisterForceExemplarsFull(const VacuumBoundaryConditions*){
    // Huh, are you really sure you'll like to add such forces for Vacumm?
  }
}
