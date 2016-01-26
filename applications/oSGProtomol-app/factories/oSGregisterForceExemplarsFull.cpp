#include "oSGregisterForceExemplarsFull.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CubicCellManager.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
//#include "oSGCoulombForce.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "oSGLennardJonesForce.h"
#include "NonbondedFullSystemForce.h"
#include "oSGOneAtomPairFull.h"
#include "oSGOneAtomPairTwoFull.h"
#include "Parameter.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "UniversalSwitchingFunction.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void oSGregisterForceExemplarsFull(const PeriodicBoundaryConditions*){

    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,oSGCoulombForce> >());
    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,oSGCoulombForce> >());
    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGCoulombForce> >());
    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairFull<PeriodicBoundaryConditions,ShiftSwitchingFunction,oSGLennardJonesForce> >());


    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairTwoFull<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,ShiftSwitchingFunction,oSGCoulombForce> >());
    //ForceFactory::registerExemplar(new NonbondedFullSystemForce<oSGOneAtomPairTwoFull<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombForce> >());

  }
   
  void oSGregisterForceExemplarsFull(const VacuumBoundaryConditions*){
    // Huh, are you really sure you'll like to add such forces for Vacumm?
  }
}
