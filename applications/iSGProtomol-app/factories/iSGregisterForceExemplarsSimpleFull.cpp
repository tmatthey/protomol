#include "iSGregisterForceExemplarsSimpleFull.h"
#include "ScalarStructure.h"
#include "Topology.h"
#include "iSGOneAtomPair.h"
#include "iSGOneAtomPairTwo.h"
#include "iSGCoulombForce.h"
#include "ForceFactory.h"
#include "iSGLennardJonesForce.h"
#include "NonbondedSimpleFullSystemForce.h"
#include "PeriodicBoundaryConditions.h"
#include "UniversalSwitchingFunction.h"
#include "ComplementSwitchingFunction.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CutoffSwitchingFunction.h"
#include "ShiftSwitchingFunction.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void iSGregisterForceExemplarsSimpleFull(const PeriodicBoundaryConditions*){

    // NonbondedSimpleFullSystemForce ISGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,iSGCoulombForce> >());
    
    // NonbondedSimpleFullSystemForce ISGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,iSGLennardJonesForce> >());

    // NonbondedSimpleFullSystemForce iSGLennardJonesForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,iSGLennardJonesForce,UniversalSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGLennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,iSGCoulombForce> >());
                                       
  }

  void iSGregisterForceExemplarsSimpleFull(const VacuumBoundaryConditions*){

   // NonbondedSimpleFullSystemForce ISGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<iSGOneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,iSGCoulombForce> >());

    // NonbondedSimpleFullSystemForce ISGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   iSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,iSGLennardJonesForce> >());
                                   
    // NonbondedSimpleFullSystemForce iSGLennardJonesForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,iSGLennardJonesForce,UniversalSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,iSGLennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,iSGCoulombForce> >());
  }
}
