#include "oSGregisterForceExemplarsSimpleFull.h"
#include "ScalarStructure.h"
#include "Topology.h"
#include "oSGOneAtomPair.h"
#include "oSGOneAtomPairTwo.h"
//#include "oSGCoulombForce.h"
#include "ForceFactory.h"
#include "oSGLennardJonesForce.h"
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
  void oSGregisterForceExemplarsSimpleFull(const PeriodicBoundaryConditions*){

/*
    // NonbondedSimpleFullSystemForce ISGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,oSGCoulombForce> >());
*/
    
    // NonbondedSimpleFullSystemForce ISGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,oSGLennardJonesForce> >());
/*
    // NonbondedSimpleFullSystemForce oSGLennardJonesForce oSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,oSGLennardJonesForce,UniversalSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGLennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,oSGCoulombForce> >());
*/                                       
  }

  void oSGregisterForceExemplarsSimpleFull(const VacuumBoundaryConditions*){
/*
   // NonbondedSimpleFullSystemForce ISGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<oSGOneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,oSGCoulombForce> >());
*/

    // NonbondedSimpleFullSystemForce ISGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   oSGOneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,oSGLennardJonesForce> >());
/*                                   
    // NonbondedSimpleFullSystemForce oSGLennardJonesForce oSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,oSGLennardJonesForce,UniversalSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,oSGLennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,oSGCoulombForce> >());
*/  
  }
}
