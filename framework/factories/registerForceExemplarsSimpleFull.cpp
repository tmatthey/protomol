#include "registerForceExemplarsSimpleFull.h"
#include "ScalarStructure.h"
#include "Topology.h"
#include "OneAtomPair.h"
#include "OneAtomPairTwo.h"
#include "CoulombForce.h"
#include "CoulombForceDiElec.h"
#include "CoulombSCPISMForce.h"
#include "ForceFactory.h"
#include "GravitationForce.h"
#include "LennardJonesForce.h"
#include "MagneticDipoleForce.h"
#include "NonbondedSimpleFullSystemForce.h"
#include "PeriodicBoundaryConditions.h"
#include "UniversalSwitchingFunction.h"
#include "ComplementSwitchingFunction.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CnSwitchingFunction.h"
#include "CutoffSwitchingFunction.h"
#include "ShiftSwitchingFunction.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void registerForceExemplarsSimpleFull(const PeriodicBoundaryConditions*){
    // NonbondedSimpleFullSystemForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,CoulombForce> >());
    
    // NonbondedSimpleFullSystemForce LennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,LennardJonesForce> >());
    
    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<PeriodicBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<PeriodicBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombForce> >());

    //CoulombSCPISM


    // MagneticDipole 
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,MagneticDipoleForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,MagneticDipoleForce> >());
    // GravitationForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,GravitationForce> >());

  }

  void registerForceExemplarsSimpleFull(const VacuumBoundaryConditions*){
   // NonbondedSimpleFullSystemForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,CoulombForce> >());

   // NonbondedSimpleFullSystemForce CoulombForceDiElec
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,CoulombForceDiElec> >());

    //CoulombSCPISM
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombSCPISMForce> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<CutoffSwitchingFunction>,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<ShiftSwitchingFunction>,LennardJonesForce> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());

    // NonbondedSimpleFullSystemForce LennardJonesForce CoulombForceDiElec
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,ComplementSwitchingFunction<CnSwitchingFunction>,LennardJonesForce,ComplementSwitchingFunction<C1SwitchingFunction>,CoulombForceDiElec> >());

    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());


    // MagneticDipole 
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,MagneticDipoleForce> >());
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,ComplementSwitchingFunction<C1SwitchingFunction>,MagneticDipoleForce> >());
    // GravitationForce
    ForceFactory::registerExemplar(new NonbondedSimpleFullSystemForce<
				   OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,GravitationForce> >());

  }
}
