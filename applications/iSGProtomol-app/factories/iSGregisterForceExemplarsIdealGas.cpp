#include "iSGregisterForceExemplarsIdealGas.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "iSGIdealGasCoulombForce.h"
#include "CubicCellManager.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "iSGIdealGasLennardJonesForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "iSGOneAtomPair.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void iSGregisterForceExemplarsIdealGas(const PeriodicBoundaryConditions*, const CubicCellManager*){

    // NonbondedCutoffSystemForce iSGIdealGasCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGIdealGasCoulombForce> >());

    // NonbondedCutoffSystemForce ISGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGIdealGasLennardJonesForce> >());
    
  }
  void iSGregisterForceExemplarsIdealGas(const VacuumBoundaryConditions*, const CubicCellManager*){
 
    // NonbondedCutoffSystemForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGIdealGasCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,iSGIdealGasCoulombForce> >());

    // NonbondedCutoffSystemForce iSGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    				   iSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    			           iSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGIdealGasLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,iSGIdealGasLennardJonesForce> >());

  }
}
