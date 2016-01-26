#include "iSGregisterForceExemplarsCutoff.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "iSGCoulombForce.h"
#include "iSGCoulombEwaldRealForce.h"
#include "iSGCoulombEwaldRealTableForce.h"
#include "CubicCellManager.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "iSGLennardJonesForce.h"
#include "iSGLennardJonesTableForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "iSGOneAtomPair.h"
#include "iSGOneAtomPairTwo.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void iSGregisterForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*){

    // NonbondedCutoffSystemForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGCoulombForce> >());

    // NonbondedCutoffSystemForce iSGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,iSGLennardJonesForce> >());

    // NonbondedCutoffSystemForce iSGLennardJonesForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,C1SwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,ShiftSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,ShiftSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >());    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >()); 
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());
                                  
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >());

    
  }
  void iSGregisterForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*){
 
    // NonbondedCutoffSystemForce iSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,iSGCoulombForce> >());

    // NonbondedCutoffSystemForce iSGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    				   iSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    			           iSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   iSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,iSGLennardJonesForce> >());

    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,C1SwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,ShiftSwitchingFunction,iSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,ShiftSwitchingFunction,iSGCoulombForce> >());
    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >());   
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >());    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());

    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   iSGOneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,iSGLennardJonesForce,CutoffSwitchingFunction,iSGCoulombEwaldRealForce> >());

  }
}
