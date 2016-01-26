#include "oSGregisterForceExemplarsCutoff.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
//#include "oSGCoulombForce.h"
//#include "oSGCoulombEwaldRealForce.h"
//#include "oSGCoulombEwaldRealTableForce.h"
#include "CubicCellManager.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "oSGLennardJonesForce.h"
//#include "oSGLennardJonesTableForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "oSGOneAtomPair.h"
#include "oSGOneAtomPairTwo.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"

namespace ProtoMol {
  void oSGregisterForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*){
/*
    // NonbondedCutoffSystemForce oSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,oSGCoulombForce> >());
*/
    // NonbondedCutoffSystemForce oSGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,oSGLennardJonesForce> >());
/*
    // NonbondedCutoffSystemForce oSGLennardJonesForce oSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,C1SwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,ShiftSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,ShiftSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >());    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >()); 
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());
                                  
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >());
*/
    
  }
  void oSGregisterForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*){
/* 
    // NonbondedCutoffSystemForce oSGCoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,oSGCoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,oSGCoulombForce> >());
*/
    // NonbondedCutoffSystemForce oSGLennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    				   oSGOneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
    			           oSGOneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
				   oSGOneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,oSGLennardJonesForce> >());
/*
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,C1SwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,ShiftSwitchingFunction,oSGCoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,ShiftSwitchingFunction,oSGCoulombForce> >());
    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >());   
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >());    
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());

    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealTableForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
                                   oSGOneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,oSGLennardJonesForce,CutoffSwitchingFunction,oSGCoulombEwaldRealForce> >());
*/
  }
}
