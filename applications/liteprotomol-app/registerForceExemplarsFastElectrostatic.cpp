#include "registerForceExemplarsFastElectrostatic.h"
#include "ForceFactory.h"
#include "NonbondedFullEwaldSystemForce.h"
#include "NonbondedPMEwaldSystemForce.h"
#include "NonbondedMultiGridSystemForce.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "CubicCellManager.h"
#include "CoulombForce.h"
#include "CutoffSwitchingFunction.h"
#include "ShiftSwitchingFunction.h"
#include "C1SwitchingFunction.h"

#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"

#include "CubicCellManager.h"

#include "BSpline.h"
#include "Lagrange.h"
#include "Hermite.h"

#include "Vector.h"

namespace ProtoMol {
  void registerForceExemplarsFastElectrostatic(const PeriodicBoundaryConditions*, const CubicCellManager*){
#ifndef LITE_PERIODIC
    // Full Ewald
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,CutoffSwitchingFunction>(),Vector<std::string>("CoulombEwald"));
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#endif

    // PME
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>(),Vector<std::string>("CoulombPME"));
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());

    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>());
#endif

    // MultiGrid
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,true,true>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,false,true,true>());
#endif
#endif
  }

  void registerForceExemplarsFastElectrostatic(const VacuumBoundaryConditions*, const CubicCellManager*){ 
#ifndef LITE_SWITCH
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,CutoffSwitchingFunction>(),Vector<std::string>("CoulombEwald"));
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#endif
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,C1SwitchingFunction>());
#ifndef LITE_SWITCH
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#endif
#endif


#ifndef LITE_SWITCH
    // PME
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>(),Vector<std::string>("CoulombPME"));
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());
#endif
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>());
#ifndef LITE_SWITCH
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());

    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new NonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>());
#endif
#endif

#ifndef LITE_SWITCH
    // MultiGrid
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true,true>());
#endif
#ifndef LITE_SWITCH
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,true,true>());
#endif
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,false,true,true>());
#ifndef LITE_SWITCH
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,false,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,false,true,true>());
#endif
#endif
  }

}
