#include "iSGregisterForceExemplarsFastElectrostatic.h"
#include "ForceFactory.h"
#include "iSGNonbondedFullEwaldSystemForce.h"
#include "iSGNonbondedPMEwaldSystemForce.h"
//#include "iSGNonbondedMultiGridSystemForce.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "CubicCellManager.h"
#include "iSGCoulombForce.h"
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
  void iSGregisterForceExemplarsFastElectrostatic(const PeriodicBoundaryConditions*, const CubicCellManager*) {
    // iSGFullEwald
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,CutoffSwitchingFunction>(),Vector<std::string>("iSGCoulombEwald"));
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#endif

    // iSGPME
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>(),Vector<std::string>("iSGCoulombPME"));
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());

    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<PeriodicBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>());
#endif


/*
    // MultiGrid
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,false>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<PeriodicBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,false>());
#endif
*/
  }

  void iSGregisterForceExemplarsFastElectrostatic(const VacuumBoundaryConditions*, const CubicCellManager*) { 
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,CutoffSwitchingFunction>(),Vector<std::string>("iSGCoulombEwald"));
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedFullEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,CutoffSwitchingFunction>());
#endif


    // PME
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,CutoffSwitchingFunction>(),Vector<std::string>("iSGCoulombPME"));
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,ShiftSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,BSpline,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,BSpline,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,BSpline,CutoffSwitchingFunction>());

    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,true,false,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,true,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new iSGNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,true,false,Hermite,CutoffSwitchingFunction>());
    ForceFactory::registerExemplar(new vNonbondedPMEwaldSystemForce<VacuumBoundaryConditions,CubicCellManager,false,false,true,Hermite,CutoffSwitchingFunction>());
#endif


/*
    // MultiGrid
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C1,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C2,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C3,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,Hermite,CoulombForce::C4,true,false>());
#ifdef USE_FULL_FORCE_SET
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C1,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C2,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C3,true,false>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,false,true>());
    ForceFactory::registerExemplar(new NonbondedMultiGridSystemForce<VacuumBoundaryConditions,CubicCellManager,BSpline,CoulombForce::C4,true,false>());
#endif
*/
  }
}
