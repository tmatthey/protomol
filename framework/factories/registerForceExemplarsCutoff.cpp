#include "registerForceExemplarsCutoff.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CnSwitchingFunction.h"
#include "CmpCnCnSwitchingFunction.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "CoulombForce.h"
#include "CoulombForceDiElec.h"
#include "CoulombSCPISMForce.h"
#include "CoulombBornRadiiForce.h"
#include "CoulombEwaldRealForce.h"
#include "CoulombEwaldRealTableForce.h"
#include "CoulombTableForce.h"
#include "CoulombMultiGridDirectForce.h"
#include "CoulombMultiGridDirectTableForce.h"
#include "CubicCellManager.h"
#include "CutoffSwitchingFunction.h"
#include "ForceFactory.h"
#include "LennardJonesForce.h"
#include "LennardJonesTableForce.h"
#include "MagneticDipoleForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "NonbondedCutoffBornForce.h"
#include "OneAtomPair.h"
#include "OneAtomPairTwo.h"
#include "PeriodicBoundaryConditions.h"
#include "ShiftSwitchingFunction.h"
#include "UniversalSwitchingFunction.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"



namespace ProtoMol {

  void registerForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*){
    // NonbondedCutoffSystemForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CnSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,CoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,CoulombForce> >());
    
    // NonbondedCutoffSystemForce LennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,ShiftSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<C2SwitchingFunction,7,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CnSwitchingFunction,7,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CmpCnCnSwitchingFunction,7,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce> >());

    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C1SwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,C2SwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CnSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CmpCnCnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<C2SwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<C2SwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CmpCnCnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<PeriodicBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CmpCnCnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());

    // MagneticDipole 
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<PeriodicBoundaryConditions,C1SwitchingFunction,MagneticDipoleForce> >());
  }
  void registerForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*){

    // NonbondedCutoffSystemForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForce> >());

    // NonbondedCutoffSystemForce CoulombForceDiElec
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombForceDiElec> >());

    
    // NonbondedCutoffSystemForce LennardJonesForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,LennardJonesForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<C2SwitchingFunction,7,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CnSwitchingFunction,7,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CmpCnCnSwitchingFunction,7,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesForce> >());
    
    // NonbondedCutoffSystemForce LennardJonesForce CoulombForce
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    // ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());

    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C1> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C2> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C3> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectForce<CoulombForce::C4> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C1,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C2,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C3,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombMultiGridDirectTableForce<CoulombForce::C4,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    
	ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombEwaldRealForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CutoffSwitchingFunction,LennardJonesForce,UniversalSwitchingFunction,CoulombEwaldRealTableForce<CutoffSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<C2SwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,UniversalSwitchingFunction,LennardJonesTableForce<CmpCnCnSwitchingFunction,7,Real>,UniversalSwitchingFunction,CoulombTableForce<ShiftSwitchingFunction,6,Real> > >());

    // NonbondedCutoffSystemForce LennardJonesForce CoulombForceDiElec
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C2SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombForceDiElec> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombForceDiElec> >());

    // MagneticDipole
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,MagneticDipoleForce> >());


    // SCPISM stuff
    // NonbondedCutoffSystemForce CoulombSCPISMForce
     ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C1SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,C2SwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,C2SwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CnSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,C1SwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,CutoffSwitchingFunction,CoulombSCPISMForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPairTwo<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,LennardJonesForce,ShiftSwitchingFunction,CoulombSCPISMForce> >());

    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombBornRadiiForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombBornRadiiForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombBornRadiiForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombBornRadiiForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombBornRadiiForce> >());
    ForceFactory::registerExemplar(new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombBornRadiiForce> >());



  }
}
