#include "registerForceExemplarsOther.h"
#include "ForceFactory.h"
#include "CellListEnumerator_periodicBoundaries.h"
#include "CellListEnumerator_standard.h"
#include "PeriodicBoundaryConditions.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PaulTrapExtendedForce.h"
#include "ElectricFieldSystemForce.h"
#include "SphericalSystemForce.h"
#include "SphericalRestraintSystemForce.h"
#include "HapticSystemForce.h"
#include "FrictionExtendedForce.h"
#include "WrapperMetaForce.h"
#include "OneAtomPair.h"
#include "OneAtomPairTwo.h"
#include "OneMollyPair.h"
#include "OneMollyPairTwo.h"
#include "RangeSwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "C1SwitchingFunction.h"
#include "LennardJonesForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "NonbondedCutoffMollyForce.h"
#include "BondSystemForce.h"
#include "AngleSystemForce.h"
#include "oneAtomContraints.h"
#include "ExternalGravitationSystemForce.h"
#include "ExternalMagneticFieldExtendedForce.h"
#include "MagneticDipoleMirrorSystemForce.h"

namespace ProtoMol {
  void registerForceExemplarsOther(const PeriodicBoundaryConditions*, const CubicCellManager*){

//     ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
// 				   OneAtomPair<PeriodicBoundaryConditions,
// 				   RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce> >());
//
//     ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
// 				   OneAtomPair<PeriodicBoundaryConditions,
// 				   RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());
//
//     ForceFactory::registerExemplar(new NonbondedCutoffSystemForce<CubicCellManager,
// 				   OneAtomPairTwo<PeriodicBoundaryConditions,
// 				   RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
// 				   RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >());

    ForceFactory::registerExemplar(new WrapperMetaForce("MollyLennardJones",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPair<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce> >(),
							"MollyLennardJones",							
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPair<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							HBondConstraint> >(),
							"HBondLennardJones"));


    ForceFactory::registerExemplar(new WrapperMetaForce("MollyCoulomb",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPair<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >(),
							"MollyCoulomb",							
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPair<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce,
							HBondConstraint> >(),
							"HBondCoulomb"));


    ForceFactory::registerExemplar(new WrapperMetaForce("MollyLennardJonesCoulomb",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPairTwo<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >(),
							"MollyLennardJonesCoulomb",
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPairTwo<PeriodicBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce,
							HBondConstraint> >(),
							"HBondLennardJonesCoulomb"));

  }

  void registerForceExemplarsOther(const PeriodicBoundaryConditions*){

    ForceFactory::registerExemplar(new WrapperMetaForce("MollyBond",true,new BondSystemForce<PeriodicBoundaryConditions>(),"MollyBond"));

    ForceFactory::registerExemplar(new WrapperMetaForce("MollyAngle",true,new AngleSystemForce<PeriodicBoundaryConditions>(),"MollyAngle"));

    ForceFactory::registerExemplar(new ElectricFieldSystemForce<PeriodicBoundaryConditions,C1SwitchingFunction>());
    ForceFactory::registerExemplar(new MagneticDipoleMirrorSystemForce<PeriodicBoundaryConditions>());
  }

  void registerForceExemplarsOther(const VacuumBoundaryConditions*, const CubicCellManager*){ 

    ForceFactory::registerExemplar(new WrapperMetaForce("MollyLennardJones",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPair<VacuumBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce> >(),
							"MollyLennardJones",							
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPair<VacuumBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							HBondConstraint> >(),
							"HBondLennardJones"));


    ForceFactory::registerExemplar(new WrapperMetaForce("MollyCoulomb",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPair<VacuumBoundaryConditions,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >(),
							"MollyCoulomb",							
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPair<VacuumBoundaryConditions,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce,
							HBondConstraint> >(),
							"HBondCoulomb"));


    ForceFactory::registerExemplar(new WrapperMetaForce("MollyLennardJonesCoulomb",true, 
							new NonbondedCutoffMollyForce<CubicCellManager,
							OneMollyPairTwo<VacuumBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce> >(),
							"MollyLennardJonesCoulomb",
							new NonbondedCutoffSystemForce<CubicCellManager,
							OneAtomPairTwo<VacuumBoundaryConditions,
							RangeSwitchingFunction<C2SwitchingFunction>,LennardJonesForce,
							RangeSwitchingFunction<C1SwitchingFunction>,CoulombForce,
							HBondConstraint> >(),
							"HBondLennardJonesCoulomb"));




    ForceFactory::registerExemplar(new ElectricFieldSystemForce<VacuumBoundaryConditions,C1SwitchingFunction>());
  }

  void registerForceExemplarsOther(){ 
    ForceFactory::registerExemplar(new ExternalMagneticFieldExtendedForce());
    ForceFactory::registerExemplar(new ExternalGravitationSystemForce());
    ForceFactory::registerExemplar(new HapticSystemForce());
    ForceFactory::registerExemplar(new FrictionExtendedForce());
  }

  void registerForceExemplarsOther(const VacuumBoundaryConditions*){
    // PaulTrap for Coulomb Crystals
    ForceFactory::registerExemplar(new PaulTrapExtendedForce<VacuumBoundaryConditions>());

    // Spherical boundary conditions
    ForceFactory::registerExemplar(new SphericalSystemForce());
    ForceFactory::registerExemplar(new SphericalRestraintSystemForce());
    ForceFactory::registerExemplar(new WrapperMetaForce("MollyBond",true,new BondSystemForce<VacuumBoundaryConditions>(),"MollyBond"));

    ForceFactory::registerExemplar(new WrapperMetaForce("MollyAngle",true,new AngleSystemForce<VacuumBoundaryConditions>(),"MollyAngle"));
    ForceFactory::registerExemplar(new MagneticDipoleMirrorSystemForce<VacuumBoundaryConditions>());


  }



}
