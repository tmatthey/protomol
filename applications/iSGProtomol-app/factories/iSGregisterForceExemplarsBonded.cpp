#include "iSGregisterForceExemplarsBonded.h"
#include "ForceFactory.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"
#include "iSGBondSystemForce.h"
#include "iSGAngleSystemForce.h"
#include "iSGImproperSystemForce.h"
#include "iSGDihedralSystemForce.h"
namespace ProtoMol {

  void iSGregisterForceExemplarsBonded(const PeriodicBoundaryConditions*){
    ForceFactory::registerExemplar(new iSGBondSystemForce<PeriodicBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGAngleSystemForce<PeriodicBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGImproperSystemForce<PeriodicBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGDihedralSystemForce<PeriodicBoundaryConditions>());
  }
  void iSGregisterForceExemplarsBonded(const VacuumBoundaryConditions*){
    ForceFactory::registerExemplar(new iSGBondSystemForce<VacuumBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGAngleSystemForce<VacuumBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGImproperSystemForce<VacuumBoundaryConditions>());
    ForceFactory::registerExemplar(new iSGDihedralSystemForce<VacuumBoundaryConditions>());
  }
}
