#include "registerForceExemplarsBonded.h"
#include "ForceFactory.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"
#include "BondSystemForce.h"
#include "AngleSystemForce.h"
#include "ImproperSystemForce.h"
#include "DihedralSystemForce.h"
#include "HarmDihedralSystemForce.h"

namespace ProtoMol
{
	void registerForceExemplarsBonded(const PeriodicBoundaryConditions*)
	{
		ForceFactory::registerExemplar(new BondSystemForce<PeriodicBoundaryConditions>());
		ForceFactory::registerExemplar(new AngleSystemForce<PeriodicBoundaryConditions>());
		ForceFactory::registerExemplar(new ImproperSystemForce<PeriodicBoundaryConditions>());
		ForceFactory::registerExemplar(new DihedralSystemForce<PeriodicBoundaryConditions>());
		ForceFactory::registerExemplar(new HarmDihedralSystemForce<PeriodicBoundaryConditions>());
	}

	void registerForceExemplarsBonded(const VacuumBoundaryConditions*)
	{
		ForceFactory::registerExemplar(new BondSystemForce<VacuumBoundaryConditions>());
		ForceFactory::registerExemplar(new AngleSystemForce<VacuumBoundaryConditions>());
		ForceFactory::registerExemplar(new ImproperSystemForce<VacuumBoundaryConditions>());
		ForceFactory::registerExemplar(new DihedralSystemForce<VacuumBoundaryConditions>());
		ForceFactory::registerExemplar(new HarmDihedralSystemForce<VacuumBoundaryConditions>());
	}
}
