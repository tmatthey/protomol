#include "registerForceExemplars.h"
#include "registerForceExemplarsBonded.h"
#include "registerForceExemplarsCutoff.h"
#include "registerForceExemplarsOther.h"
#include "registerForceExemplarsSimpleFull.h"
#include "registerForceExemplarsFull.h"
#include "registerForceExemplarsFastElectrostatic.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"
#include "CubicCellManager.h"

namespace ProtoMol
{
	template <typename BC, typename CM>
	inline void registerForceExemplarsDispatch(const Topology<BC, CM>* topo)
	{
		if (topo == NULL)
			return;
		const BC* bc = NULL;
		const CM* cm = NULL;

		registerForceExemplarsCutoff(bc, cm);

		registerForceExemplarsFull(bc);

		registerForceExemplarsSimpleFull(bc);

		registerForceExemplarsBonded(bc);

		registerForceExemplarsOther();
		registerForceExemplarsOther(bc);
		registerForceExemplarsOther(bc, cm);

		registerForceExemplarsFastElectrostatic(bc, cm);
	}

	void registerForceExemplars(const GenericTopology* topo)
	{
		registerForceExemplarsDispatch(dynamic_cast<const Topology<PeriodicBoundaryConditions, CubicCellManager>*>(topo));
		registerForceExemplarsDispatch(dynamic_cast<const Topology<VacuumBoundaryConditions, CubicCellManager>*>(topo));
	}
}
