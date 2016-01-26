#include "iSGregisterForceExemplars.h"
#include "iSGregisterForceExemplarsBonded.h"
#include "iSGregisterForceExemplarsCutoff.h"
#include "iSGregisterForceExemplarsSimpleFull.h"
#include "iSGregisterForceExemplarsFull.h"
#include "iSGregisterForceExemplarsFastElectrostatic.h"
#include "iSGregisterForceExemplarsIdealGas.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"
#include "CubicCellManager.h"

namespace ProtoMol {

  template<typename BC,typename CM>
  inline void iSGregisterForceExemplarsDispatch(const Topology<BC,CM>* topo){
    if(topo == NULL)
      return;
    const BC* bc=NULL;
    const CM* cm=NULL;

    iSGregisterForceExemplarsCutoff(bc,cm);

    iSGregisterForceExemplarsFull(bc);

    iSGregisterForceExemplarsSimpleFull(bc);

    iSGregisterForceExemplarsBonded(bc);

    iSGregisterForceExemplarsFastElectrostatic(bc,cm);

    iSGregisterForceExemplarsIdealGas(bc,cm);
  }

  void iSGregisterForceExemplars(const GenericTopology* topo){
    iSGregisterForceExemplarsDispatch(dynamic_cast<const Topology<PeriodicBoundaryConditions,CubicCellManager>*>(topo));
    iSGregisterForceExemplarsDispatch(dynamic_cast<const Topology<VacuumBoundaryConditions,CubicCellManager>*>(topo));
  }
}
