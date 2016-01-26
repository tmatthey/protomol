#include "oSGregisterForceExemplars.h"
#include "registerForceExemplarsBonded.h"
#include "oSGregisterForceExemplarsCutoff.h"
#include "oSGregisterForceExemplarsSimpleFull.h"
#include "oSGregisterForceExemplarsFull.h"
//#include "oSGregisterForceExemplarsFastElectrostatic.h"
#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "PeriodicBoundaryConditions.h"
#include "CubicCellManager.h"

namespace ProtoMol {

  template<typename BC,typename CM>
  inline void oSGregisterForceExemplarsDispatch(const Topology<BC,CM>* topo){
    if(topo == NULL)
      return;
    const BC* bc=NULL;
    const CM* cm=NULL;

    oSGregisterForceExemplarsCutoff(bc,cm);

    oSGregisterForceExemplarsFull(bc);

    oSGregisterForceExemplarsSimpleFull(bc);

    registerForceExemplarsBonded(bc);

    //oSGregisterForceExemplarsFastElectrostatic(bc,cm);
  }

  void oSGregisterForceExemplars(const GenericTopology* topo){
    oSGregisterForceExemplarsDispatch(dynamic_cast<const Topology<PeriodicBoundaryConditions,CubicCellManager>*>(topo));
    oSGregisterForceExemplarsDispatch(dynamic_cast<const Topology<VacuumBoundaryConditions,CubicCellManager>*>(topo));
  }
}
