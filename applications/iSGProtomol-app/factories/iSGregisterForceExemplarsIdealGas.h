/* -*- c++ -*- */
#ifndef ISGREGISTERFORCEEXEMPLARSIDEALGAS_H
#define ISGREGISTERFORCEEXEMPLARSIDEALGAS_H


namespace ProtoMol {

  class CubicCellManager;
  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void iSGregisterForceExemplarsIdealGas(const PeriodicBoundaryConditions*, const CubicCellManager*);
  void iSGregisterForceExemplarsIdealGas(const VacuumBoundaryConditions*, const CubicCellManager*);

}
#endif /* ISGREGISTERFORCEEXEMPLARSIDEALGAS_H */
