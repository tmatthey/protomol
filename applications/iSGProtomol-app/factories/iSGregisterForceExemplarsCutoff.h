/* -*- c++ -*- */
#ifndef ISGREGISTERFORCEEXEMPLARSCUTOFF_H
#define ISGREGISTERFORCEEXEMPLARSCUTOFF_H


namespace ProtoMol {

  class CubicCellManager;
  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void iSGregisterForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*);
  void iSGregisterForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*);

}
#endif /* REGISTERFORCEEXEMPLARS_H */
