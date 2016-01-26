/* -*- c++ -*- */
#ifndef OSGREGISTERFORCEEXEMPLARSCUTOFF_H
#define OSGREGISTERFORCEEXEMPLARSCUTOFF_H


namespace ProtoMol {

  class CubicCellManager;
  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void oSGregisterForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*);
  void oSGregisterForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*);

}
#endif /* REGISTERFORCEEXEMPLARS_H */
