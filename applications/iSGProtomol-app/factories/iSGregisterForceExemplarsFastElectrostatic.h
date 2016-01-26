/* -*- c++ -*- */
#ifndef ISGREGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H
#define ISGREGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H


namespace ProtoMol {

  class CubicCellManager;
  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void iSGregisterForceExemplarsFastElectrostatic(const PeriodicBoundaryConditions*, const CubicCellManager*);
  void iSGregisterForceExemplarsFastElectrostatic(const VacuumBoundaryConditions*, const CubicCellManager*);


}
#endif /* REGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H */
