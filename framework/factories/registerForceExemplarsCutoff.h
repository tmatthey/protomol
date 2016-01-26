/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSCUTOFF_H
#define REGISTERFORCEEXEMPLARSCUTOFF_H


namespace ProtoMol {

  class CubicCellManager;
  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  /// registers cutoff force method prototypes (periodic boundary conditions) to be recognized by the parser
  void registerForceExemplarsCutoff(const PeriodicBoundaryConditions*, const CubicCellManager*);
  /// registers cutoff force method  prototypes(vacuum) to be recognized by the parser
  void registerForceExemplarsCutoff(const VacuumBoundaryConditions*, const CubicCellManager*);

}
#endif /* REGISTERFORCEEXEMPLARS_H */
