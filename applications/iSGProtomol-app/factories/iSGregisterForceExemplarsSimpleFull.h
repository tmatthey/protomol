/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSSIMPLEFULL_H
#define REGISTERFORCEEXEMPLARSSIMPLEFULL_H


namespace ProtoMol {

  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void iSGregisterForceExemplarsSimpleFull(const PeriodicBoundaryConditions*);
  void iSGregisterForceExemplarsSimpleFull(const VacuumBoundaryConditions*);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
