/* -*- c++ -*- */
#ifndef ISGREGISTERFORCEEXEMPLARSBONDED_H
#define ISGREGISTERFORCEEXEMPLARSBONDED_H


namespace ProtoMol {

  class PeriodicBoundaryConditions;
  class VacuumBoundaryConditions;

  void iSGregisterForceExemplarsBonded(const PeriodicBoundaryConditions*);
  void iSGregisterForceExemplarsBonded(const VacuumBoundaryConditions*);

}
#endif /* REGISTERFORCEEXEMPLARS_H */
