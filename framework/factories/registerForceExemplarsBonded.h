/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSBONDED_H
#define REGISTERFORCEEXEMPLARSBONDED_H


namespace ProtoMol
{
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;

	/// registers bond force method prototypes (periodic boundary conditions) to be recognized by the parser
	void registerForceExemplarsBonded(const PeriodicBoundaryConditions*);
	/// registers bond force method prototypes (vacuum) to be recognized by the parser
	void registerForceExemplarsBonded(const VacuumBoundaryConditions*);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
