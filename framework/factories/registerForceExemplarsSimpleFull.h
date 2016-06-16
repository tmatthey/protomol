/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSSIMPLEFULL_H
#define REGISTERFORCEEXEMPLARSSIMPLEFULL_H


namespace ProtoMol
{
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;

	/// registers direct force method prototypes (periodic boundary conditions) to be recognized by the parser
	void registerForceExemplarsSimpleFull(const PeriodicBoundaryConditions*);
	/// registers direct force method prototypes (vacuum) to be recognized by the parser
	void registerForceExemplarsSimpleFull(const VacuumBoundaryConditions*);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
