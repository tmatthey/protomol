/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSOTHER_H
#define REGISTERFORCEEXEMPLARSOTHER_H


namespace ProtoMol
{
	class CubicCellManager;
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;

	/// registers special force prototypes (independent of boundary condtions and cell list manager) to be recognized by the parser
	void registerForceExemplarsOther();
	/// registers special force prototypes (periodic boundary conditionsand cell list manager) to be recognized by the parser
	void registerForceExemplarsOther(const PeriodicBoundaryConditions*, const CubicCellManager*);
	/// registers special force prototypes (periodic boundary conditions) to be recognized by the parser
	void registerForceExemplarsOther(const PeriodicBoundaryConditions*);
	/// registers special force prototypes (vacuum and cell list manager) to be recognized by the parser
	void registerForceExemplarsOther(const VacuumBoundaryConditions*, const CubicCellManager*);
	/// registers special force prototypes (vacuum) to be recognized by the parser
	void registerForceExemplarsOther(const VacuumBoundaryConditions*);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
