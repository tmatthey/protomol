/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H
#define REGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H


namespace ProtoMol
{
	class CubicCellManager;
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;

	/// registers fast electrostatic force method prototypes (periodic boundary conditions) to be recognized by the parser
	void registerForceExemplarsFastElectrostatic(const PeriodicBoundaryConditions*, const CubicCellManager*);
	/// registers fast electrostatic force method prototypes (vacuum) to be recognized by the parser
	void registerForceExemplarsFastElectrostatic(const VacuumBoundaryConditions*, const CubicCellManager*);
}
#endif /* REGISTERFORCEEXEMPLARSFASTELECTROSTATIC_H */
