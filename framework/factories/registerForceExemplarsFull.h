/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARSFULL_H
#define REGISTERFORCEEXEMPLARSFULL_H


namespace ProtoMol
{
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;

	/// registers direct force method prototypes for multiple images (periodic boundary conditions) prototypes to be recognized by the parser
	void registerForceExemplarsFull(const PeriodicBoundaryConditions*);
	/// registers direct force method prototypes for multiple images (vacuum) to be recognized by the parser
	/// should performe as SimpleFull!
	void registerForceExemplarsFull(const VacuumBoundaryConditions*);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
