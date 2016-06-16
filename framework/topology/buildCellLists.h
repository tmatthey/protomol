/* -*- c++ -*- */
#ifndef BUILDCELLLIST_H
#define BUILDCELLLIST_H

namespace ProtoMol
{
	class PeriodicBoundaryConditions;
	class VacuumBoundaryConditions;
	class CubicCellManager;
	class Vector3DBlock;
	template <class TBoundaryConditions, class TCellManager>
	class Topology;

	/// builds the cell list for periodic boundary conditions
	void buildCellLists(const Topology<PeriodicBoundaryConditions, CubicCellManager>* topo,
	                    const Vector3DBlock* positions);

	/// build the cell list for vacuum
	void buildCellLists(const Topology<VacuumBoundaryConditions, CubicCellManager>* topo,
	                    const Vector3DBlock* positions);
}
#endif /* BUILDCELLLIST_H */
