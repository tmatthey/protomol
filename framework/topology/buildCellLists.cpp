#include "buildCellLists.h"
#include "Topology.h"
#include "CubicCellManager.h"
#include "PeriodicBoundaryConditions.h"
#include "VacuumBoundaryConditions.h"
#include "simpleTypes.h"

namespace ProtoMol
{
	//_________________________________________________________________ buildCellLists
	void buildCellLists(const Topology<PeriodicBoundaryConditions, CubicCellManager>* topo,
	                    const Vector3DBlock* positions)
	{
		topo->cellManager.initialize(topo->cellLists, topo->min, topo->max, true);
		Vector3D delta(topo->boundaryConditions.origin() - topo->min);
		CubicCellManager::Cell myCell;
		CubicCellManager::CellListStructure::iterator myCellList;
		CubicCellManager::CellListStructure::iterator end = topo->cellLists.end();
		for (int i = (int)topo->atoms.size() - 1; i >= 0; i--)
		{
			myCell = topo->cellManager.findCell(delta +
				topo->boundaryConditions.minimalPosition((*positions)[i]));
			//Report::report << myCell.x<<","<<myCell.y<<","<<myCell.z<<Report::endr;
			myCellList = topo->cellLists.find(myCell);
			if (myCellList == end)
			{
				// This atom is the first on its cell list, so make a new list for it.
				topo->atoms[i].cellListNext = -1;
				topo->cellLists[myCell] = i;
			}
			else
			{
				topo->atoms[i].cellListNext = myCellList->second;
				myCellList->second = i;
			}
		}
		topo->cellManager.updateCache(topo->cellLists);
	}

	//_________________________________________________________________ buildCellLists
	void buildCellLists(const Topology<VacuumBoundaryConditions, CubicCellManager>* topo,
	                    const Vector3DBlock* positions)
	{
		topo->cellManager.initialize(topo->cellLists, topo->min, topo->max, false);
		Vector3D delta(topo->boundaryConditions.origin() - topo->min);
		CubicCellManager::Cell myCell;
		CubicCellManager::CellListStructure::iterator myCellList;
		CubicCellManager::CellListStructure::iterator end = topo->cellLists.end();
		for (int i = (int)topo->atoms.size() - 1; i >= 0; i--)
		{
			myCell = topo->cellManager.findCell(delta +
				topo->boundaryConditions.minimalPosition((*positions)[i]));
			myCellList = topo->cellLists.find(myCell);
			if (myCellList == end)
			{
				// This atom is the first on its cell list, so make a new list for it.
				topo->atoms[i].cellListNext = -1;
				topo->cellLists[myCell] = i;
			}
			else
			{
				topo->atoms[i].cellListNext = myCellList->second;
				myCellList->second = i;
			}
		}
		topo->cellManager.updateCache(topo->cellLists);
	}
}
