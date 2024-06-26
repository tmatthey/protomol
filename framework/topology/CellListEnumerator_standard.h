/* -*- c++ -*- */
#ifndef CELLLISTENUMERATOR_STANDARD_H
#define CELLLISTENUMERATOR_STANDARD_H

#include "CellListEnumerator.h"

#include "CubicCellManager.h"
#include "VacuumBoundaryConditions.h"
#include "Topology.h"
#include <vector>

//#define DEBUG_CELLLISTENUMERATOR_STANDARD
namespace ProtoMol
{
	//_________________________________________________________________ CellListEnumerator_standard

	/**
	 * Specialization of the cell enumerator for vacuum and cubic cell manager
	 */
	template <>
	class CellListEnumerator<VacuumBoundaryConditions, CubicCellManager>
	{
	public:
		//typedef pair<int,int> CellPair; // The first atom from each cell list in the pair
		struct CellPair
		{
			int first;
			int second;
		};

	public:
		CellListEnumerator(): myCutoff(-1.0),
		                      myCellSize(Vector3D(-1.0, -1.0, -1.0)),
		                      myMax(CubicCellManager::Cell(-1, -1, -1))
		{
		}

		void initialize(const Topology<VacuumBoundaryConditions, CubicCellManager>* topo,
		                Real cutoff)
		{
			the_beginning = topo->cellLists.begin();
			the_end = topo->cellLists.end();

			i = the_beginning;
			j = i;

			the_first = 0;
			the_counter = the_first;
			myCellListStruct = &(topo->cellLists);

			if (myCutoff != cutoff ||
				topo->cellManager.getRealCellSize() != myCellSize ||
				!(myMax == topo->cellManager.findCell(topo->max - topo->min)))
			{
				myCellSize = topo->cellManager.getRealCellSize();
				myCutoff = cutoff;

				myMax = topo->cellManager.findCell(topo->max - topo->min);

				// 	Report::report <<"New: "<<myMax.x<<","<< myMax.y<<","<< myMax.z<<Report::endr;
				Real cutoff2 = cutoff * cutoff;

				CubicCellManager::Cell zero(0, 0, 0);
				myDeltaList.clear();
				int nx = (int)(cutoff / myCellSize.x + 1.0 + Constant::EPSILON);
				int ny = (int)(cutoff / myCellSize.y + 1.0 + Constant::EPSILON);
				int nz = (int)(cutoff / myCellSize.z + 1.0 + Constant::EPSILON);
				// Do not consider deltas bigger than the dimesion of the simulation box
				int n0 = std::min(nx, topo->cellLists.getDimX() - 1);
				int n1 = std::min(ny, topo->cellLists.getDimY() - 1);
				int n2 = std::min(nz, topo->cellLists.getDimZ() - 1);
				Real xx = myCellSize.x * myCellSize.x;
				Real yy = myCellSize.y * myCellSize.y;
				Real zz = myCellSize.z * myCellSize.z;
				for (int k = -n0; k <= n0; k++)
				{
					int x = abs(k) - 1;
					if (x < 0) x = 0;
					Real d0 = x * x * xx;
					for (int l = -n1; l <= n1; l++)
					{
						int y = abs(l) - 1;
						if (y < 0) y = 0;
						Real d1 = d0 + y * y * yy;
						for (int m = -n2; m <= n2; m++)
						{
							int z = abs(m) - 1;
							if (z < 0) z = 0;
							if (d1 + z * z * zz < cutoff2)
							{
								CubicCellManager::Cell delta(k, l, m);
								if (zero < delta)
									myDeltaList.push_back(delta);
							}
						}
					}
				}
				std::sort(myDeltaList.begin(), myDeltaList.end());
				// 	for(unsigned int i=0;i<myDeltaList.size();i++){
				// 	  Report::report <<"Delta["<<i<<"] "<< myDeltaList[i].x<<","<< myDeltaList[i].y<<","<< myDeltaList[i].z<<Report::endr;
				// 	}

				//Report::report << Report::debug(2) << "Cell list algorithm: "<< topo->cellLists.getDimX()*topo->cellLists.getDimY()*topo->cellLists.getDimZ()<<", "<<myDeltaList.size()<<", "<< topo->cellLists.getDimX()*topo->cellLists.getDimY()*topo->cellLists.getDimZ()*myCellSize.x*myCellSize.y*myCellSize.z << ", "<< myDeltaList.size()*myCellSize.x*myCellSize.y*myCellSize.z <<Report::endr;
			}
			the_last = myDeltaList.size();
		}

		/// retrieve the current cell pair
		void get(CellPair& cp)
		{
			cp.first = i->second;
			cp.second = j->second;
		}

		/// retrieve the current cell pair
		void get(int& a, int& b)
		{
			a = i->second;
			b = j->second;
		}

		/// if the cells of the current cell pair are the same
		bool notSameCell()
		{
			return (i != j);
		}

		/// reached the end of the list of cell pairs
		bool done()
		{
			return (i == the_end);
		}

		/// goto the end of pair of cells with the same first cell
		void gotoEndPair()
		{
			j = the_end;
		};

		/// advance by inc in the cell list in respect to first cell 
		void nextNewPair(int inc)
		{
			if (inc < 1)
				return;
			while (inc > 0 && i != the_end)
			{
				++i;
				--inc;
				while (i != the_end && i->second < 0)
					++i;
			}
			j = i;
		}

		/// get next pair
		void next()
		{
			if (i != the_end)
			{
				j = the_end;
				while (the_counter != the_last &&
					the_end == (j = myCellListStruct->find(i->first + myDeltaList[the_counter++])));

				if (j == the_end)
				{
					++i;
					while (i != the_end && i->second < 0)
						++i;
					j = i;
					the_counter = the_first;
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		CubicCellManager::CellListStructure::const_iterator i, j, the_beginning, the_end;
		Real myCutoff;
		std::vector<CubicCellManager::Cell> myDeltaList;
		int the_first, the_last, the_counter;
		const CubicCellManager::CellListStructure* myCellListStruct;
		Vector3D myCellSize;
		CubicCellManager::Cell myMax;
	};
}
#endif /* CELLLISTENUMERATOR_STANDARD_H */
