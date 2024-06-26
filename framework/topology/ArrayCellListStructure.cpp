#include "ArrayCellListStructure.h"
#include "Report.h"
using namespace ProtoMol::Report;

//#define DEBUG_ARRAYCELLLISTSTRUCTURE

namespace ProtoMol
{
	//_________________________________________________________________ ArrayCellListStructure
	ArrayCellListStructure::ArrayCellListStructure(): valid(false),
	                                                  myArray(ArraySizes(0)(0)(0)),
	                                                  myCellSize(Vector3D(0.0, 0.0, 0.0)),
	                                                  myNX(0),
	                                                  myNY(0),
	                                                  myNZ(0),
	                                                  myBegin(myArray.begin()),
	                                                  myEnd(myArray.end()),
	                                                  myBeginConst(myArray.begin()),
	                                                  myEndConst(myArray.end()),
	                                                  mySize(0)
	{
	}

	void ArrayCellListStructure::initialize(const Vector3D& max, Vector3D cellSize)
	{
		int nx = std::max(1, (int)floor(max.x / cellSize.x + Constant::EPSILON));
		int ny = std::max(1, (int)floor(max.y / cellSize.y + Constant::EPSILON));
		int nz = std::max(1, (int)floor(max.z / cellSize.z + Constant::EPSILON));

		if (nx != myNX || ny != myNY || nz != myNZ || myCellSize != cellSize)
		{
#ifdef DEBUG_ARRAYCELLLISTSTRUCTURE
      report << hint << "CubicCellManager: re-size from ("<<myNX<<","<<myNY<<","<<myNZ<<") to ("<<nx<<","<<ny<<","<<nz<<")."<<endr;
      report << hint << "CubicCellManager: N_x "<<toString(max.x/cellSize.x)<<","<<(int)floor(max.x/cellSize.x+Constant::EPSILON)<<","<<(int)floor(max.x/cellSize.x+1-Constant::EPSILON)<<","<<(int)floor(max.x/cellSize.x+1+Constant::EPSILON)<<endr;
#endif

			if (nx * ny * nz * sizeof(T) >= power<sizeof(size_t) * 8>(2.0))
			{
				report << error
					<< "Your systems is expanding such that the "
					<< "the ratio simulation box / cell size is to big. "
					<< "You may decrease your timestep in your integrator or "
					<< "increase your cell size. End of advice."
					<< endr;
			}
			myCellSize = cellSize;
			myNX = nx;
			myNY = ny;
			myNZ = nz;

			//report << "About to check myArray.resize in ArrayCellListStructure" << endr;
			if (!myArray.resize(ArraySizes(myNX)(myNY)(myNZ)))
				report << error
					<< "[ArrayCellListStructure::initialize] Could not allocate memory for CellListStructure["
					<< myNX << "][" << myNY << "][" << myNZ << "]." << endr;

			//report << "Finished checking myArray.resize in ArrayCellListStructure" << endr;

			int max2 = Constant::MAX_INT_2;
			myMaxNX = max2 - max2 % nx;
			myMaxNY = max2 - max2 % ny;
			myMaxNZ = max2 - max2 % nz;
			myNX1 = -myNX / 2;
			myMaxNX1 = myMaxNX - myNX1;
			myNY1 = -myNY / 2;
			myMaxNY1 = myMaxNY - myNY1;
			myNZ1 = -myNZ / 2;
			myMaxNZ1 = myMaxNZ - myNZ1;

			for (int i = 0; i < myNX; i++)
			{
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
				Array<T, 3>::RefArray<2> z2 = myArray[i];
#else
	RefArray<T,2> z2=myArray[i];
#endif
				for (int j = 0; j < myNY; j++)
				{
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
					Array<T, 3>::RefArray<1> z1 = z2[j];
#else
	  RefArray<T,1> z1=z2[j];
#endif
					for (int k = 0; k < myNZ; k++)
					{
						z1[k].first = T1(i, j, k);
						z1[k].second = -1;
					}
				}
			}
		}
		else
		{
			T* a = myArray.begin();
			for (unsigned int i = 0; i < myArray.size(); i++)
				a[i].second = -1;
			//report << "Cell dimensions seemed to be okay" << endr;
		}

		valid = false;
		myBegin = myArray.begin();
		myEnd = myArray.end();
		myBeginConst = myArray.begin();
		myEndConst = myArray.end();
		mySize = myArray.size();


#ifdef DEBUG_ARRAYCELLLISTSTRUCTURE
    report << plain 
	   <<"[ArrayCellListStructure::initialize] cellsize="<<cellSize<<", box("<<myNX<<","<<myNY<<","<<myNZ<<")."
	   <<endr;
#endif
	}

	void ArrayCellListStructure::updateCache()
	{
		T* a = myArray.begin();
		mySize = 0;
		myBegin = myArray.begin();
		myEnd = myArray.begin();
		myBeginConst = myArray.begin();
		myEndConst = myArray.begin();
		unsigned int count = myArray.size();
		bool first = true;
		for (unsigned int i = 0; i < count; i++)
		{
			if (a[i].second >= 0)
			{
				if (first)
				{
					myBegin = myArray.begin() + i;
					myBeginConst = myArray.begin() + i;
				}
				myEnd = myArray.begin() + i + 1;
				myEndConst = myArray.begin() + i + 1;
				first = false;
				mySize++;
			}
		}
		valid = true;
	}
}
