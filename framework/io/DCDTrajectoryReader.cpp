#include "DCDTrajectoryReader.h"

#include "Report.h"
#include "systemutilities.h"
#include "typeSelection.h"

using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________DCDTrajectoryReader

	DCDTrajectoryReader::DCDTrajectoryReader(): Reader(std::ios::binary), myCoords(NULL), mySwapEndian(false), myFirst(true)
	{
	}

	DCDTrajectoryReader::DCDTrajectoryReader(const std::string& filename): Reader(std::ios::binary, filename), myCoords(NULL), mySwapEndian(false), myFirst(true)
	{
	}

	DCDTrajectoryReader::~DCDTrajectoryReader()
	{
		if (myCoords != NULL)
			delete myCoords;
	}

	bool DCDTrajectoryReader::tryFormat()
	{
		if (!open())
			return false;

		myFile.seekg(0, std::ios::end);
		std::ios::pos_type size = myFile.tellg();
		myFile.seekg(0, std::ios::beg);

		int32 n = 0;
		File::read(reinterpret_cast<char*>(&n), sizeof(int32));
		char coord[5];
		File::read(coord, 4);
		coord[4] = '\0';
		close();

		int32 m = n;
		swapBytes(m);

		if (static_cast<long>(size) >= 104 && (n == 84 || m == 84) && string(coord) == "CORD")
		{
			return !myFile.fail();
		}

		myFile.setstate(std::ios::failbit);
		return false;
	}

	bool DCDTrajectoryReader::read()
	{
		if (myCoords == NULL)
			myCoords = new Vector3DBlock();
		return read(*myCoords);
	}

	bool DCDTrajectoryReader::read(Vector3DBlock& coords)
	{
		if (myFirst)
		{
			// First time ...
			myX.resize(0);
			myY.resize(0);
			myZ.resize(0);

			if (!open())
				return false;
			myFirst = false;

			myFile.seekg(0, std::ios::end);
			std::ios::pos_type size = myFile.tellg();
			myFile.seekg(0, std::ios::beg);

			int32 n = 0;
			File::read(reinterpret_cast<char*>(&n), sizeof(int32));
			char coord[5];
			File::read(coord, 4);
			coord[4] = '\0';

			int32 m = n;
			swapBytes(m);

			// Check endianess and if the header looks ok out ...
			if (static_cast<long>(size) >= 104 && (n == 84 || m == 84) && string(coord) == "CORD")
			{
				if (m == 84)
				{
					mySwapEndian = true;
					report << hint << "[DCDTrajectoryReader::read] Reading " << (ISLITTLEENDIAN ? "big" : "little")
						<< "endian input on " << (ISLITTLEENDIAN ? "little" : "big") << "endian machine." << endr;
				}
			}
			else
			{
				myFile.setstate(std::ios::failbit);
				close();
				return false;
			}

			//       int32 i32;
			//       float4 f32;
			//       myFile.seekg ( 8, std::ios::beg);File::read(reinterpret_cast<char*>(&i32),4); if(mySwapEndian) swapBytes(i32);
			//       report << "NumSets     :"<<i32<<endr;
			//       myFile.seekg (12, std::ios::beg);File::read(reinterpret_cast<char*>(&i32),4); if(mySwapEndian) swapBytes(i32);
			//       report << "Start       :"<<i32<<endr;
			//       myFile.seekg (16, std::ios::beg);File::read(reinterpret_cast<char*>(&i32),4); if(mySwapEndian) swapBytes(i32);
			//       report << "Timesteps   :"<<i32<<endr;
			//       myFile.seekg (20, std::ios::beg);File::read(reinterpret_cast<char*>(&i32),4); if(mySwapEndian) swapBytes(i32);
			//       report << "NAMD writes :"<<i32<<endr;
			//       myFile.seekg (44, std::ios::beg);File::read(reinterpret_cast<char*>(&f32),4); if(mySwapEndian) swapBytes(f32);
			//       report << "length      :"<<f32<<endr;
			//       myFile.seekg (48, std::ios::beg);File::read(reinterpret_cast<char*>(&i32),4); if(mySwapEndian) swapBytes(i32);
			//       report << "Unit cell   :"<<i32<<endr;


			int32 freeIndexes = 0;
			myFile.seekg(40, std::ios::beg);
			File::read(reinterpret_cast<char*>(&freeIndexes), sizeof(int32));
			if (mySwapEndian) swapBytes(freeIndexes);

			// Skip header
			myFile.seekg(84 + 4, std::ios::beg);
			n = 0;
			File::read(reinterpret_cast<char*>(&n), sizeof(int32));
			if (mySwapEndian) swapBytes(n);
			if (n != 84 || myFile.fail())
			{
				myFile.setstate(std::ios::failbit);
				close();
				return false;
			}
			// Skip titles
			int32 l = 0;
			n = 0;
			m = 0;
			File::read(reinterpret_cast<char*>(&n), sizeof(int32));
			if (mySwapEndian) swapBytes(n);
			File::read(reinterpret_cast<char*>(&m), sizeof(int32));
			if (mySwapEndian) swapBytes(m);
			//myFile.seekg (m*80, std::ios::cur);
			myComment.resize(m * 80 + m - 1);
			for (unsigned int i = 0, j = 0; i < (unsigned int)m; i++ , j += 81)
			{
				File::read(&myComment[j], 80);
				myComment[j + 80] = '\n';
			}
			File::read(reinterpret_cast<char*>(&l), sizeof(int32));
			if (mySwapEndian) swapBytes(l);
			if ((n - 4) % 80 != 0 || myFile.fail() || l != n)
			{
				myFile.setstate(std::ios::failbit);
				close();
				return false;
			}


			n = 0;
			l = 0;
			int32 count = 0;
			File::read(reinterpret_cast<char*>(&n), sizeof(int32));
			if (mySwapEndian) swapBytes(n); // 4
			File::read(reinterpret_cast<char*>(&count), sizeof(int32));
			if (mySwapEndian) swapBytes(count); // number of atoms
			File::read(reinterpret_cast<char*>(&l), sizeof(int32));
			if (mySwapEndian) swapBytes(l); // 4
			if (n != 4 || l != 4 || myFile.fail())
			{
				myFile.setstate(std::ios::failbit);
				close();
				return false;
			}

			// Skip free indexes
			if (freeIndexes > 0)
				myFile.seekg(4 * (count - freeIndexes + 2), std::ios::cur);

			myX.resize(count);
			myY.resize(count);
			myZ.resize(count);
		}

		// Read next frame
		// X-dim
		int32 count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count); // number of atoms
		count /= sizeof(int32);
		if ((unsigned int)count != myX.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		File::read(reinterpret_cast<char*>(&(myX[0])), sizeof(float4) * count);
		count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count);
		count /= sizeof(int32);
		if ((unsigned int)count != myX.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		// Y-dim
		count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count); // number of atoms
		count /= sizeof(int32);
		if ((unsigned int)count != myY.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		File::read(reinterpret_cast<char*>(&(myY[0])), sizeof(float4) * count);
		count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count);
		count /= sizeof(int32);
		if ((unsigned int)count != myY.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		// Z-dim
		count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count); // number of atoms
		count /= sizeof(int32);
		if ((unsigned int)count != myZ.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		File::read(reinterpret_cast<char*>(&(myZ[0])), sizeof(float4) * count);
		count = 0;
		File::read(reinterpret_cast<char*>(&count), sizeof(int32));
		if (mySwapEndian) swapBytes(count);
		count /= sizeof(int32);
		if ((unsigned int)count != myZ.size() || myFile.fail())
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}


		// Copy back to right structure
		coords.resize(count);
		for (int i = 0; i < count; ++i)
		{
			if (mySwapEndian) swapBytes(myX[i]);
			coords[i].x = myX[i];
			if (mySwapEndian) swapBytes(myY[i]);
			coords[i].y = myY[i];
			if (mySwapEndian) swapBytes(myZ[i]);
			coords[i].z = myZ[i];
		}
		return !myFile.fail();
	}

	XYZ DCDTrajectoryReader::getXYZ() const
	{
		XYZ res;
		if (myCoords != NULL)
			res.coords = (*myCoords);
		res.names.resize(res.coords.size(), "NONAME");
		return res;
	}

	Vector3DBlock* DCDTrajectoryReader::orphanCoords()
	{
		Vector3DBlock* tmp = myCoords;
		myCoords = NULL;
		return tmp;
	}

	DCDTrajectoryReader& operator>>(DCDTrajectoryReader& dcdTrajectoryReader, XYZ& xyz)
	{
		dcdTrajectoryReader.read(xyz.coords);
		if (xyz.coords.size() != xyz.names.size())
			xyz.names.resize(xyz.coords.size(), "NONAME");
		return dcdTrajectoryReader;
	}

	DCDTrajectoryReader& operator>>(DCDTrajectoryReader& dcdTrajectoryReader, Vector3DBlock& coords)
	{
		dcdTrajectoryReader.read(coords);
		return dcdTrajectoryReader;
	}
}
