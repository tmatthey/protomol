/*  -*- c++ -*-  */
#ifndef DCDTRAJECTORYREADER_H
#define DCDTRAJECTORYREADER_H

#include "Reader.h"
#include "XYZ.h"
#include "typeSelection.h"

namespace ProtoMol
{
	//_________________________________________________________________DCDTrajectoryReader
	/**
	 * Reads a DCD trajectory file, frame by frame. Automatic endianess
	 * detection.
	 */
	class DCDTrajectoryReader : public Reader
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Typedef
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		typedef TypeSelection::Int<4>::type int32;
		typedef TypeSelection::Float<4>::type float4;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		DCDTrajectoryReader();
		explicit DCDTrajectoryReader(const std::string& filename);
		virtual ~DCDTrajectoryReader();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class File
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool open()
		{
			myFirst = true;
			return File::open();
		}

		virtual bool open(const std::string& filename)
		{
			myFirst = true;
			return File::open(filename);
		}

		virtual bool open(const char* filename)
		{
			myFirst = true;
			return File::open(filename);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Reader
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool tryFormat();
		virtual bool read();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DCDTrajectory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:

		bool read(Vector3DBlock& coords);
		XYZ getXYZ() const;

		Vector3DBlock* orphanCoords();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend DCDTrajectoryReader& operator>>(DCDTrajectoryReader& DCDTrajectoryReader, XYZ& xyz);
		friend DCDTrajectoryReader& operator>>(DCDTrajectoryReader& DCDTrajectoryReader, Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		std::vector<float4> myX;
		std::vector<float4> myY;
		std::vector<float4> myZ;
		Vector3DBlock* myCoords;
		bool mySwapEndian;
		bool myFirst;
	};

	//____________________________________________________________________________INLINES
}
#endif /* DCDTRAJECTORYREADER_H */
