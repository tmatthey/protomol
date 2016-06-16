/*  -*- c++ -*-  */
#ifndef XYZTRAJECTORYREADER_H
#define XYZTRAJECTORYREADER_H

#include "Reader.h"
#include "XYZ.h"

namespace ProtoMol
{
	//_________________________________________________________________XYZTrajectoryReader
	/**
	 * Reads a XYY trajectory files (ASCII). 
	 */
	class XYZTrajectoryReader : public Reader
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		XYZTrajectoryReader();
		explicit XYZTrajectoryReader(const std::string& filename);
		virtual ~XYZTrajectoryReader();

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
		// New methods of class XYZ
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool read(XYZ& xyz);
		bool read(Vector3DBlock& coords, std::vector<std::string>& names);

		XYZ getXYZ() const;
		Vector3DBlock* orphanCoords();
		std::vector<std::string>* orphanNames();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend XYZTrajectoryReader& operator>>(XYZTrajectoryReader& xyzReader, XYZ& xyz);
		friend XYZTrajectoryReader& operator>>(XYZTrajectoryReader& xyzReader, Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Vector3DBlock* myCoords;
		std::vector<std::string>* myNames;
		bool myFirst;
	};

	//____________________________________________________________________________INLINES
}
#endif /* XYZTRAJECTORYREADER_H */
