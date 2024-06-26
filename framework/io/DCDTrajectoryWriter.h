/*  -*- c++ -*-  */
#ifndef DCDTRAJECTORYWRITER_H
#define DCDTRAJECTORYWRITER_H

#include "Writer.h"
#include "XYZ.h"
#include "systemutilities.h"
#include "typeSelection.h"

namespace ProtoMol
{
	//_________________________________________________________________DCDTrajectoryWriter
	/**
	 * Writes DCD trajectories and updates the number of coordinate sets 
	 * after each write, no need to know the final number of sets.
	 */
	class DCDTrajectoryWriter : public Writer
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
		explicit DCDTrajectoryWriter(Real timestep = 1.0, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		explicit DCDTrajectoryWriter(const std::string& filename, Real timestep = 1.0, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		explicit DCDTrajectoryWriter(const char* filename, Real timestep = 1.0, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		explicit DCDTrajectoryWriter(std::ios::openmode mode, const std::string& filename, Real timestep = 1.0, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		explicit DCDTrajectoryWriter(std::ios::openmode mode, const char* filename, Real timestep = 1.0, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		// Need this implementation, otherwise const char* will bee converted to bool or int ...

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class File
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool open()
		{
			myFirst = false;
			return File::open();
		}

		virtual bool open(const std::string& filename)
		{
			myFirst = false;
			return File::open(filename);
		}

		virtual bool open(const char* filename)
		{
			myFirst = false;
			return File::open(filename);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DCD
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool open(Real timestep, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		bool open(const std::string& filename, Real timestep, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);
		bool open(const char* filename, Real timestep, unsigned int firststep = 1, bool isLittleEndian = ISLITTLEENDIAN);

		bool write(const Vector3DBlock& coords);

		void setLittleEndian(bool littleEndian);
		void setTimestep(Real timestep);
		void setFirststep(unsigned int firststep);

		bool reopen(unsigned int numAtoms);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend DCDTrajectoryWriter& operator<<(DCDTrajectoryWriter& dcdWriter, const Vector3DBlock& coords);
		friend DCDTrajectoryWriter& operator<<(DCDTrajectoryWriter& dcdWriter, const XYZ& xyz);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		bool myIsLittleEndian;
		unsigned int myFirstStep;
		Real myTimeStep;
		std::vector<float4> myX;
		std::vector<float4> myY;
		std::vector<float4> myZ;
		bool myFirst; // flag if the file has to be opend and cleared
	};

	//____________________________________________________________________________INLINES
}
#endif /* DCDTRAJECTORYWRITER_H */
