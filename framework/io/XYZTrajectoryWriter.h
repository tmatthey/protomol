/*  -*- c++ -*-  */
#ifndef XYZTRAJECTORYWRITER_H
#define XYZTRAJECTORYWRITER_H

#include "Writer.h"
#include "XYZ.h"
#include "Atom.h"
#include "AtomType.h"

namespace ProtoMol
{
	//_________________________________________________________________XYZTrajectoryWriter
	/**
	 * Writes XYZ trajectories (ASCII) and updates the number of coordinate sets 
	 * after each write, no need to know the final number of sets.
	 */
	class XYZTrajectoryWriter : public Writer
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		XYZTrajectoryWriter();
		explicit XYZTrajectoryWriter(const std::string& filename);
		virtual ~XYZTrajectoryWriter();

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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Writer
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual bool write();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class XYZ
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool open(const XYZ& xyz);
		bool open(const std::vector<std::string>& names);
		bool open(const std::vector<Atom>& atoms, const std::vector<AtomType>& atomTypes);

		bool open(const std::string& filename, const XYZ& xyz);
		bool open(const std::string& filename, const std::vector<std::string>& names);
		bool open(const std::string& filename, const std::vector<Atom>& atoms, const std::vector<AtomType>& atomTypes);

		bool write(const XYZ& xyz);
		bool write(const Vector3DBlock& coords);
		bool write(const Vector3DBlock& coords, const std::vector<std::string>& names);
		bool write(const Vector3DBlock& coords, const std::vector<Atom>& atoms, const std::vector<AtomType>& atomTypes);

		void setNames(const XYZ& xyz);
		void setNames(const std::vector<std::string>& names);
		void setNames(const std::vector<Atom>& atoms, const std::vector<AtomType>& atomTypes);

	private:
		bool reopen();
		void setCoords(const Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend XYZTrajectoryWriter& operator<<(XYZTrajectoryWriter& xyzWriter, const XYZ& xyz);
		friend XYZTrajectoryWriter& operator<<(XYZTrajectoryWriter& xyzWriter, const Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		const Vector3DBlock* myCoords;
		const std::vector<std::string>* myNames;
		const std::vector<Atom>* myAtoms;
		const std::vector<AtomType>* myAtomTypes;
		bool myFirst;
	};

	//____________________________________________________________________________INLINES
}
#endif /* XYZTRAJECTORYWRITER_H */
