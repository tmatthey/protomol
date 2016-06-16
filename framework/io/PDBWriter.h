/*  -*- c++ -*-  */
#ifndef PDBWRITER_H
#define PDBWRITER_H

#include "Writer.h"
#include "PDB.h"

namespace ProtoMol
{
	//_________________________________________________________________PDBWriter
	/**
	 * Writes a PDB (ASCII) file, ATOM record only.
	 */
	class PDBWriter : public Writer
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PDBWriter();
		explicit PDBWriter(const std::string& filename);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class File
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool open()
		{
			return File::open();
		}

		virtual bool open(const std::string& filename)
		{
			return File::open(filename);
		}

		virtual bool open(const char* filename)
		{
			return File::open(filename);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class PDB
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool write(const PDB& pdb);
		bool write(const Vector3DBlock& coords, const PDB& pdb);
		bool write(const Vector3DBlock& coords, const std::vector<PDB::Atom>& atoms);
		bool write(const Vector3DBlock& coords, const std::vector<PDB::Atom>& atoms, const std::vector<PDB::Ter>& ters);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend PDBWriter& operator<<(PDBWriter& pdbWriter, const PDB& pdb);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
	};

	//____________________________________________________________________________INLINES
}
#endif /* PDBWRITER_H */
