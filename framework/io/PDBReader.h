/*  -*- c++ -*-  */
#ifndef PDBREADER_H
#define PDBREADER_H

#include "Reader.h"
#include "PDB.h"
#include "XYZ.h"

namespace ProtoMol
{
	//_________________________________________________________________PDBReader
	/**
	 * Reads a PDB (ASCII) file, ATOM, HETATOM and TER record only.
	 */
	class PDBReader : public Reader
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PDBReader();
		explicit PDBReader(const std::string& filename);
		virtual ~PDBReader();

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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Reader
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool tryFormat();
		virtual bool read();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class PDB
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool read(PDB& pdb);
		bool read(Vector3DBlock& coords, std::vector<PDB::Atom>& atoms, std::vector<PDB::Ter>& ters);

		PDB getPDB() const;
		Vector3DBlock* orphanCoords();
		std::vector<PDB::Atom>* orphanAtoms();
		std::vector<PDB::Ter>* orphanTers();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend PDBReader& operator>>(PDBReader& pdbReader, PDB& pdb);
		friend PDBReader& operator>>(PDBReader& pdbReader, XYZ& xyz);
		friend PDBReader& operator>>(PDBReader& pdbReader, Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Vector3DBlock* myCoords;
		std::vector<PDB::Atom>* myAtoms;
		std::vector<PDB::Ter>* myTers;
	};

	//____________________________________________________________________________INLINES
}
#endif /* PDBREADER_H */
