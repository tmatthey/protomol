/*  -*- c++ -*-  */
#ifndef PGMWRITER_H
#define PGMWRITER_H

#include "Writer.h"
#include "PGM.h"
#include "PPM.h"

namespace ProtoMol
{
	//_________________________________________________________________PGMWriter
	/**
	 * Writes a PGM raw (binary) image from PPM or PGM image
	 */
	class PGMWriter : public Writer
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PGMWriter();
		explicit PGMWriter(const std::string& filename);

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
		// New methods of class PGM
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool write(const PGM& pgm);
		bool write(const PPM& ppm);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend PGMWriter& operator<<(PGMWriter& pgmWriter, const PGM& pgm);
		friend PGMWriter& operator<<(PGMWriter& pgmWriter, const PPM& ppm);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
	};

	//____________________________________________________________________________INLINES
}
#endif /* PGMWRITER_H */
