/*  -*- c++ -*-  */
#ifndef PPMWRITER_H
#define PPMWRITER_H

#include "Writer.h"
#include "PPM.h"
#include "PGM.h"

namespace ProtoMol
{
	//_________________________________________________________________PPMWriter
	/**
	 * Writes a PPM raw (binary) image from a PPM (color) or PGM (gray) image
	 */
	class PPMWriter : public Writer
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PPMWriter();
		explicit PPMWriter(const std::string& filename);

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
		// New methods of class PPM
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		/// write PPM from PPM (color) image
		bool write(const PPM& ppm);
		/// write PPM from PGM (gray) image
		bool write(const PGM& pgm);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend PPMWriter& operator<<(PPMWriter& ppmWriter, const PPM& ppm);
		friend PPMWriter& operator<<(PPMWriter& ppmWriter, const PGM& pgm);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
	};

	//____________________________________________________________________________INLINES
}
#endif /* PPMWRITER_H */
