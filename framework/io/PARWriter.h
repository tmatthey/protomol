/*  -*- c++ -*-  */
#ifndef PARWRITER_H
#define PARWRITER_H

#include "Writer.h"
#include "PAR.h"

namespace ProtoMol
{
	//_________________________________________________________________PARWriter
	/**
	 * Writes PAR-Charmm 19/27 files (ASCII)
	 */
	class PARWriter : public Writer
	{
	public:
		static const PAR::CharmmTypeEnum defaultCharmmType;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		explicit PARWriter(PAR::CharmmTypeEnum charmmType = defaultCharmmType);
		explicit PARWriter(const std::string& filename, PAR::CharmmTypeEnum charmmType = defaultCharmmType);
		explicit PARWriter(const char* filename, PAR::CharmmTypeEnum charmmType = defaultCharmmType);
		// Need this implementation, otherwise const char* will bee converted to bool or int ...

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
		// New methods of class PAR
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool open(const std::string& filename, PAR::CharmmTypeEnum charmmType);
		bool open(const char* filename, PAR::CharmmTypeEnum charmmType);
		bool open(PAR::CharmmTypeEnum charmmType);

		bool write(const PAR& par);

		void setCharmmType(PAR::CharmmTypeEnum charmmType);
		PAR::CharmmTypeEnum getCharmmType() const;

	private:
		void setPAR(const PAR& par);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend PARWriter& operator<<(PARWriter& parWriter, const PAR& par);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		PAR::CharmmTypeEnum myCharmmType;
	};

	//____________________________________________________________________________INLINES
}
#endif /* PARWRITER_H */
