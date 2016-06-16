/*  -*- c++ -*-  */
#ifndef CONFIGURATIONREADER_H
#define CONFIGURATIONREADER_H

#include "Reader.h"
#include "Configuration.h"

namespace ProtoMol
{
	//_________________________________________________________________ConfigurationReader
	/*
	 * Reads and parses a ProtoMol configuaration file. The acutal parsing is
	 * delegated to each entry of the configuration container. The entries have
	 * a keyowrd (identifier) and a value with associated type and constraint.
	 * The parsing is implemented in the traits of the types.
	 */
	class ConfigurationReader : public Reader
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ConfigurationReader();
		explicit ConfigurationReader(const std::string& filename);
		virtual ~ConfigurationReader();

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
		// New methods of class Configuration
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool read(Configuration& config);

		Configuration* orphanConfiguration();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend ConfigurationReader& operator>>(ConfigurationReader& configReader, Configuration& config);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Configuration* myConfig;
	};

	//____________________________________________________________________________INLINES
}
#endif /* CONFIGURATIONXYZREADER_H */
