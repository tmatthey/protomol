/*  -*- c++ -*-  */
#ifndef READER_H
#define READER_H

#include "File.h"

namespace ProtoMol
{
	//_________________________________________________________________ Reader
	/**
	 * Base class of readers
	 */
	class Reader : public File
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		Reader();
		explicit Reader(const std::string& filename);
		/// To open with special file flags, std::ios::in is set
		explicit Reader(std::ios::openmode mode);
		/// To open with special file flags, std::ios::in is set
		Reader(std::ios::openmode mode, const std::string& filename);
	public:
		virtual ~Reader();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Reader
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool tryFormat() =0;
		/// Simple test, true if it the format might be correct/readable
		virtual bool read() =0;
		const std::string& getComment() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
	};

	//______________________________________________________________________ INLINES
}

#endif /* READER_H */
