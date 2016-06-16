/*  -*- c++ -*-  */
#ifndef FILE_H
#define FILE_H

#include <fstream>
#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ File
	/** 
	 * Abstract base class for all I/O; readers and writers. The readers and writes
	 * are intend to act STL alike to stream into or from a supported structure or 
	 * container.
	 *
	 *
	 *
	 * NB:
	 * - New writer or reader should never inherit directly from File, but 
	 *   from Reader or Writer.
	 * - Reading binaries one should always use File::read(), rather directly 
	 *   myFile.read, since some compilers like Sun WorkShop CC have problems.
	 * - Be careful when defining constructors or methods with default values;
	 *   especially implicit conversion const char* to string
	 * - File objects can be used as ios_base objects inside expression 
	 *   (e.g., while(dcdReader >> xyz){ ... } )
	 */
	class File
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		File()
		{
		} // Force to use one of the constructors below
	protected:
		explicit File(std::ios::openmode mode);
		File(std::ios::openmode mode, const std::string& filename);

	public:
		virtual ~File();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class File
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void setFilename(const std::string& filename);
		std::string getFilename() const;

		virtual bool open() =0;
		virtual bool open(const std::string& filename) =0;
		virtual bool open(const char* filename) =0;
		// We give a default implementation, but force to
		// choose the default implementation or to reimplement

		bool isAccessible();

		void close();

		operator void*() const;
		bool operator!() const;
		// enable expression testing

	protected:
		std::fstream& read(char* c, std::streamsize count);
		// Redirect of fstream::read  (Sun WorkShop CC does not properly read more than one char ...) 
		std::string getline();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		std::ios::openmode myMode;
	protected:
		std::fstream myFile;
		std::string myFilename;
		std::string myComment; // Reader has get and Writer has set method
	};

	//______________________________________________________________________ INLINES

	inline void File::setFilename(const std::string& filename)
	{
		myFilename = filename;
	}

	inline std::string File::getFilename() const
	{
		return myFilename;
	}
}

#endif /* FILE_H */
