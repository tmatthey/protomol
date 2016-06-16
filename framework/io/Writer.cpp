#include "Writer.h"

using std::string;

namespace ProtoMol
{
	//_________________________________________________________________ Writer
	Writer::Writer(): File(std::ios::out | std::ios::trunc)
	{
	}

	Writer::Writer(const string& filename): File(std::ios::out | std::ios::trunc, filename)
	{
		File::open();
	}

	Writer::Writer(std::ios::openmode mode): File(std::ios::out | mode)
	{
	}

	Writer::Writer(std::ios::openmode mode, const string& filename): File(std::ios::out | mode, filename)
	{
		File::open();
	}

	Writer::~Writer()
	{
	}

	void Writer::setComment(const std::string& comment)
	{
		myComment = comment;
	}
}
