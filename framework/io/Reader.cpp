#include "Reader.h"

using std::string;

namespace ProtoMol
{
	//_________________________________________________________________ Reader
	Reader::Reader(): File(std::ios::in)
	{
	}

	Reader::Reader(const string& filename): File(std::ios::in, filename)
	{
		File::open();
	}

	Reader::Reader(std::ios::openmode mode): File(std::ios::in | mode)
	{
	}

	Reader::Reader(std::ios::openmode mode, const string& filename): File(std::ios::in | mode, filename)
	{
		File::open();
	}

	Reader::~Reader()
	{
	}

	const string& Reader::getComment() const
	{
		return myComment;
	}
}
