#include "TempFileReader.h"
#include <iostream>

using std::cout;
using std::endl;

namespace ProtoMol
{
	bool TempFileReader::tryFormat()
	{
		if (!open())
			return false;
		Real num = 0;
		myFile >> num;
		if (myFile.good() && num != 0)
		{
			close();
			return true;
		}
		return false;
	}

	bool TempFileReader::read()
	{
		return read(myTemp);
	}

	bool TempFileReader::read(Real& temp)
	{
		if (!tryFormat())
			return false;
		open();
		for (int i = 0; i <= node; i++)
		{
			if (myFile.eof())
			{
				close();
				return false;
			}
			myFile >> temp;
		}
		close();
		return true;
	}

	void TempFileReader::setNodeNum(int num)
	{
		node = num;
	}

	TempFileReader& operator>>(TempFileReader& reader, Real& temp)
	{
		reader.read(temp);
		return reader;
	}
}
