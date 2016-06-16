#include "EigenvectorTextReader.h"

#include "stringutilities.h"
#include "Report.h"
#include "mathutilities.h"

#include "systemutilities.h"

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::ios;
using std::endl;
using std::find;
using std::stringstream;
using std::fstream;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________EigenvectorTextReader
	EigenvectorTextReader::EigenvectorTextReader():
		Reader(),
		myEigenvectorInfo(NULL)
	{
	}

	EigenvectorTextReader::EigenvectorTextReader(const std::string& filename):
		Reader(filename),
		myEigenvectorInfo(NULL)
	{
	}

	EigenvectorTextReader::EigenvectorTextReader(const char* filename):
		Reader(string(filename)),
		myEigenvectorInfo(NULL)
	{
	}

	EigenvectorTextReader::~EigenvectorTextReader()
	{
		if (myEigenvectorInfo != NULL)
			delete myEigenvectorInfo;
	}

	bool EigenvectorTextReader::read()
	{
		if (myEigenvectorInfo == NULL)
			myEigenvectorInfo = new EigenvectorInfo();
		return read(*myEigenvectorInfo);
	}

	bool EigenvectorTextReader::read(EigenvectorInfo& ei)
	{
		if (!open())
			return false;

		int num, num1;
		double ev;
		myFile >> num;
		report << plain << num << endr;
		myFile >> num1;
		report << plain << num1 << endr;

		myFile >> ev;
		report << plain << ev << endr;

		string str;
		for (unsigned int i = 0; i < 4; i++)
		{
			myFile >> str;
			report << plain << str << endr;
		}

		ei.myEigenvectorLength = num;
		ei.myNumEigenvectors = num1;
		ei.myMaxEigenvalue = ev;
		ei.initializeEigenvectors();


		for (unsigned int i = 0; i < num; i++)
		{
			int x;
			double y;
			for (unsigned int j = 0; j < num1; j++)
			{
				myFile >> x;
				//report << plain << x <<endr;
				for (unsigned int k = 0; k < 3; k++)
				{
					myFile >> y;
					//report << plain << y <<endr;
					ei.myEigenvectors[i * 3 * num1 + j * 3 + k] = y;
				}
			}
		}

		return true;
	}


	EigenvectorTextReader& operator>>(EigenvectorTextReader& eigenvectorReader, EigenvectorInfo& info)
	{
		eigenvectorReader.read(info);
		return eigenvectorReader;
	}
}
