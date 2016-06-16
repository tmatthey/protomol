#include "PGMReader.h"
#include "stringutilities.h"
#include "Report.h"

using namespace ProtoMol::Report;
using std::string;

namespace ProtoMol
{
	//_________________________________________________________________PGMReader

	PGMReader::PGMReader(): Reader(), myPGM(NULL)
	{
	}

	PGMReader::PGMReader(const std::string& filename): Reader(filename), myPGM(NULL)
	{
	}

	PGMReader::~PGMReader()
	{
		if (myPGM != NULL)
			delete myPGM;
	}

	bool PGMReader::tryFormat()
	{
		if (!open())
			return false;
		string p5, str;

		myFile >> p5;
		getline();
		unsigned int w = 0, h = 0;
		while (myFile >> str)
		{
			if (str != "#")
			{
				w = toInt(str);
				myFile >> h;
				break;
			}
			getline();
		}
		close();
		return (!myFile.fail() && w * h > 0 && (p5 == string("P5") || p5 == string("P3")));
	}

	bool PGMReader::read()
	{
		if (myPGM == NULL)
			myPGM = new PGM();
		return read(*myPGM);
	}

	bool PGMReader::read(PGM& pgm)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;
		string p5, str;
		myFile >> p5;
		getline();
		unsigned int w = 0, h = 0;
		while (myFile >> str)
		{
			if (str != "#")
			{
				w = toInt(str);
				myFile >> h;
				break;
			}
			myComment = str + " " + getline();
		}
		getline();
		getline();
		pgm.resize(w, h);
		if (p5 == string("P5"))
		{
			myFile.read(reinterpret_cast<char*>(pgm.begin()), pgm.size());
		}
		else
		{
			int c;
			for (unsigned char* p = pgm.begin(); p != pgm.end(); p++)
			{
				myFile >> c;
				(*p) = static_cast<unsigned char>(c);
			}
		}
		close();
		return !myFile.fail();
	}

	bool PGMReader::read(PPM& ppm)
	{
		PGM pgm;
		bool res = read(pgm);
		ppm = pgm;
		return res;
	}

	PGM* PGMReader::orphanPGM()
	{
		PGM* tmp = myPGM;
		myPGM = NULL;
		return tmp;
	}

	PGMReader& operator>>(PGMReader& pgmReader, PGM& pgm)
	{
		pgmReader.read(pgm);
		return pgmReader;
	}

	PGMReader& operator>>(PGMReader& pgmReader, PPM& ppm)
	{
		pgmReader.read(ppm);
		return pgmReader;
	}
}
