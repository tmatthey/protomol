#include "PPMReader.h"
#include "Report.h"
#include "stringutilities.h"

using namespace ProtoMol::Report;
using std::string;

namespace ProtoMol
{
	//_________________________________________________________________PPMReader

	PPMReader::PPMReader(): Reader(), myPPM(NULL)
	{
	}

	PPMReader::PPMReader(const std::string& filename): Reader(filename), myPPM(NULL)
	{
	}

	PPMReader::~PPMReader()
	{
		if (myPPM != NULL)
			delete myPPM;
	}

	bool PPMReader::tryFormat()
	{
		if (!open())
			return false;
		string p6, str;

		myFile >> p6;
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
		return (!myFile.fail() && w * h > 0 && (p6 == string("P6") || p6 == string("P3") || p6 == string("P5") || p6 == string("P2")));
	}

	bool PPMReader::read()
	{
		if (myPPM == NULL)
			myPPM = new PPM();
		return read(*myPPM);
	}

	bool PPMReader::read(PPM& ppm)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;
		string p6, str;
		myFile >> p6;
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
		ppm.resize(w, h);
		if (p6 == string("P6"))
		{
			myFile.read(reinterpret_cast<char*>(ppm.begin()), ppm.size());
		}
		else if (p6 == string("P3"))
		{
			int c;
			for (unsigned char* p = ppm.begin(); p != ppm.end(); p++)
			{
				myFile >> c;
				(*p) = static_cast<unsigned char>(c);
			}
		}
		else if (p6 == string("P5"))
		{
			myFile.read(reinterpret_cast<char*>(ppm.begin()), ppm.size() / 3);
			unsigned char* p = ppm.begin();
			for (int i = ppm.size() / 3 - 1; i >= 0; i--)
			{
				p[i * 3 + 2] = p[i];
				p[i * 3 + 1] = p[i];
				p[i * 3 + 0] = p[i];
			}
		}
		else if (p6 == string("P2"))
		{
			int c;
			unsigned char* p = ppm.begin();
			for (unsigned int i = 0; i < ppm.size() / 3; i++)
			{
				myFile >> c;
				p[0] = static_cast<unsigned char>(c);
				p[1] = static_cast<unsigned char>(c);
				p[2] = static_cast<unsigned char>(c);
				p += 3;
			}
		}

		close();
		return !myFile.fail();
	}

	bool PPMReader::read(PGM& pgm)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;
		string p6, str;
		myFile >> p6;
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
		if (p6 == string("P6"))
		{
			PPM ppm;
			myFile.read(reinterpret_cast<char*>(ppm.begin()), ppm.size());
			pgm = ppm;
		}
		else if (p6 == string("P3"))
		{
			int r, g, b;
			for (unsigned char* p = pgm.begin(); p != pgm.end(); p++)
			{
				myFile >> r >> g >> b;
				(*p) = static_cast<unsigned char>((r + g + b) / 3);
			}
		}
		else if (p6 == string("P5"))
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


	PPM* PPMReader::orphanPPM()
	{
		PPM* tmp = myPPM;
		myPPM = NULL;
		return tmp;
	}

	PPMReader& operator>>(PPMReader& ppmReader, PPM& ppm)
	{
		ppmReader.read(ppm);
		return ppmReader;
	}

	PPMReader& operator>>(PPMReader& ppmReader, PGM& pgm)
	{
		ppmReader.read(pgm);
		return ppmReader;
	}
}
