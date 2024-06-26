#include "XYZReader.h"

#include "Report.h"
#include "stringutilities.h"

using std::string;
using std::vector;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________XYZReader

	XYZReader::XYZReader(): Reader(), myCoords(NULL), myNames(NULL)
	{
	}

	XYZReader::XYZReader(const std::string& filename): Reader(filename), myCoords(NULL), myNames(NULL)
	{
	}

	XYZReader::~XYZReader()
	{
		if (myCoords != NULL)
			delete myCoords;
		if (myNames != NULL)
			delete myNames;
	}

	bool XYZReader::tryFormat()
	{
		if (!open())
			return false;

		// Number of atoms
		unsigned int count = 0;
		myFile >> count;
		if (myFile.good() && count == 0)
		{
			close();
			return true;
		}

		string str(getline());

		// Comment
		str = getline();
		bool testTrajectory = isUInt(str) && (toUInt(str) > 1);

		// Atoms
		Real x;
		if (testTrajectory)
		{
			unsigned int frames = toUInt(str);
			for (unsigned int i = 0; i < count && !myFile.fail(); ++i)
			{
				myFile >> str >> x >> x >> x;
			}
			if (frames > 1)
				myFile >> count >> str >> x >> x >> x;
		}
		else
		{
			myFile >> str >> x >> x >> x;
		}
		myFile.close();
		return !myFile.fail();
	}

	bool XYZReader::read()
	{
		if (myCoords == NULL)
			myCoords = new Vector3DBlock();
		if (myNames == NULL)
			myNames = new vector<string>();
		return read(*myCoords, *myNames);
	}

	bool XYZReader::read(XYZ& xyz)
	{
		return read(xyz.coords, xyz.names);
	}

	bool XYZReader::read(Vector3DBlock& coords, std::vector<std::string>& names)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;

		// Number of atoms
		unsigned int n = 0;
		myFile >> n;
		string str(getline());

		// Comment
		str = getline();
		if (myFile.fail())
		{
			close();
			return false;
		}

		if (isInt(str))
			report << hint << "[XYZReader::read] Does also match XYZ trajectory format." << endr;

		myComment = str;

		coords.resize(n);
		names.resize(n);

		// Read atoms
		for (unsigned int i = 0; i < n && !myFile.fail(); ++i)
		{
			myFile >> names[i] >> coords[i].x >> coords[i].y >> coords[i].z;
		}
		close();
		return !myFile.fail();
	}

	XYZ XYZReader::getXYZ() const
	{
		XYZ res;
		if (myCoords != NULL)
			res.coords = (*myCoords);
		if (myNames != NULL)
			res.names = (*myNames);
		return res;
	}

	Vector3DBlock* XYZReader::orphanCoords()
	{
		Vector3DBlock* tmp = myCoords;
		myCoords = NULL;
		return tmp;
	}

	vector<string>* XYZReader::orphanNames()
	{
		vector<string>* tmp = myNames;
		myNames = NULL;
		return tmp;
	}

	XYZReader& operator>>(XYZReader& xyzReader, XYZ& xyz)
	{
		xyzReader.read(xyz.coords, xyz.names);
		return xyzReader;
	}

	XYZReader& operator>>(XYZReader& xyzReader, Vector3DBlock& coords)
	{
		if (xyzReader.myNames == NULL)
			xyzReader.myNames = new vector<string>();
		xyzReader.read(coords, *xyzReader.myNames);
		return xyzReader;
	}
}
