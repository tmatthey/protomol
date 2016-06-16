#include "PGM.h"
#include "PPM.h"
#include <vector> // NULL definition

namespace ProtoMol
{
	//_____________________________________________________________________ PGM
	PGM::~PGM()
	{
		delete [] p;
	}

	PGM::PGM(): w(0), h(0), p(NULL)
	{
	}

	PGM::PGM(int x, int y): w(x), h(y), p(new unsigned char[x * y])
	{
		clear();
	}

	void PGM::clear()
	{
		for (unsigned int i = 0; i < w * h; i++)
			p[i] = static_cast<unsigned char>(0);
	}

	void PGM::resize(unsigned int width, unsigned int height)
	{
		if (p != NULL && width * height != w * h)
		{
			delete [] p;
			p = NULL;
		}
		if (p == NULL && width * height > 0)
			p = new unsigned char[width * height];
		w = width;
		h = height;
		clear();
	}

	PGM& PGM::operator=(const PPM& ppm)
	{
		resize(ppm.width(), ppm.height());
		unsigned char* p0 = ppm.begin();
		for (unsigned char* p1 = begin(); p1 != end(); p1++)
		{
			(*p1) = static_cast<unsigned char>((static_cast<unsigned int>(p0[0]) + static_cast<unsigned int>(p0[1]) + static_cast<unsigned int>(p0[2])) / 3);
			p0 += 3;
		}
		return (*this);
	}
}
