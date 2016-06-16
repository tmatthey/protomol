/*  -*- c++ -*-  */
#ifndef PGM_H
#define PGM_H

namespace ProtoMol
{
	class PPM;

	//_____________________________________________________________________ PGM
	/**
	  Container for PGM binary image.@n @n
  
	  PGM memory layout:@n @n
  
	  h                @n
	  |  1. ----->     @n
	  |  2. ----->     @n
	  |  ...           @n
	  0  N. ----->     @n
	  y                @n
	     x  0-----w    @n
	*/
	class PGM
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		~PGM();
		PGM();
		PGM(int x, int y);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class PGM
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void set(unsigned int x, unsigned int y, unsigned char a)
		{
			p[(h - y - 1) * w + x] = a;
		}

		void set(unsigned int x, unsigned int y, unsigned char r, unsigned char g, unsigned char b)
		{
			p[(h - y - 1) * w + x] = static_cast<unsigned char>((static_cast<unsigned int>(r) + static_cast<unsigned int>(g) + static_cast<unsigned int>(b)) / 3);
		}

		unsigned char get(unsigned int x, unsigned int y) const
		{
			return p[(h - y - 1) * w + x];
		}

		unsigned char getRed(unsigned int x, unsigned int y) const
		{
			return p[(h - y - 1) * w + x];
		}

		unsigned char getGreen(unsigned int x, unsigned int y) const
		{
			return p[(h - y - 1) * w + x];
		}

		unsigned char getBlue(unsigned int x, unsigned int y) const
		{
			return p[(h - y - 1) * w + x];
		}

		unsigned int width() const
		{
			return w;
		}

		unsigned int height() const
		{
			return h;
		}

		unsigned int size() const
		{
			return w * h;
		}

		void resize(unsigned int width, unsigned int height);

		unsigned char* begin() const
		{
			return &p[0];
		}

		unsigned char* end() const
		{
			return &p[w * h];
		}

		void clear();
		PGM& operator=(const PPM& ppm);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		unsigned int w, h;
		unsigned char* p;
	};
}
#endif /* PGM_H */
