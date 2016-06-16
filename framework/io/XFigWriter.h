/*  -*- c++ -*-  */
#ifndef XFIGWRITER_H
#define XFIGWRITER_H

#include "Writer.h"
#include "XYZ.h"
#include "Matrix3by3.h"

namespace ProtoMol
{
	//_________________________________________________________________XFigWriter
	/**
	 * Writes a XFig file
	 */
	class XFigWriter : public Writer
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		XFigWriter();
		explicit XFigWriter(const std::string& filename);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class File
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool open()
		{
			return File::open();
		}

		virtual bool open(const std::string& filename)
		{
			return File::open(filename);
		}

		virtual bool open(const char* filename)
		{
			return File::open(filename);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class XFigWriter
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		bool write(const XYZ& xyz);
		bool write(const Vector3DBlock& coords, const std::vector<std::string>& names);

		Real setRadius(Real radius);

		Real getRadius() const
		{
			return myRadius;
		}

		Real setFromZ(Real FromZ);

		Real getFromZ() const
		{
			return myFromZ;
		}

		Real setToZ(Real ToZ);

		Real getToZ() const
		{
			return myToZ;
		}

		bool setColor(bool color);

		bool getColor() const
		{
			return myColor;
		}

		bool setAxes(bool axes);

		bool getAxes() const
		{
			return myAxes;
		}

		bool setLegend(bool legend);

		bool getLegend() const
		{
			return myLegend;
		}

		bool setAutofit(bool autofit);

		bool getAutofit() const
		{
			return myAutofit;
		}

		Matrix3by3 setTransformation(const Matrix3by3& mat);

		Matrix3by3 getTransformation() const
		{
			return myMat;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend XFigWriter& operator<<(XFigWriter& xfigWriter, const XYZ& xyz);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Real myRadius;
		Real myFromZ;
		Real myToZ;
		bool myColor;
		bool myAxes;
		bool myLegend;
		bool myAutofit;
		Matrix3by3 myMat;
	};

	//____________________________________________________________________________INLINES
}
#endif /* XFIGWRITER_H */
