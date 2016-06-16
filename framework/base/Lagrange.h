/*  -*- c++ -*-  */
#ifndef LAGRANGE_H
#define LAGRANGE_H

#include "Real.h"
#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ Lagrange
	/**
	 * Lagrange interpolation.
	 * theta are the lk, where dTheta are the derivatives of lk
	 */
	class Lagrange
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Lagrange();
		Lagrange(unsigned int order);
		Lagrange(unsigned int order, Real w);
		~Lagrange();
		Lagrange(const Lagrange& Lagrange);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Lagrange
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		unsigned int getOrder() const
		{
			return myInterOrder;
		}

		void setOrder(unsigned int order);
		void set(Real w);

		static bool isSigma(unsigned int order)
		{
			return (Lagrange(order, 0.0).theta[order - 1] == 0.0 && Lagrange(order, 0.0).dTheta[order - 1] == 0.0);
		}

		///< true iff one theta is 1 and all other 0 for w=0

		static const std::string& getKeyword()
		{
			return keyword;
		}

		///< Returns the keyword/name of this interpolation

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		unsigned int myInterOrder;
	public:
		Real* theta; ///< lk
		Real* dTheta; ///< derivatives of lk
	public:
		static const std::string keyword;
	};

	//______________________________________________________________________ INLINES
}
#endif /* LAGRANGE_H */
