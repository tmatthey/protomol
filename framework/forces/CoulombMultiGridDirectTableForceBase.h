/*  -*- c++ -*-  */
#ifndef COULOMBMULTIGRIDDIRECTTABLEFORCEBASE_H
#define COULOMBMULTIGRIDDIRECTTABLEFORCEBASE_H

#include<string>
#include "mathutilities.h"

namespace ProtoMol
{
	//_________________________________________________________________ CoulombMultiGridDirectTableForceBase

	class CoulombMultiGridDirectTableForceBase
	{
	public:
		//_________________________________________________________________ LookUpValues
		/**
		 * Defines the function(s) to be tabulated.
		 */
		template <class TKernel>
		class LookUpValues
		{
		public:
			enum
			{
				ENTRIES=1
			};

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Constructors, destructors, assignment
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			LookUpValues(Real a): myS(a)
			{
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// New methods of class LookUpValues
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			template <typename TReal>
			void assign(Real, Real r1, int, Real v, Real d, TReal* val) const
			{
				Real e = (TKernel::kernel(r1) - TKernel::smooth(r1, myS, 1.0 / myS));
				Real f = (-TKernel::dKernel(r1) + TKernel::dSmooth(r1, myS, 1.0 / myS)) / r1;
				val[0] = e * v;
				val[1] = -0.5 * (f * v - e * d);
			}

		private:
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My data members
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Real myS;
		};

	public:
		static const std::string keyword;
	};
}
#endif /* COULOMBMULTIGRIDDIRECTTABLEFORCEBASE_H */
