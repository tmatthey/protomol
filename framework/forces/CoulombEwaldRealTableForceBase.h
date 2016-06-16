/*  -*- c++ -*-  */
#ifndef COULOMBEWALDREALTABLEFORCEBASE_H
#define COULOMBEWALDREALTABLEFORCEBASE_H

#include<string>
#include "mathutilities.h"

namespace ProtoMol
{
	//_________________________________________________________________ CoulombEwaldRealTableForceBase

	class CoulombEwaldRealTableForceBase
	{
	public:
		//_________________________________________________________________ LookUpValues
		/**
		 * Defines the function(s) to be tabulated.
		 */
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
			LookUpValues(Real a): myAlpha(a)
			{
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// New methods of class LookUpValues
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		public:
			template <typename TReal>
			void assign(Real r, Real r1, int ex, Real v, Real d, TReal* val) const
			{
				Real r2 = (ex == 2 ? r : power<2>(r1));
				Real e = erfc(myAlpha * r1) / r1;
				Real f = (e + 2.0 * myAlpha / sqrt(Constant::M_PI) * exp(-myAlpha * myAlpha * r2)) / r2;
				val[0] = e * v;
				val[1] = -0.5 * (f * v - e * d);
			}

		private:
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// My data members
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Real myAlpha;
		};

	public:
		static const std::string keyword;
	};
}
#endif /* COULOMBEWALDREALTABLEFORCEBASE_H */
