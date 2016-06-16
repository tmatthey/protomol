/*  -*- c++ -*-  */
#ifndef REDUCEDHESSCOULOMB_H
#define REDUCEDHESSCOULOMB_H

#include "Matrix3by3.h"
#include "ExclusionTable.h"

namespace ProtoMol
{
	class GenericTopology;

	class ReducedHessCoulomb
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class ReducedHessCoulomb
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Matrix3by3 operator()(const Real rawEnergy,
		                      const Real rawForce,
		                      Real distSquared,
		                      Real rDistSquared,
		                      const Vector3D& diff,
		                      const GenericTopology* topo,
		                      int atom1, int atom2,
		                      const Real switchingValue,
		                      const Real switchingDeriv,
		                      const Matrix3by3& switchingHess,
		                      ExclusionClass excl) const;
	};
}
#endif
