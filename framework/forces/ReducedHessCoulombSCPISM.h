/*  -*- c++ -*-  */
#ifndef REDUCEDHESSCOULOMBSCPISM_H
#define REDUCEDHESSCOULOMBSCPISM_H

#include "Matrix3by3.h"
#include "ExclusionTable.h"

namespace ProtoMol
{
	class GenericTopology;

	class ReducedHessCoulombSCPISM
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class ReducedHessCoulombSCPISM
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
		                      ExclusionClass excl,
		                      Real alpha_ij,
		                      Real es) const;
	};
}
#endif
