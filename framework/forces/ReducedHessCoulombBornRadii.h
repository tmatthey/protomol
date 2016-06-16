/*  -*- c++ -*-  */
#ifndef REDUCEDHESSCOULOMBBORNRADII_H
#define REDUCEDHESSCOULOMBBORNRADII_H

#include "Matrix3by3.h"
#include "ExclusionTable.h"
#include "Vector3DBlock.h"

namespace ProtoMol
{
	class GenericTopology;
	class CoulombBornRadiiForce;

	class ReducedHessCoulombBornRadii
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  Constructor,destructor
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class ReducedHessCoulombBornRadii
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
		                      const Vector3DBlock* positions,
		                      CoulombBornRadiiForce& hForce
		) const;
	};
}
#endif
