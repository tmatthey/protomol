/* -*- c++ -*- */
#ifndef EXTENDEDCOMPAREFORCE_H
#define EXTENDEDCOMPAREFORCE_H

#include "CompareForce.h"
#include "ExtendedForce.h"

namespace ProtoMol
{
	//_________________________________________________________________ ExtendedCompareForce

	class ExtendedCompareForce : public CompareForce, public ExtendedForce
	{
		// This class contains the definition of one force

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ExtendedCompareForce(Force* actualForce, CompareForce* compareForce);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class ExtendedForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* positions,
		                      const Vector3DBlock* velocities,
		                      Vector3DBlock* forces,
		                      ScalarStructure* energies);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              const Vector3DBlock* velocities,
		                              Vector3DBlock* forces,
		                              ScalarStructure* energies);
	};

	//______________________________________________________________________ INLINES
}
#endif /* EXTENDEDCOMPAREFORCE_H */
