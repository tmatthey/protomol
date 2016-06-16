/* -*- c++ -*- */
#ifndef SYSTEMTIMEFORCE_H
#define SYSTEMTIMEFORCE_H

#include "TimeForce.h"
#include "SystemForce.h"

namespace ProtoMol
{
	//_________________________________________________________________ SystemTimeForce

	class SystemTimeForce : public TimeForce, public SystemForce
	{
		// This class contains the definition of one force

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		SystemTimeForce(Force* actualForce);

		virtual ~SystemTimeForce()
		{
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class SystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* positions,
		                      Vector3DBlock* forces,
		                      ScalarStructure* energies);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              Vector3DBlock* forces,
		                              ScalarStructure* energies);
	};

	//______________________________________________________________________ INLINES
}
#endif /* SYSTEMTIMEFORCE_H */
