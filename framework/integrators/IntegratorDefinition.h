/*  -*- c++ -*-  */
#ifndef INTEGRATORDEFINITION_H
#define INTEGRATORDEFINITION_H

#include "MakeableDefinition.h"

namespace ProtoMol
{
	//________________________________________________________ IntegratorDefinition
	struct IntegratorDefinition
	{
		// Container struct for integrator definitions

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		IntegratorDefinition()
		{
		}

		IntegratorDefinition(const MakeableDefinition& i, const std::vector<MakeableDefinition>& f):
			integrator(i), forces(f)
		{
		}

		std::string print() const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		MakeableDefinition integrator;
		std::vector<MakeableDefinition> forces;
	};
}
#endif /* INTEGRATORDEFINITION_H */
