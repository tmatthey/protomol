/*  -*- c++ -*-  */
#ifndef NONSTANDARDINTEGRATOR_H
#define NONSTANDARDINTEGRATOR_H

#include "Integrator.h"

namespace ProtoMol
{
	//_________________________________________________________________ NonStandardIntegrator
	class ForceGroup;

	class NonStandardIntegrator: public Integrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NonStandardIntegrator()
		{
		};

		NonStandardIntegrator(ForceGroup* forceGroup): Integrator(forceGroup)
		{
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NonStandardIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual NonStandardIntegrator* make(std::string& errMsg, std::vector<Value> values, ForceGroup* fg) const =0;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
	};

	//______________________________________________________________________ INLINES
}
#endif
