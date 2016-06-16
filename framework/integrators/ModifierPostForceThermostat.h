/*  -*- c++ -*-  */
#ifndef MODIFIERPOSTFORCETHERMOSTAT_H
#define MODIFIERPOSTFORCETHERMOSTAT_H

#include "Modifier.h"

namespace ProtoMol
{
	//_________________________________________________________________ ModifierPostForceThermostat
	template <class TIntegrator>
	class ModifierPostForceThermostat : public Modifier
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierPostForceThermostat(TIntegrator* i): Modifier(Constant::MAX_INT - 400), myTheIntegrator(i)
		{
		}

		ModifierPostForceThermostat(TIntegrator* i, int order): Modifier(order), myTheIntegrator(i)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Modifier
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual bool isInternal() const
		{
			return true;
		}

	private:
		virtual void doExecute()
		{
			myTheIntegrator->PostForceThermostat();
		}

		virtual std::string doPrint() const
		{
			return std::string("PostForceThermostat");
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		TIntegrator* myTheIntegrator;
	};
}
#endif /* MODIFIERPOSTFORCETHERMOSTAT_H */
