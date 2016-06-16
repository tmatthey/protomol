/*  -*- c++ -*-  */
#ifndef MODIFIERFRICTION_H
#define MODIFIERFRICTION_H

#include "Modifier.h"
#include "NoseNVTLeapfrogIntegrator.h"

namespace ProtoMol
{
	//_________________________________________________________________ ModifierFriction
	class ModifierFriction : public Modifier
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierFriction(NoseNVTLeapfrogIntegrator* i): Modifier(Constant::MAX_INT - 400), myTheIntegrator(i)
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
			myTheIntegrator->friction();
		}

		virtual std::string doPrint() const
		{
			return std::string("Friction");
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		NoseNVTLeapfrogIntegrator* myTheIntegrator;
	};
}
#endif /* MODIFIER_H */
