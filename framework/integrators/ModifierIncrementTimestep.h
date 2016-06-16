/*  -*- c++ -*-  */
#ifndef MODIFIERINCREMENTTIMESTEP_H
#define MODIFIERINCREMENTTIMESTEP_H
#include "Modifier.h"
#include "GenericTopology.h"
#include "STSIntegrator.h"

namespace ProtoMol
{
	class ModifierIncrementTimestep : public Modifier
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierIncrementTimestep(STSIntegrator* i): Modifier(Constant::MAX_INT), myTheIntegrator(i)
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
			myTopology->time += myTheIntegrator->getTimestep();
		}

		virtual std::string doPrint() const
		{
			return std::string("IncrementTimestep");
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		STSIntegrator* myTheIntegrator;
	};
}
#endif /* MODIFIERINCREMENTTIMESTEP_H */
