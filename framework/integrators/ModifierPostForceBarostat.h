/*  -*- c++ -*-  */
#ifndef MODIFIERPOSTFORCEBAROSTAT_H
#define MODIFIERPOSTFORCEBAROSTAT_H

#include "Modifier.h"


namespace ProtoMol
{
	//_________________________________________________________________ ModifierPostForceBarostat
	template <class TIntegrator>
	class ModifierPostForceBarostat : public Modifier
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierPostForceBarostat(TIntegrator* i): Modifier(Constant::MAX_INT - 400), myTheIntegrator(i)
		{
		}

		ModifierPostForceBarostat(TIntegrator* i, int order): Modifier(order), myTheIntegrator(i)
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
			myTheIntegrator->PostForceBarostat();
		}

		virtual std::string doPrint() const
		{
			return std::string("PostForceBarostat");
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		TIntegrator* myTheIntegrator;
	};
}
#endif /* MODIFIERPOSTFORCEBAROSTAT_H */
