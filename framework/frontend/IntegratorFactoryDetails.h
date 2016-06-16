/* -*- c++ -*- */
#ifndef INTEGRATORFACTORYDETAILS_H
#define INTEGRATORFACTORYDETAILS_H

#include "Factory.h"

#include "Integrator.h"
#include "Value.h"

namespace ProtoMol
{
	class IntegratorFactoryDetails;

	//_____________________________________________________ FactoryTraits<Integrator>
	template <>
	class FactoryTraits<Integrator>
	{
	public:
		typedef IntegratorFactoryDetails Details;
	};


	//_____________________________________________________ IntegratorFactoryDetails
	class IntegratorFactoryDetails : public FactoryBase<Integrator>
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Types
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		struct IntegratorInput
		{
			IntegratorInput(): prototype(NULL)
			{
			}

			const Integrator* prototype;
			std::vector<Value> values;
			std::vector<std::string> forces;
		};

		typedef Factory<Integrator> TFactory;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		IntegratorFactoryDetails();
		virtual ~IntegratorFactoryDetails();
	private:
		IntegratorFactoryDetails(const IntegratorFactoryDetails&);
		IntegratorFactoryDetails& operator=(const IntegratorFactoryDetails&);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From FactoryBase<Integrator>
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual std::string doPrint() const;
		virtual void doRegisterHelpText() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DetailsIntegratorFactory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static Integrator* make(std::string& errMsg, const std::string& definition);
	private:
		Integrator* doMake(std::string& errMsg, const std::string& definition) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	};
}
#endif /* INTEGRATORFACTORYDETAILS_H */
