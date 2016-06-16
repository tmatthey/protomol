/* -*- c++ -*- */
#ifndef OUTPUTFACTORYDETAILS_H
#define OUTPUTFACTORYDETAILS_H

#include "Factory.h"

#include "Output.h"
#include "Value.h"

namespace ProtoMol
{
	class OutputFactoryDetails;
	class OutputCollection;

	//_____________________________________________________ FactoryTraits<Output>
	template <>
	class FactoryTraits<Output>
	{
	public:
		typedef OutputFactoryDetails Details;
	};


	//_____________________________________________________ OutputFactoryDetails
	class OutputFactoryDetails : public FactoryBase<Output>
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Types
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		typedef Factory<Output> TFactory;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		OutputFactoryDetails();
		virtual ~OutputFactoryDetails();
	private:
		OutputFactoryDetails(const OutputFactoryDetails&);
		OutputFactoryDetails& operator=(const OutputFactoryDetails&);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From FactoryBase<Output>
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual std::string doPrint() const;
		virtual void doRegisterHelpText() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DetailsOutputFactory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static void registerAllExemplarsConfiguration(Configuration* config);
		static Output* make(std::string& errMsg, const std::string& id, const std::vector<Value>& values);
		static OutputCollection* makeCollection(std::string& errMsg, const Configuration* config);
	private:
		void doRegisterAllExemplarsConfiguration(Configuration* config) const;
		Output* doMake(std::string& errMsg, const std::string& id, const std::vector<Value>& values) const;
		OutputCollection* doMakeCollection(std::string& errMsg, const Configuration* config) const;


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	};
}
#endif /* OUTPUTFACTORYDETAILS_H */
