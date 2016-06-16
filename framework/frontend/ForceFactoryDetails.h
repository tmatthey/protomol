/* -*- c++ -*- */
#ifndef FORCEFACTORYDETAILS_H
#define FORCEFACTORYDETAILS_H

#include "Factory.h"

#include "Force.h"
#include "Value.h"

namespace ProtoMol
{
	class ForceFactoryDetails;

	//_____________________________________________________ FactoryTraits<Force>
	template <>
	class FactoryTraits<Force>
	{
	public:
		typedef ForceFactoryDetails Details;
	};


	//_____________________________________________________ ForceFactoryDetails
	class ForceFactoryDetails : public FactoryBase<Force>
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Types
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//_____________________________________________________ ForceType
		struct ForceType
		{
			std::set<std::string, ltstrNocase> policy;
			std::set<std::string, ltstrNocase> policies;
			std::map<std::string, std::string, ltstrNocase> policiesSorted;
		};

	private:
		typedef Factory<Force> TFactory;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		ForceFactoryDetails();
		virtual ~ForceFactoryDetails();
	private:
		ForceFactoryDetails(const ForceFactoryDetails&);
		ForceFactoryDetails& operator=(const ForceFactoryDetails&);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From FactoryBase<Force>
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual std::string doPrint() const;
		virtual void doRegisterHelpText() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DetailsForceFactory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static Force* make(std::string& errMsg, const std::string& id, std::vector<Value> values = std::vector<Value>());
	private:
		Force* doMake(std::string& errMsg, const std::string& id, std::vector<Value> values) const;

	private:
		void updateCache() const;
		std::vector<std::string> splitForceString(const std::string& id) const;
		std::vector<std::string> splitForceStringSorted(const std::string& id) const;
		std::string sortForceString(const std::string& id) const;
		std::string uniqueForceString(const std::string& id) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		mutable std::map<std::string, ForceType, ltstrNocase> myForceTypes;
		mutable std::map<std::string, std::string, ltstrNocase> myForceTypesSorted;
		mutable CompareForce* myLastCompareForce;
	};
}
#endif /* FORCEFACTORYDETAILS_H */
