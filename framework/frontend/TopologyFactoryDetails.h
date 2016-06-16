/* -*- c++ -*- */
#ifndef TOPOLOGYFACTORYDETAILS_H
#define TOPOLOGYFACTORYDETAILS_H

#include "Factory.h"
#include "GenericTopology.h"

namespace ProtoMol
{
	class TopologyFactoryDetails;
	class Configuration;

	//_____________________________________________________ FactoryTraits<GenericTopology>
	template <>
	class FactoryTraits<GenericTopology>
	{
	public:
		typedef TopologyFactoryDetails Details;
	};

	//_____________________________________________________ TopologyFactoryDetails
	class TopologyFactoryDetails : public FactoryBase<GenericTopology>
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Types
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		typedef Factory<GenericTopology> TFactory;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		TopologyFactoryDetails();
		virtual ~TopologyFactoryDetails();
	private:
		TopologyFactoryDetails(const TopologyFactoryDetails&);
		TopologyFactoryDetails& operator=(const TopologyFactoryDetails&);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From FactoryBase<GenericTopology>
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual std::string doPrint() const;
		virtual void doRegisterHelpText() const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class DetailsTopologyFactory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static void registerAllExemplarsConfiguration(Configuration* config);
		static GenericTopology* make(std::string& errMsg, const Configuration* config);
		static GenericTopology* make(std::string& errMsg, const std::string& id, const std::vector<Value>& values);

	private:
		void doRegisterAllExemplarsConfiguration(Configuration* config) const;
		GenericTopology* doMake(std::string& errMsg, const std::string& id, const std::vector<Value>& values) const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	};
}
#endif /* TOPOLOGYFACTORYDETAILS_H */
