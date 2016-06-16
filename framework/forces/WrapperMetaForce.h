/* -*- c++ -*- */
#ifndef WRAPPERMETAFORCE_H
#define WRAPPERMETAFORCE_H

#include "MetaForce.h"

namespace ProtoMol
{
	class Force;

	//_________________________________________________________________ WrapperMetaForce

	class WrapperMetaForce : public MetaForce
	{
		// This class contains the definition of one Meta force 

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		WrapperMetaForce(std::string id, bool minimal, Force* force, std::string forceAlias);
		WrapperMetaForce(std::string id, bool minimal, Force* force1, std::string forceAlias1, Force* force2, std::string forceAlias2);
		virtual ~WrapperMetaForce();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const
		{
			parameters = myParameters;
		}

		virtual unsigned int getParameterSize() const
		{
			return myParameterSize;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Force
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getKeyword() const
		{
			return keyword;
		}

	private:
		virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class WrapperMetaForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void getDeepForces(std::vector<Force*>& forces) const;
	private:
		void addForce(Force* force, std::string forceAlias);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
	public:
		const std::string keyword;
	private:
		std::vector<Force*> myForces;
		bool myMinimal;
		std::vector<Parameter> myParameters;
		unsigned int myParameterSize;
	};

	//______________________________________________________________________ INLINES
}
#endif /* WRAPPERMETAFORCE_H */
