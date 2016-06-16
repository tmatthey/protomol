/*  -*- c++ -*-  */
#ifndef MAKEABLE_H
#define MAKEABLE_H

#include "MakeableDefinition.h"
#include "CheckpointInputStream.h"
#include "CheckpointOutputStream.h"

#include <iostream>
using namespace std;

namespace ProtoMol
{
	class Configuration;

	//________________________________________________________ Makeable
	/** 
	    Base class of all object, which can be create dynamically based on a prototype,
	    normally used together with a Factory. 
	*/

	class Makeable
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual ~Makeable()
		{
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	public:
		//virtual Makeable* make(std::string& errMsg, std::vector<Value> values)const=0;
		/// retrieve all parameters
		virtual void getParameters(std::vector<Parameter>& parameters) const =0;
		/// number of parameters
		virtual unsigned int getParameterSize() const =0;
		/// retrieve all parameters
		std::vector<Parameter> getParameters() const;
		virtual MakeableDefinition getDefinition() const;

		virtual void restoreState(CheckpointInputStream&)
		{
		};

		virtual void saveState(CheckpointOutputStream&)
		{
		};

		std::string getId() const;
		std::string getAlias() const;
		std::string setAlias(const std::string& id);
		virtual std::string getIdNoAlias() const =0;
		virtual std::string getScope() const =0;

		virtual std::string getText() const
		{
			return std::string();
		}

		bool checkParameters(std::string& errMsg, const std::vector<Value>& values) const;
		bool checkParameters(const std::vector<Value>& values) const;
		bool checkParameterTypes(const std::vector<Value>& values) const;

		virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* /*config*/) const
		{
			return checkParameterTypes(values);
		}

		template <typename T>
		static T* copy(T* obj)
		{
			T* clone = NULL;
			if (obj != NULL)
			{
				std::string err;
				std::vector<Parameter> p;
				obj->getParameters(p);
				std::vector<Value> v(p.size());
				for (unsigned int i = 0; i < p.size(); i++)
					v[i] = p[i].value;
				clone = obj->make(err, v);
			}
			return clone;
		}

	protected:
		template <typename T>
		T* adjustAlias(T* obj) const
		{
			if (obj != NULL)
				obj->setAlias(getId());
			return obj;
		}

	private:
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		std::string myAlias;
	};
}
#endif
