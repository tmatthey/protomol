/*  -*- c++ -*-  */
#ifndef OBJECTDEFINITION_H
#define OBJECTDEFINITION_H

#include <vector>
#include "Parameter.h"

namespace ProtoMol
{
	//________________________________________________________ ObjectDefinition
	struct ObjectDefinition
	{
		// Container struct for object definitions

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ObjectDefinition()
		{
		}

		ObjectDefinition(const std::string& i, const std::vector<Parameter>& p): id(i), parameters(p)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		std::string id;
		std::vector<Parameter> parameters;
	};
}
#endif /* OBJECTDEFINITION_H */
