#include "Force.h"
#include "CompareForce.h"

using namespace ProtoMol::Report;

using std::vector;
using std::string;

namespace ProtoMol
{
	//_________________________________________________________________ Force

	const string Force::scope("Force");

	Force* Force::make(string& errMsg, const vector<Value>& values) const
	{
		errMsg = "";
		if (!checkParameters(errMsg, values))
			return NULL;
		return adjustAlias(doMake(errMsg, values));
	}

	CompareForce* Force::makeCompareForce(Force* actualForce, CompareForce* compareForce) const
	{
		report << error << "No support to compare " << (actualForce ? actualForce->getId() : string("null pointer")) << " vs. " << (compareForce ? compareForce->getId() : string("null pointer")) << "." << endr;
		return NULL;
	}

	TimeForce* Force::makeTimeForce(Force* actualForce) const
	{
		report << error << "No support to time " << (actualForce ? actualForce->getId() : string("null pointer")) << "." << endr;
		return NULL;
	}

	void Force::setParameters(string& errMsg, vector<Value> values)
	{
		errMsg = "";
		if (checkParameters(errMsg, values))
			doSetParameters(errMsg, values);
	}
}
