#include "CutoffSwitchingFunction.h"
#include "stringutilities.h"
using std::string;
using std::vector;

namespace ProtoMol
{
	//_________________________________________________________________ CutoffSwitchingFunction
	CutoffSwitchingFunction::CutoffSwitchingFunction(): myCutoff(0.0), myCutoff2(0.0)
	{
	}

	CutoffSwitchingFunction::CutoffSwitchingFunction(Real cutoff): myCutoff(cutoff), myCutoff2(cutoff * cutoff)
	{
	}

	void CutoffSwitchingFunction::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("cutoff swf cutoff")));
	}

	CutoffSwitchingFunction CutoffSwitchingFunction::make(string&, vector<Value> values)
	{
		return CutoffSwitchingFunction(values[0]);
	}
}
