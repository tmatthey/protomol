#include "ShiftSwitchingFunction.h"
using std::string;
using std::vector;

namespace ProtoMol
{
	//_________________________________________________________________ ShiftSwitchingFunction
	ShiftSwitchingFunction::ShiftSwitchingFunction(): myCutoff(0.0),
	                                                  myCutoff2(0.0),
	                                                  myCutoff_2(0.0),
	                                                  my4Cutoff_2(0.0)
	{
	}

	ShiftSwitchingFunction::ShiftSwitchingFunction(Real cutoff): myCutoff(cutoff),
	                                                             myCutoff2(cutoff * cutoff),
	                                                             myCutoff_2(1.0 / (cutoff * cutoff)),
	                                                             my4Cutoff_2(-4.0 / (cutoff * cutoff))
	{
	}

	void ShiftSwitchingFunction::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("shift swf cutoff")));
	}

	ShiftSwitchingFunction ShiftSwitchingFunction::make(string&, vector<Value> values)
	{
		return ShiftSwitchingFunction(values[0]);
	}
}
