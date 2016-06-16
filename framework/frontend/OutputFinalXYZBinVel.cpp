#include "OutputFinalXYZBinVel.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "XYZBinWriter.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputFinalXYZBinVel
	const string OutputFinalXYZBinVel::keyword("finXYZBinVelFile");

	OutputFinalXYZBinVel::OutputFinalXYZBinVel(): Output(1), myFilename("")
	{
	}

	OutputFinalXYZBinVel::OutputFinalXYZBinVel(const string& filename): Output(1), myFilename(filename)
	{
	}

	void OutputFinalXYZBinVel::doFinalize(int)
	{
		XYZBinWriter writer;
		if (!writer.open(myFilename))
			report << error << "Can't open " << getId() << " \'" << myFilename << "\'." << endr;
		if (!writer.write(*myVelocities))
			report << error << "Could not write " << getId() << " \'" << myFilename << "\'." << endr;
	}

	Output* OutputFinalXYZBinVel::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputFinalXYZBinVel(values[0]));
	}

	void OutputFinalXYZBinVel::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myFilename, ConstraintValueType::NotEmpty())));
	}
}
