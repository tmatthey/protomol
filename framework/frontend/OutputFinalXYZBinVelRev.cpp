#include "OutputFinalXYZBinVelRev.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "XYZBinRevWriter.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputFinalXYZBinVelRev
	const string OutputFinalXYZBinVelRev::keyword("finXYZBinVelRevFile");

	OutputFinalXYZBinVelRev::OutputFinalXYZBinVelRev(): Output(1), myFilename("")
	{
	}

	OutputFinalXYZBinVelRev::OutputFinalXYZBinVelRev(const string& filename): Output(1), myFilename(filename)
	{
	}

	void OutputFinalXYZBinVelRev::doFinalize(int)
	{
		XYZBinRevWriter writer;
		if (!writer.open(myFilename))
			report << error << "Can't open " << getId() << " \'" << myFilename << "\'." << endr;
		if (!writer.write(*myVelocities))
			report << error << "Could not write " << getId() << " \'" << myFilename << "\'." << endr;
	}

	Output* OutputFinalXYZBinVelRev::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputFinalXYZBinVelRev(values[0]));
	}

	void OutputFinalXYZBinVelRev::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myFilename, ConstraintValueType::NotEmpty())));
	}
}
