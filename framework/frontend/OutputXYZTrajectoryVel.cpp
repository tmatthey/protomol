#include "OutputXYZTrajectoryVel.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZTrajectoryWriter.h"
#include "inputValueDefinitions.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputXYZTrajectoryVel
	const string OutputXYZTrajectoryVel::keyword("XYZVelFile");

	OutputXYZTrajectoryVel::OutputXYZTrajectoryVel(): Output(), myXYZ()
	{
	}

	OutputXYZTrajectoryVel::OutputXYZTrajectoryVel(const string& filename, int freq): Output(freq), myXYZ(new XYZTrajectoryWriter(filename))
	{
	}

	OutputXYZTrajectoryVel::~OutputXYZTrajectoryVel()
	{
		if (myXYZ != NULL)
			delete myXYZ;
	}

	void OutputXYZTrajectoryVel::doInitialize()
	{
		if (myXYZ == NULL || !myXYZ->open())
			report << error << " Can not open \'" << (myXYZ != NULL ? myXYZ->getFilename() : "") << "\' for " << getId() << "." << endr;
	}

	void OutputXYZTrajectoryVel::doRun(int)
	{
		if (!myXYZ->write(*myVelocities, myTopology->atoms, myTopology->atomTypes))
			report << error << "Could not write " << getId() << " \'" << myXYZ->getFilename() << "\'." << endr;
	}

	void OutputXYZTrajectoryVel::doFinalize(int)
	{
		myXYZ->close();
	}

	Output* OutputXYZTrajectoryVel::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputXYZTrajectoryVel(values[0], values[1]));
	}

	void OutputXYZTrajectoryVel::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myXYZ != NULL ? myXYZ->getFilename() : "", ConstraintValueType::NotEmpty())));
		parameter.push_back(Parameter(keyword + "OutputFreq", Value(myOutputFreq, ConstraintValueType::Positive())));
	}

	bool OutputXYZTrajectoryVel::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
	{
		if (!checkParameterTypes(values))
			return false;
		if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
			values[1] = (*config)[InputOutputfreq::keyword];
		return checkParameters(values);
	}
}
