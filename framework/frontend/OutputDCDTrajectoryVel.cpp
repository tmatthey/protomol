#include "OutputDCDTrajectoryVel.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "DCDTrajectoryWriter.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputDCDTrajectoryVel
	const string OutputDCDTrajectoryVel::keyword("DCDVELFile");

	OutputDCDTrajectoryVel::OutputDCDTrajectoryVel(): Output(), myDCD(NULL), myMinimalImage(false)
	{
	}

	OutputDCDTrajectoryVel::OutputDCDTrajectoryVel(const string& filename, int freq, bool minimal): Output(freq), myDCD(new DCDTrajectoryWriter(filename)), myMinimalImage(minimal)
	{
	}

	OutputDCDTrajectoryVel::~OutputDCDTrajectoryVel()
	{
		if (myDCD != NULL)
			delete myDCD;
	}

	void OutputDCDTrajectoryVel::doInitialize()
	{
		if (myDCD == NULL || !myDCD->open())
			report << error << " Can not open \'" << (myDCD != NULL ? myDCD->getFilename() : "") << "\' for " << getId() << "." << endr;
	}

	void OutputDCDTrajectoryVel::doRun(int)
	{
		const Vector3DBlock* vel = myVelocities;
		if (!myDCD->write(*vel))
			report << error << "Could not write " << getId() << " \'" << myDCD->getFilename() << "\'." << endr;
	}

	void OutputDCDTrajectoryVel::doFinalize(int)
	{
		myDCD->close();
	}

	Output* OutputDCDTrajectoryVel::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputDCDTrajectoryVel(values[0], values[1], values[2]));
	}

	void OutputDCDTrajectoryVel::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myDCD != NULL ? myDCD->getFilename() : "", ConstraintValueType::NotEmpty())));
		parameter.push_back(Parameter(keyword + "OutputFreq", Value(myOutputFreq, ConstraintValueType::Positive())));
		parameter.push_back(Parameter(keyword + "MinimalImage", Value(myMinimalImage), Text("whether the coordinates should be transformed to minimal image or not (NA for Vel)")));
	}

	bool OutputDCDTrajectoryVel::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
	{
		if (!checkParameterTypes(values))
			return false;
		if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
			values[1] = (*config)[InputOutputfreq::keyword];
		if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
			values[2] = (*config)[InputMinimalImage::keyword];
		return checkParameters(values);
	}
}
