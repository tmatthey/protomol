#include "OutputDCDTrajectoryNoWater.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "DCDTrajectoryWriter.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::cout;

namespace ProtoMol
{
	//________________________________________________________ OutputDCDTrajectoryNoWater
	const string OutputDCDTrajectoryNoWater::keyword("DCDFileNoWater");

	OutputDCDTrajectoryNoWater::OutputDCDTrajectoryNoWater(): Output(), myDCD(NULL), myMinimalImage(false)
	{
	}

	OutputDCDTrajectoryNoWater::OutputDCDTrajectoryNoWater(const string& filename, int freq, bool minimal): Output(freq), myDCD(new DCDTrajectoryWriter(filename)), myMinimalImage(minimal)
	{
	}

	OutputDCDTrajectoryNoWater::~OutputDCDTrajectoryNoWater()
	{
		if (myDCD != NULL)
			delete myDCD;
	}

	void OutputDCDTrajectoryNoWater::doInitialize()
	{
		if (myDCD == NULL || !myDCD->open())
			report << error << " Can not open \'" << (myDCD != NULL ? myDCD->getFilename() : "") << "\' for " << getId() << "." << endr;
	}

	void OutputDCDTrajectoryNoWater::doRun(int)
	{
		//report << "Going to call the cache" << endr;
		const Vector3DBlock* pos = myCache->PositionsNoWater();
		//cout << "Cache called, Vector3DBlock size is: " << (*pos).size() << endl;
		//const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
		if (!myDCD->write(*pos))
			report << error << "Could not write " << getId() << " \'" << myDCD->getFilename() << "\'." << endr;
	}

	void OutputDCDTrajectoryNoWater::doFinalize(int)
	{
		myDCD->close();
	}

	Output* OutputDCDTrajectoryNoWater::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputDCDTrajectoryNoWater(values[0], values[1], values[2]));
	}

	void OutputDCDTrajectoryNoWater::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myDCD != NULL ? myDCD->getFilename() : "", ConstraintValueType::NotEmpty())));
		parameter.push_back(Parameter(keyword + "OutputFreq", Value(myOutputFreq, ConstraintValueType::Positive())));
		parameter.push_back(Parameter(keyword + "MinimalImage", Value(myMinimalImage), Text("whether the coordinates should be transformed to minimal image or not")));
	}

	bool OutputDCDTrajectoryNoWater::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
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
