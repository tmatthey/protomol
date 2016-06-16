#include "OutputXYZBinTrajectoryPos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZBinWriter.h"
#include "inputValueDefinitions.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputXYZBinTrajectoryPos
	const string OutputXYZBinTrajectoryPos::keyword("XYZBinPosFile");

	OutputXYZBinTrajectoryPos::OutputXYZBinTrajectoryPos(): Output(), myXYZ(NULL), myMinimalImage(false)
	{
	}

	OutputXYZBinTrajectoryPos::OutputXYZBinTrajectoryPos(const string& filename, int freq, bool minimal): Output(freq), myXYZ(new XYZBinWriter(filename)), myMinimalImage(minimal)
	{
	}

	OutputXYZBinTrajectoryPos::~OutputXYZBinTrajectoryPos()
	{
		if (myXYZ != NULL)
			delete myXYZ;
	}

	void OutputXYZBinTrajectoryPos::doInitialize()
	{
		if (myXYZ == NULL || !myXYZ->open())
			report << error << " Can not open \'" << (myXYZ != NULL ? myXYZ->getFilename() : "") << "\' for " << getId() << "." << endr;
	}

	void OutputXYZBinTrajectoryPos::doRun(int)
	{
		//const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
		//no minimal image for a binary restart file
		const Vector3DBlock* pos = myPositions;

		if (!myXYZ->write(*pos))
			report << error << "Could not write " << getId() << " \'" << myXYZ->getFilename() << "\'." << endr;
	}

	void OutputXYZBinTrajectoryPos::doFinalize(int)
	{
		myXYZ->close();
	}

	Output* OutputXYZBinTrajectoryPos::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputXYZBinTrajectoryPos(values[0], values[1], values[2]));
	}

	void OutputXYZBinTrajectoryPos::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myXYZ != NULL ? myXYZ->getFilename() : "", ConstraintValueType::NotEmpty())));
		parameter.push_back(Parameter(keyword + "OutputFreq", Value(myOutputFreq, ConstraintValueType::Positive())));
		parameter.push_back(Parameter(keyword + "MinimalImage", Value(myMinimalImage), Text("whether the coordinates should be transformed to minimal image or not")));
	}

	bool OutputXYZBinTrajectoryPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
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
