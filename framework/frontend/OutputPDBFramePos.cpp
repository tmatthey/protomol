#include "OutputPDBFramePos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "PDBWriter.h"
#include "inputValueDefinitions.h"

#include <algorithm>
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//________________________________________________________ OutputPDBFramePos
	const string OutputPDBFramePos::keyword("PDBPosFile");

	OutputPDBFramePos::OutputPDBFramePos(): Output(), myFilename(""), myMinimalImage(false)
	{
	}

	OutputPDBFramePos::OutputPDBFramePos(const string& filename, int freq, bool minimal): Output(freq), myFilename(filename), myMinimalImage(minimal)
	{
	}

	void OutputPDBFramePos::doRun(int step)
	{
		const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);

		const vector<PDB::Atom>& pdbAtoms = myCache->pdb();
		if (pdbAtoms.size() != myPositions->size())
		{
			if (pdbAtoms.empty())
				report << error << " " << getId() << " can't write PDB frame without PDB atoms." << endr;
			report << error << " Number of PDB atoms (" << pdbAtoms.size() << ") does not correspond to number of atoms (" << myPositions->size() << ")." << endr;
		}
		unsigned int n = std::max(toString(getFirstStep()).size(), toString(getLastStep()).size());
		string filename = myFilename + "." + getEnd(string(n, '0') + toString(step), n) + ".pdb";
		PDBWriter pdb(filename);
		pdb.setComment(keyword + ": time=" + toString(myCache->time()) + "[fs], step : " + toString(step) + (myMinimalImage ? ", minimal Image" : ""));
		if (!pdb)
			report << error << " Can not open \'" << filename << "\' for " << getId() << "." << endr;

		if (!pdb.write(*pos, pdbAtoms))
			report << error << " Can not write frame \'" << filename << "\' for " << getId() << "." << endr;
		//report << debug << (*myPositions)[0]<<endr;
	}

	Output* OutputPDBFramePos::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputPDBFramePos(values[0], values[1], values[2]));
	}

	void OutputPDBFramePos::getParameters(vector<Parameter>& parameter) const
	{
		parameter.push_back(Parameter(getId(), Value(myFilename, ConstraintValueType::NotEmpty())));
		parameter.push_back(Parameter(keyword + "OutputFreq", Value(myOutputFreq, ConstraintValueType::Positive()), 1));
		parameter.push_back(Parameter(keyword + "MinimalImage", Value(myMinimalImage), Text("whether the coordinates should be transformed to minimal image or not")));
	}

	bool OutputPDBFramePos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
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
