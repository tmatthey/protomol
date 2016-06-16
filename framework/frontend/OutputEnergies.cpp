#include "OutputEnergies.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "OutputCache.h"
#include "inputValueDefinitions.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::setw;
using std::endl;
using std::flush;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;

namespace ProtoMol
{
	//________________________________________________________ Output
	const string OutputEnergies::keyword("allEnergiesFile");

	OutputEnergies::OutputEnergies(): OutputFile(), myDoMolecularTemperature(false), myDoShadow(false)
	{
	}

	OutputEnergies::OutputEnergies(const string& filename, int freq, int cacheFreq, int cacheSize, Real closeTime, bool doMolTemp, bool doShadow): OutputFile(filename, freq, cacheFreq, cacheSize, closeTime), myDoMolecularTemperature(doMolTemp), myDoShadow(doShadow)
	{
	}

	void OutputEnergies::doInitialize()
	{
		ofstream allEnergiesHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
		if (!allEnergiesHeaderFile)
			report << error << " Can not open \'" << myFilename << ".header\' for " << getId() << "." << endr;

		allEnergiesHeaderFile << setw(14)
			<< "Time(fs)" << " "
			<< setw(14)
			<< "E_potential" << " "
			<< setw(14)
			<< "E_kinetic" << " "
			<< setw(14)
			<< "E_total" << " "
			<< setw(14)
			<< "Temperature" << " "
			<< setw(14)
			<< "E_bond" << " "
			<< setw(14)
			<< "E_angle" << " "
			<< setw(14)
			<< "E_dihedral" << " "
			<< setw(14)
			<< "E_improper" << " "
			<< setw(14)
			<< "E_VdW" << " "
			<< setw(14)
			<< "E_coulomb" << " "
			<< setw(14)
			<< "E_other" << " "
			<< setw(14)
			<< "Volume(A^3)";
		if (myEnergies->virial())
			allEnergiesHeaderFile << " "
				<< setw(14) << "Pressure(bar)";
		if (myEnergies->molecularVirial())
			allEnergiesHeaderFile << " "
				<< setw(14) << "Mol_Pres(bar)";
		if (myDoMolecularTemperature)
			allEnergiesHeaderFile << " "
				<< setw(14) << "Mol_Temp(K)";
		if (myDoShadow)
			allEnergiesHeaderFile << " "
				<< setw(20) << "E_shadow";

		allEnergiesHeaderFile << endl;
		allEnergiesHeaderFile.close();
		open();
		close();
	}

	void OutputEnergies::doRunCached(int)
	{
		myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setw(14)
			<< setprecision(2)
			<< setiosflags(std::ios::showpoint | std::ios::fixed)
			<< myCache->time() << " "
			<< resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setiosflags(std::ios::floatfield)
			<< setprecision(8)
			<< setw(14)
			<< myCache->potentialEnergy() << " "
			<< setw(14)
			<< myCache->kineticEnergy() << " "
			<< setw(14)
			<< myCache->totalEnergy() << " "
			<< setw(14)
			<< myCache->temperature() << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::BOND] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::ANGLE] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::DIHEDRAL] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::IMPROPER] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::LENNARDJONES] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::COULOMB] << " "
			<< setw(14)
			<< (*myEnergies)[ScalarStructure::OTHER] << " "
			<< setw(14)
			<< myCache->volume();
		if (myEnergies->virial())
			myBuffer << " " << setw(14)
				<< myCache->pressure();
		if (myEnergies->molecularVirial())
			myBuffer << " " << setw(14)
				<< myCache->molecularPressure();
		if (myDoMolecularTemperature)
			myBuffer << " " << setw(14)
				<< myCache->molecularTemperature();
		if (myDoShadow)
			myBuffer << " " << setw(20)
				<< setprecision(16) //  High precision needed.
				<< (*myEnergies)[ScalarStructure::SHADOW];

		myBuffer << endl;
	}

	Output* OutputEnergies::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputEnergies(values[0], values[1], values[2], values[3], values[4], values[5], values[6]));
	}

	void OutputEnergies::getParameters(vector<Parameter>& parameter) const
	{
		OutputFile::getParameters(parameter);
		parameter.push_back(Parameter("molecularTemperature", Value(myDoMolecularTemperature), false));
		parameter.push_back(Parameter("shadowEnergy", Value(myDoShadow), false));
	}
}
