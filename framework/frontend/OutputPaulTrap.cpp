#include "OutputPaulTrap.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "stringutilities.h"
#include "OutputCache.h"
#include "XYZWriter.h"
#include "Integrator.h"
#include "OutputScreen.h"
#include "inputValueDefinitions.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;

namespace ProtoMol
{
	//________________________________________________________ Output
	const string OutputPaulTrap::keyword("paulFile");

	OutputPaulTrap::OutputPaulTrap()
		: OutputFile(),
		  myOmegaR(0.0),
		  myOmegaZ(0.0),
		  myFilenameLow(""),
		  myDoPaulLow(false),
		  myScreen(false),
		  myLowOut(""),
		  myUnit("fs"),
		  myFactor(1.0)
	{
	}

	OutputPaulTrap::OutputPaulTrap(const std::string& filename, int freq, int cacheFreq, int cacheSize,
	                               Real closeTime, Real omegar, Real omegaz, const std::string& filenameLow,
	                               bool doLow, bool screen)
		: OutputFile(filename, freq, cacheFreq, cacheSize, closeTime),
		  myOmegaR(omegar),
		  myOmegaZ(omegaz),
		  myFilenameLow(filenameLow),
		  myDoPaulLow(doLow),
		  myScreen(screen),
		  myLowOut(""),
		  myUnit("fs"),
		  myFactor(1.0)
	{
	}


	void OutputPaulTrap::doInitialize()
	{
		Real step = myIntegrator->getTimestep() * std::max(1, std::min(myOutputFreq, (int)(*myConfig)[InputNumsteps::keyword]));
		if (step >= 1e13)
		{
			myUnit = "s";
			myFactor = 1e-15;
		}
		else if (step >= 1e10)
		{
			myUnit = "ms";
			myFactor = 1e-12;
		}
		else if (step >= 1e7)
		{
			myUnit = "us";
			myFactor = 1e-9;
		}
		else if (step >= 1e4)
		{
			myUnit = "ns";
			myFactor = 1e-6;
		}
		else if (step >= 1e1)
		{
			myUnit = "ps";
			myFactor = 1e-3;
		}

		ofstream paulTrapHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
		if (!paulTrapHeaderFile)
			report << error << " Can not open \'" << myFilename << ".header\' for " << getId() << "." << endr;

		paulTrapHeaderFile << setw(20)
			<< "Time(fs)" << " "
			<< setw(20)
			<< "E_total" << " "
			<< setw(20)
			<< "E_potential" << " "
			<< setw(20)
			<< "E_potential[Nq^2/a]" << " "
			<< setw(20)
			<< "E_kinetic" << " "
			<< setw(20)
			<< "E_coulomb" << " "
			<< setw(20)
			<< "E_paulTrap" << " "
			<< setw(20)
			<< "U_cohesive" << " "
			<< setw(20)
			<< "Temperature" << " "
			<< setw(20)
			<< "Gamma" << " "
			<< setw(20)
			<< "Coulomb-2*Paul" << " "
			<< setw(20)
			<< "Max_r[um]";

		paulTrapHeaderFile << endl;
		paulTrapHeaderFile.close();
		open();
		close();

		const unsigned int numberOfAtoms = myPositions->size();
		vector<int> n(myTopology->atomTypes.size(), 0);
		for (unsigned int i = 0; i < numberOfAtoms; i++)
		{
			n[myTopology->atoms[i].type]++;
		}

		myPaulQ = 0.0;
		myPaulM = 0.0;
		myPaulOmega = 0.0;
		for (unsigned int i = 0; i < n.size(); i++)
		{
			Real f = myTopology->atomTypes[i].charge / myTopology->atomTypes[0].charge * myTopology->atomTypes[0].mass / myTopology->atomTypes[i].mass;
			Real r = (myOmegaR * myOmegaR + 0.5 * myOmegaZ * myOmegaZ) * f * f - 0.5 * myOmegaZ * myOmegaZ * f;
			Real z = myOmegaZ * myOmegaZ * f;
			myPaulOmega += (Real)n[i] * (2.0 * r + z) / 3.0;
			if (numberOfAtoms > 1)
				myPaulQ += myTopology->atomTypes[i].charge * myTopology->atomTypes[i].charge * n[i] * (n[i] - 1.0) / 2.0;
			else
				myPaulQ += myTopology->atomTypes[i].charge * myTopology->atomTypes[i].charge * n[i];
			myPaulM += myTopology->atomTypes[i].mass * n[i];
			for (unsigned int j = i + 1; j < n.size(); j++)
			{
				myPaulQ += myTopology->atomTypes[i].charge * myTopology->atomTypes[j].charge * n[i] * n[j];
			}
		}

		myPaulOmega = sqrt(myPaulOmega / (Real)numberOfAtoms);

		if (numberOfAtoms > 1)
			myPaulQ = Constant::SQRTCOULOMBCONSTANT * sqrt(myPaulQ / ((Real)numberOfAtoms * ((Real)numberOfAtoms - 1.0) / 2.0));
		else
			myPaulQ = Constant::SQRTCOULOMBCONSTANT * sqrt(myPaulQ);

		if (!myTopology->atomTypes.empty() && myTopology->atomTypes[0].charge < 0)
			myPaulQ = -myPaulQ;

		myPaulK = myPaulM / (Real)numberOfAtoms * myPaulOmega * myPaulOmega * 1.0e7 / 4184.0;
		myPaulA = pow(myPaulQ * myPaulQ / myPaulK, 1.0 / 3.0);
		myPaulR = pow((Real)numberOfAtoms * myPaulQ * myPaulQ / myPaulK, 1.0 / 3.0);
		myPaulUHom = 9.0 / 10.0 * pow((Real)numberOfAtoms, 5.0 / 3.0) * myPaulQ * myPaulQ / myPaulA;
		myPaulF = 1.0 / ((Real)numberOfAtoms * myPaulQ * myPaulQ / myPaulA);
		myPaulF2 = myPaulQ * myPaulQ / myPaulA / Constant::BOLTZMANN;
		myPaulLow = myCache->potentialEnergy();

		report.precision(13);
		report << plain
			<< "Paul Trap O : w           = " << myPaulOmega << "[fs-1]\n"
			<< "            : a           = " << myPaulA * 1e-4 << "[um]\n"
			<< "            : f           = " << pow((*myEnergies)[ScalarStructure::COULOMB] / (2.0 * (*myEnergies)[ScalarStructure::OTHER]), 1.0 / 3.0) << "\n"
			<< "            : 1[kcal/mol] = " << myPaulF * pow(2.0, 1.0 / 3.0) << "[E0]\n"
			<< "            : 1[mu]       = " << 1e4 / myPaulA * pow(0.5, 1.0 / 3.0) << "[r0]\n"
			<< "            : 1[AA]       = " << 1.0 / myPaulA * pow(0.5, 1.0 / 3.0) << "[r0]" << endr;

		myPaulLow = myCache->potentialEnergy();
	}

	void OutputPaulTrap::doRunCached(int step)
	{
		if (myScreen)
		{
			report << plain << "Step : ";
			report.setf(std::ios::right);
			report.width(10);
			report << step << ", Time : ";
			report.width(18);
			report.setf(std::ios::showpoint | std::ios::fixed);
			report.precision(2);
			report << myCache->time() * myFactor << " [" << myUnit << "], TE : ";
			report.precision(4);
			report.width(16);
			report << myCache->totalEnergy() << " [kcal/mol]";
			report << ", E : ";
			report.precision(4);
			report.width(16);
			report << myCache->potentialEnergy() * myPaulF << " [eps]";
			report << ", T : ";
			report.precision(4);
			report.width(10);
			report << myCache->temperature() << " [K]" << endr;
			report.reset();
		}


		const unsigned int numberOfAtoms = myPositions->size();
		Real maxR = 0.0;
		for (unsigned int i = 0; i < numberOfAtoms; i++)
		{
			Real r = (*myPositions)[i].norm();
			if (maxR < r)
				maxR = r;
		}

		stringstream out;
		out << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setw(20)
			<< setprecision(2)
			<< setiosflags(std::ios::showpoint | std::ios::fixed)
			<< myCache->time() << " "
			<< resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setiosflags(std::ios::floatfield)
			<< setprecision(14)
			<< setw(20)
			<< myCache->totalEnergy() << " "
			<< setw(20)
			<< myCache->potentialEnergy() << " "
			<< setw(20)
			<< myCache->potentialEnergy() * myPaulF << " "
			<< setw(20)
			<< myCache->kineticEnergy() << " "
			<< setw(20)
			<< (*myEnergies)[ScalarStructure::COULOMB] << " "
			<< setw(20)
			<< (*myEnergies)[ScalarStructure::OTHER] << " "
			<< setw(20)
			<< (myCache->potentialEnergy() - myPaulUHom) * myPaulF << " "
			<< setw(20)
			<< myCache->temperature() << " "
			<< setw(20)
			<< (myCache->temperature() != 0.0 ? myPaulF2 / myCache->temperature() : Constant::REAL_INFINITY) << " "
			<< setw(20)
			<< (*myEnergies)[ScalarStructure::COULOMB] - 2 * (*myEnergies)[ScalarStructure::OTHER] << " "
			<< setw(20)
			<< maxR * 1e-4 << " " << endl;
		myBuffer << out.str();

		if (myPaulLow > myCache->potentialEnergy() && myDoPaulLow && !myFilenameLow.empty())
		{
			//report.precision(15);
			//report << hint << "Lower state found, U="<<myCache->potentialEnergy()<< " ("<<myCache->potentialEnergy()-myPaulLow<<")."<<endr;

			myLowOut = out.str();
			myLowXYZ = *myPositions;
			myPaulLow = myCache->potentialEnergy();
			myLowComment = string("Time : " + toString(myCache->time()) + ", step : " + toString(step) + ".");
		}
	}

	void OutputPaulTrap::doFlushCache()
	{
		if (!myLowOut.empty())
		{
			XYZWriter pos(myFilenameLow);
			if (!pos.open(myFilenameLow))
				report << error << "Can't open Paul Trap low XYZ file \'" << myFilenameLow << "\'." << endr;
			pos.setComment(myLowComment);
			if (!pos.write(myLowXYZ, myTopology->atoms, myTopology->atomTypes))
				report << error << "Could not write Paul Trap low XYZ file \'" << myFilenameLow << "\'." << endr;

			ofstream paulFile(string(myFilename + ".low").c_str(), std::ios::out);
			if (!paulFile.good())
				report << error << "Can't open Paul Trap low file \'" << myFilename << ".low\'." << endr;

			paulFile << myLowOut;
			paulFile.close();
		}


		myLowOut = "";
	}

	Output* OutputPaulTrap::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputPaulTrap(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]));
	}

	void OutputPaulTrap::getParameters(vector<Parameter>& parameter) const
	{
		OutputFile::getParameters(parameter);
		parameter.push_back(Parameter("paulOmega", Value(myOmegaR, ConstraintValueType::NotNegative())));
		parameter.push_back(Parameter("paulOmegaZ", Value(myOmegaZ, ConstraintValueType::NotNegative())));
		parameter.push_back(Parameter("paulLowFile", Value(myFilenameLow), string("")));
		parameter.push_back(Parameter("doPaulLowFile", (myFilenameLow.empty() ? Value(true, Value::undefined) : Value(myDoPaulLow))));
		parameter.push_back(Parameter("ScreenPaulTrap", Value(myScreen)));
	}

	bool OutputPaulTrap::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
	{
		if (!checkParameterTypes(values))
			return false;

		// Default output frequenxy
		if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
			values[1] = (*config)[InputOutputfreq::keyword];

		// Retrieve omega from Paul Trap force definition
		if (!values[5].valid())
		{
			string str;
			stringstream ss(config->get(InputIntegrator::keyword).getString());
			while (ss >> str)
			{
				Real w;
				if (equalNocase(str, "-omegaR") && ss >> w)
				{
					values[5] = w;
					break;
				}
			}
		}

		// Retrieve omega from Paul Trap force definition
		if (!values[6].valid())
		{
			string str;
			stringstream ss(config->get(InputIntegrator::keyword).getString());
			while (ss >> str)
			{
				Real w;
				if (equalNocase(str, "-omegaZ") && ss >> w)
				{
					values[6] = w;
					break;
				}
			}
		}

		if (values[5].valid() && !values[6].defined())
			values[6] = values[5];
		values[8] = (!values[7].getString().empty() && (!values[7].valid() || (bool)values[8] == true));
		if (!values[9].valid())
			values[9] = config->empty(OutputScreen::keyword);
		return checkParameters(values);
	}
}
