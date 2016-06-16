#include "OutputTemperatures.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "inputValueDefinitions.h"
#include "topologyutilities.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::set;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;
using std::ifstream;

namespace ProtoMol
{
	//________________________________________________________ OutputTemperatures
	const string OutputTemperatures::keyword("TempFile");

	OutputTemperatures::OutputTemperatures(): OutputFile(), mySeparateWater(true)
	{
	}

	OutputTemperatures::OutputTemperatures(const string& filename,
	                                       int freq,
	                                       int cacheFreq,
	                                       int cacheSize,
	                                       Real closeTime,
	                                       bool separateWater): OutputFile(filename, freq, cacheFreq, cacheSize, closeTime),
	                                                            mySeparateWater(separateWater)
	{
	}

	OutputTemperatures::~OutputTemperatures()
	{
	}

	void OutputTemperatures::doInitialize()
	{
		// Create a header file
		ofstream temperatureHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
		if (!temperatureHeaderFile)
			report << error << " Can not open \'" << myFilename << ".header\' for " << getId() << "." << endr;

		// Determine the atom types of interest
		int atomNum;
		int atomType;
		if (mySeparateWater == false)
		{
			myNonWaterAtomTypes.resize(myTopology->atomTypes.size());
			for (unsigned int i = 0; i < myTopology->atomTypes.size(); i++)
			{
				myNonWaterAtomTypes[i] = i;
			}
		}
		else
		{
			// Examine each molecule
			bool examinedWater = false; // Only look at water once

			for (unsigned int i = 0; i < myTopology->molecules.size(); i++)
			{
				if (myTopology->molecules[i].water == true)
				{
					if (examinedWater == false)
					{
						examinedWater = true;
						for (unsigned int j = 0; j < myTopology->molecules[i].atoms.size(); j++)
						{
							atomNum = myTopology->molecules[i].atoms[j];
							atomType = myTopology->atoms[atomNum].type;

							bool alreadyAdded = false;
							unsigned int k = 0;
							while ((alreadyAdded == false) && (k < myWaterAtomTypes.size()))
							{
								if (myWaterAtomTypes[k] == atomType)
								{
									alreadyAdded = true;
								}
								k++;
							}

							if (alreadyAdded == false)
							{
								myWaterAtomTypes.push_back(atomType);
							}
						}
					}
				}
				else
				{
					for (unsigned int j = 0; j < myTopology->molecules[i].atoms.size(); j++)
					{
						atomNum = myTopology->molecules[i].atoms[j];
						atomType = myTopology->atoms[atomNum].type;

						bool alreadyAdded = false;
						unsigned int k = 0;
						while ((alreadyAdded == false) && (k < myNonWaterAtomTypes.size()))
						{
							if (myNonWaterAtomTypes[k] == atomType)
							{
								alreadyAdded = true;
							}
							k++;
						}

						if (alreadyAdded == false)
						{
							myNonWaterAtomTypes.push_back(atomType);
						}
					}
				}
			}
		}

		if (myWaterAtomTypes.size() == 0)
		{
			mySeparateWater = false;
		}

		if (mySeparateWater == false)
		{
			temperatureHeaderFile << setw(18) << "Time(fs)" << " ";
			temperatureHeaderFile << setw(24) << "Temp(K)" << " ";

			for (unsigned int i = 0; i < myNonWaterAtomTypes.size(); i++)
			{
				string strHeader = myTopology->atomTypes[myNonWaterAtomTypes[i]].name + " Temp(K)";
				temperatureHeaderFile << setw(24) << strHeader << " ";
			}
		}
		else
		{
			temperatureHeaderFile << setw(18) << "" << " ";
			temperatureHeaderFile << setw(24) << "" << " ";
			temperatureHeaderFile << setw(24) << "" << " ";
			temperatureHeaderFile << setw(24) << "" << " ";
			temperatureHeaderFile << setw(24) << "WATER" << " ";
			for (unsigned int i = 1; i < myWaterAtomTypes.size(); i++)
			{
				temperatureHeaderFile << setw(24) << "" << " ";
			}
			temperatureHeaderFile << setw(24) << "NON-WATER" << " ";
			for (unsigned int i = 1; i < myNonWaterAtomTypes.size(); i++)
			{
				temperatureHeaderFile << setw(24) << "" << " ";
			}
			temperatureHeaderFile << endl;

			temperatureHeaderFile << setw(18) << "Time(fs)" << " ";
			temperatureHeaderFile << setw(24) << "Temp(K)" << " ";
			temperatureHeaderFile << setw(24) << "Water Temp(K)" << " ";
			temperatureHeaderFile << setw(24) << "Non-Water Temp(K)" << " ";
			for (unsigned int i = 0; i < myWaterAtomTypes.size(); i++)
			{
				string strHeader = myTopology->atomTypes[myWaterAtomTypes[i]].name + " Temp(K)";
				temperatureHeaderFile << setw(24) << strHeader << " ";
			}
			for (unsigned int i = 0; i < myNonWaterAtomTypes.size(); i++)
			{
				string strHeader = myTopology->atomTypes[myNonWaterAtomTypes[i]].name + " Temp(K)";
				temperatureHeaderFile << setw(24) << strHeader << " ";
			}
		}

		temperatureHeaderFile.close();

		open();
		close();
	}

	void OutputTemperatures::doRunCached(int)
	{
		// Output the current time
		myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setw(18)
			<< setprecision(2)
			<< setiosflags(std::ios::showpoint | std::ios::fixed)
			<< myCache->time() << " ";

		// Output the current temperature
		myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setw(24)
			<< setprecision(8)
			<< setiosflags(std::ios::showpoint | std::ios::fixed)
			<< myCache->temperature() << " ";

		// If requested, output the separate water
		// and non-water temperatures
		if (mySeparateWater == true)
		{
			myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
				<< setw(24)
				<< setprecision(8)
				<< setiosflags(std::ios::showpoint | std::ios::fixed)
				<< myCache->temperatureForWater() << " ";
			myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
				<< setw(24)
				<< setprecision(8)
				<< setiosflags(std::ios::showpoint | std::ios::fixed)
				<< myCache->temperatureForNonWater() << " ";
			for (unsigned int i = 0; i < myWaterAtomTypes.size(); i++)
			{
				myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
					<< setw(24)
					<< setprecision(8)
					<< setiosflags(std::ios::showpoint | std::ios::fixed)
					<< ProtoMol::temperatureForAtomType(myTopology, myVelocities, myWaterAtomTypes[i], ONLY_WATER) << " ";
			}
			for (unsigned int i = 0; i < myNonWaterAtomTypes.size(); i++)
			{
				myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
					<< setw(24)
					<< setprecision(8)
					<< setiosflags(std::ios::showpoint | std::ios::fixed)
					<< ProtoMol::temperatureForAtomType(myTopology, myVelocities, myNonWaterAtomTypes[i], IGNORE_WATER) << " ";
			}
		}
		else
		{
			for (unsigned int i = 0; i < myNonWaterAtomTypes.size(); i++)
			{
				myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
					<< setw(24)
					<< setprecision(8)
					<< setiosflags(std::ios::showpoint | std::ios::fixed)
					<< ProtoMol::temperatureForAtomType(myTopology, myVelocities, myNonWaterAtomTypes[i], ALL) << " ";
			}
		}

		myBuffer << endl;
	}

	Output* OutputTemperatures::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputTemperatures(values[0], values[1], values[2], values[3], values[4], values[5]));
	}

	void OutputTemperatures::getParameters(vector<Parameter>& parameter) const
	{
		OutputFile::getParameters(parameter);
		parameter.push_back(Parameter("TemperatureSeparateWater", Value(mySeparateWater)));
	}
}
