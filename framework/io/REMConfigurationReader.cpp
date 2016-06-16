#include "REMConfigurationReader.h"

#include "Report.h"
#include "stringutilities.h"

using std::string;
using std::vector;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ConfigurationReader

	REMConfigurationReader::REMConfigurationReader(Real temp): Reader(), myConfig(NULL), myTemp(temp)
	{
	}

	REMConfigurationReader::REMConfigurationReader(const std::string& filename): Reader(filename), myConfig(NULL)
	{
	}

	REMConfigurationReader::~REMConfigurationReader()
	{
		if (myConfig != NULL)
			delete myConfig;
	}

	bool REMConfigurationReader::tryFormat()
	{
		open();
		return !myFile.fail();
	}

	bool REMConfigurationReader::read()
	{
		if (myConfig == NULL)
			myConfig = new Configuration();
		return read(*myConfig);
	}

	bool REMConfigurationReader::read(Configuration& config)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;

		// Remove comments and reformat
		stringstream all;
		int a, b;
		while (!myFile.eof() && !myFile.fail())
		{
			string line(getline());
			//if (line.find("File ", 0) != string::npos) { // begin additions
			//  line = line + string(fcvt(myTemp, 0, &a, &b));
			//} 
			if (line.find("temperature ", 0) != string::npos)
			{
				line = line.substr(0, line.find("temperature ", 0) + 12);
				line = line + string(fcvt(myTemp, 0, &a, &b));
			} // end additions
			stringstream ss(string(line.begin(), std::find(line.begin(), line.end(), '#')));
			string str;
			while (ss >> str)
			{
				all << (all.str().empty() ? "" : " ") << str;
			}
		}

		close();
		if (myFile.fail())
			return false;


		// Nothing to do ...
		if (all.str().empty())
			return true;


		// First get the keyword and then let Value read from istream ...
		string str;
		string bad;
		bool res = true;
		while (all >> str)
		{
			if (!config.empty(str))
			{
				if (!bad.empty())
				{
					report << recoverable << "Ignoring:" << bad << endr;
					bad = "";
				}
				std::ios::pos_type start = all.tellg();
				all >> config[str];
				if (!config[str].valid())
				{
					std::ios::pos_type end = all.tellg();
					if (end > start)
					{
						std::streamsize len = static_cast<std::streamsize>(end - start);
						string tmp(len, ' ');
						all.seekg(start);
						all.read(&(tmp[0]), len);
						report << recoverable << "Could not parse \'" << removeBeginEndBlanks(tmp) << "\' for keyword \'" << str
							<< "\', expecting type " << config[str].getDefinitionTypeString() << "." << endr;
					}
				}
			}
			else
			{
				res = false;
				bad += " " + str;
			}
		}
		if (!bad.empty())
			report << recoverable << "Ignoring:" << bad << endr;
		return res;
	}

	Configuration* REMConfigurationReader::orphanConfiguration()
	{
		Configuration* tmp = myConfig;
		myConfig = NULL;
		return tmp;
	}

	REMConfigurationReader& operator>>(REMConfigurationReader& configReader, Configuration& config)
	{
		configReader.read(config);
		return configReader;
	}
}
