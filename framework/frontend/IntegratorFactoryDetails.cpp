#include "IntegratorFactoryDetails.h"
#include "STSIntegrator.h"
#include "MTSIntegrator.h"
#include "NonStandardIntegrator.h"
#include "stringutilities.h"
#include "ForceFactory.h"
#include "ForceGroup.h"
#include "Report.h"
#include "HelpTextFactory.h"

using std::vector;
using std::map;
using std::string;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_____________________________________________________ IntegratorFactoryDetails
	IntegratorFactoryDetails::IntegratorFactoryDetails(): FactoryBase<Integrator>()
	{
		//report << plain <<"[IntegratorFactoryDetails::IntegratorFactoryDetails()]"<< endr;
	}

	IntegratorFactoryDetails::~IntegratorFactoryDetails()
	{
		//report << plain <<"[IntegratorFactoryDetails::~IntegratorFactoryDetails()]"<< endr;
	}

	IntegratorFactoryDetails::IntegratorFactoryDetails(const IntegratorFactoryDetails&)
	{
	}

	IntegratorFactoryDetails& IntegratorFactoryDetails::operator=(const IntegratorFactoryDetails&)
	{
		return *this;
	}

	Integrator* IntegratorFactoryDetails::make(string& errMsg, const string& definition)
	{
		return TFactory::instance().doMake(errMsg, definition);
	}

	Integrator* IntegratorFactoryDetails::doMake(string& errMsg, const string& definition) const
	{
		errMsg = "";
		string str;
		vector<IntegratorInput> integratorInput;
		stringstream ss(definition);
		while (ss >> str)
		{
			// Parse level and integrator type
			string levelStr, integrator, d;
			ss >> levelStr >> integrator >> d;
			if (!(equalNocase(str, "level") && isUInt(levelStr) && !integrator.empty() && d == "{"))
			{
				errMsg += " Integration definition mismatch, expecting \'level <number> <integrator> { ... }.";
				return NULL;
			}

			const Integrator* prototype = getPrototype(integrator);
			if (prototype == NULL)
			{
				errMsg += " Could not find any match for \'" + integrator + "\' in " + Integrator::scope + "Factory. Possible integrators are:\n" + doPrint();
				return NULL;
			}

			// Read first integrator parameters and then force definitions
			string parameterStr, forceStr;
			while (ss >> str)
			{
				if (str == "}")
					break;
				if (equalNocase(str, "Force") || !forceStr.empty())
					forceStr += (forceStr.empty() ? "" : " ") + str;
				else
					parameterStr += (parameterStr.empty() ? "" : " ") + str;
			}
			parameterStr += " "; // some compiler need this

			// Expand vector
			unsigned int level = toUInt(levelStr);
			if (integratorInput.size() <= level)
				integratorInput.resize(level + 1);

			// Return if already defined
			if (integratorInput[level].prototype != NULL)
			{
				errMsg += " Level " + toString(level) + " already defined with " + integratorInput[level].prototype->getId() + ".";
				return NULL;
			}

			// Parse integrator parameter
			vector<Parameter> parameters;
			prototype->getParameters(parameters);
			integratorInput[level].values.resize(parameters.size());

			stringstream ssp(parameterStr);
			bool foundLast = false;

			for (unsigned int i = 0; i < parameters.size(); ++i)
			{
				// parse the the parameters one by one
				string strp;
				integratorInput[level].values[i] = parameters[i].value;
				integratorInput[level].values[i].clear();

				//report << debug<< i << ":"<<parameters[i].keyword<<","<<(long)ssp.tellg()<<","<<integratorInput[level].values[i].debug()<<":";
				bool found = false;
				bool retry = true;
				if (!(ssp))
				{
					//report <<"R0,";
					ssp.str(parameterStr);
					ssp.seekg(std::ios::beg);
					ssp.clear();
				}
				//report <<(long)ssp.tellg()<<",";
				while (ssp >> strp || retry)
				{
					//report << strp <<",";
					if (!(ssp) && retry)
					{
						//report <<"R1,";
						ssp.str(parameterStr);
						ssp.seekg(std::ios::beg);
						ssp.clear();
						retry = false;
						strp = "";
						continue;
					}
					if (equalNocase(strp, parameters[i].keyword) && !parameters[i].keyword.empty())
					{
						//report <<"A,";
						ssp >> integratorInput[level].values[i];
						found = true;
						break;
					}
					else if (foundLast && parameters[i].keyword.empty())
					{
						//report <<"B,"<<strp<<",";
						ssp.seekg((-1) * static_cast<int>(strp.size()), std::ios::cur);
						ssp.clear();
						ssp >> integratorInput[level].values[i];
						found = true;
						break;
					}
				}
				foundLast = found;
				//report << (found?"found,":"");

				// If still undefined take default value is available
				if (!found && parameters[i].defaultValue.valid())
				{
					//report << "default,";
					integratorInput[level].values[i].set(parameters[i].defaultValue);
				}
				//report << integratorInput[level].values[i].valid()<<","<<integratorInput[level].values[i].getString()<<endr;
			}

			if (!prototype->checkParameters(errMsg, integratorInput[level].values))
			{
				errMsg = " Level " + toString(level) + " " + errMsg;
				return NULL;
			}
			// Parse forces
			stringstream ssf(forceStr);
			Value force(ValueType::Force(""));
			while (ssf >> force)
			{
				integratorInput[level].forces.push_back(force.getString());
			}

			// Set prototype
			integratorInput[level].prototype = prototype;
		}

		// Check if we have a definition for each level
		bool ok = true;
		for (unsigned int i = 0; i < integratorInput.size(); ++i)
		{
			if (integratorInput[i].prototype == NULL)
			{
				if (ok)
					errMsg += " Missing integrator definitions of level(s):";
				errMsg += " " + toString(i);
				ok = false;
			}
		}
		if (!ok)
		{
			errMsg += ".";
			return NULL;
		}

		// Check if the chain is ok
		ok = true;
		for (unsigned int i = 0; i < integratorInput.size(); ++i)
		{
			const Integrator* prototype = integratorInput[i].prototype;
			if (dynamic_cast<const StandardIntegrator*>(prototype))
			{
				if (!((i == 0 && dynamic_cast<const STSIntegrator*>(prototype)) ||
					i > 0 && dynamic_cast<const MTSIntegrator*>(prototype)))
				{
					if (i > 0)
						errMsg += " Integrator " + toString(prototype->getId()) + " at level " + toString(i) + " is a STS integrator, expected MTS.";
					else

						errMsg += " Integrator " + toString(prototype->getId()) + " at level " + toString(i) + " is a MTS integrator, expected STS.";
					ok = false;
				}
			}
			else if (dynamic_cast<const NonStandardIntegrator*>(prototype))
			{
				errMsg += " NonStandardIntegrator (level " + toString(i) + " " + toString(prototype->getId()) +
					") are not supported by " + Force::scope + "Factory yet.";
				ok = false;
			}
			else
			{
				report << error << "[IntegratorFactoryDetails::doMake] Found an integrator \'" << prototype->getId() << "\' neither a StandardIntegrator nor NonStandardIntegrator." << endr;
			}
		}
		if (!ok)
			return NULL;

		// Now make the integrator chain ... with all forces 
		StandardIntegrator* integrator = NULL;
		for (unsigned int i = 0; i < integratorInput.size(); ++i)
		{
			ForceGroup* forceGroup = new ForceGroup();
			for (unsigned int j = 0; j < integratorInput[i].forces.size(); ++j)
			{
				Force* force = ForceFactory::make(errMsg, integratorInput[i].forces[j]);
				if (force == NULL)
				{
					delete forceGroup;
					delete integrator;
					return NULL;
				}
				forceGroup->addForce(force);
			}

			StandardIntegrator* newIntegrator = NULL;
			if (i > 0)
			{
				newIntegrator = dynamic_cast<const MTSIntegrator*>(integratorInput[i].prototype)->make(errMsg, integratorInput[i].values, forceGroup, integrator);
			}
			else
			{
				newIntegrator = dynamic_cast<const STSIntegrator*>(integratorInput[i].prototype)->make(errMsg, integratorInput[i].values, forceGroup);
			}
			if (newIntegrator == NULL)
			{
				delete forceGroup;
				delete integrator;
				return NULL;
			}
			integrator = newIntegrator;
		}

		return integrator;
	}

	string IntegratorFactoryDetails::doPrint() const
	{
		string res;
		for (map<string, const Integrator*, ltstrNocase>::const_iterator i = myExemplars.begin(); i != myExemplars.end(); ++i)
		{
			res += (i == myExemplars.begin() ? "" : "\n") + i->first;
			vector<Parameter> parameter(i->second->getParameters());
			for (unsigned int k = 0; k < parameter.size(); k++)
			{
				if (!parameter[k].keyword.empty())
				{
					res += "\n" + Constant::PRINTINDENT + Constant::PRINTINDENT + getRightFill(parameter[k].keyword, Constant::PRINTMAXWIDTH);
				}
				res += " " + (parameter[k].defaultValue.valid() ? parameter[k].defaultValue.getDefinitionTypeString() : parameter[k].value.getDefinitionTypeString());
				if (!parameter[k].text.empty())
					res += "\t # " + parameter[k].text;
			}
		}
		res += "\nAlias:";
		for (map<string, const Integrator*, ltstrNocase>::const_iterator j = myAliasExemplars.begin(); j != myAliasExemplars.end(); ++j)
			res += "\n" + j->first + " : " + j->second->getId() + " (" + j->second->getIdNoAlias() + ")";
		return res;
	}

	void IntegratorFactoryDetails::doRegisterHelpText() const
	{
		for (map<string, const Integrator*, ltstrNocase>::const_iterator i = myExemplars.begin(); i != myExemplars.end(); ++i)
		{
			HelpText helpText;
			i->second->getParameters(helpText.parameters);
			helpText.id = i->second->getIdNoAlias();
			helpText.text = i->second->getText();
			helpText.scope = i->second->getScope();
			HelpTextFactory::registerExemplar(i->second->getId(), helpText);

			HelpText alias;
			alias.text = "alias for \'" + i->second->getId() + "\'";
			alias.scope = i->second->getScope();
			for (map<string, const Integrator*, ltstrNocase>::const_iterator j = myAliasExemplars.begin(); j != myAliasExemplars.end(); ++j)
			{
				if (j->second->getIdNoAlias() == i->second->getIdNoAlias())
				{
					alias.id = j->first;
					HelpTextFactory::registerExemplar(alias.id, alias);
				}
			}
		}
	}
}
