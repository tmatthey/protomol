#include "OutputFactoryDetails.h"
#include "OutputCollection.h"
#include "Configuration.h"
#include "stringutilities.h"
#include "OutputFile.h"
#include "HelpTextFactory.h"
#include <algorithm>

using std::copy;
using std::sort;
using std::vector;
using std::map;
using std::set;
using std::string;
using std::stringstream;
using namespace ProtoMol::Report;


namespace ProtoMol
{
	//_____________________________________________________ OutputFactoryDetails
	OutputFactoryDetails::OutputFactoryDetails(): FactoryBase<Output>()
	{
		//report << plain <<"[OutputFactoryDetails::OutputFactoryDetails()]"<< endr;
	}

	OutputFactoryDetails::~OutputFactoryDetails()
	{
		//report << plain <<"[OutputFactoryDetails::~OutputFactoryDetails()]"<< endr;
	}

	OutputFactoryDetails::OutputFactoryDetails(const OutputFactoryDetails&)
	{
	}

	OutputFactoryDetails& OutputFactoryDetails::operator=(const OutputFactoryDetails&)
	{
		return *this;
	}

	void OutputFactoryDetails::registerAllExemplarsConfiguration(Configuration* config)
	{
		TFactory::instance().doRegisterAllExemplarsConfiguration(config);
	}

	void OutputFactoryDetails::doRegisterAllExemplarsConfiguration(Configuration* config) const
	{
		for (set<const Output*>::const_iterator itr = myPointers.begin(); itr != myPointers.end(); ++itr)
		{
			const Output* prototype = (*itr);
			vector<Parameter> parameter(prototype->getParameters());
			for (unsigned int i = 0; i < parameter.size(); i++)
			{
				config->registerKeyword(parameter[i].keyword, parameter[i].defaultValue);
				if (!parameter[i].text.empty())
					config->setText(parameter[i].keyword, parameter[i].text);
			}
			if (prototype->addDoKeyword())
				config->registerKeyword("do" + prototype->getId(), Value(true, Value::undefined));
		}
	}


	Output* OutputFactoryDetails::make(string& errMsg, const string& id, const vector<Value>& values)
	{
		return TFactory::instance().doMake(errMsg, id, values);
	}

	Output* OutputFactoryDetails::doMake(string& errMsg, const string& id, const vector<Value>& values) const
	{
		errMsg = "";

		// Find prototype
		const Output* prototype = getPrototype(id);

		if (prototype == NULL)
		{
			errMsg += " Could not find any match for \'" + id + "\' in " + Output::scope + "Factory. Possible outputs are:\n" + doPrint();
			return NULL;
		}

		// Make
		Output* newObj = prototype->make(errMsg, values);
		if (newObj == NULL)
			return NULL;

		// Adjust external alias
		newObj->setAlias(id);
		return newObj;
	}

	OutputCollection* OutputFactoryDetails::makeCollection(string& errMsg, const Configuration* config)
	{
		return TFactory::instance().doMakeCollection(errMsg, config);
	}

	OutputCollection* OutputFactoryDetails::doMakeCollection(string& errMsg, const Configuration* config) const
	{
		errMsg = "";
		OutputCollection* res = new OutputCollection();
		for (Configuration::const_iterator i = config->begin(); i != config->end(); ++i)
		{
			if ((*i).second.valid())
			{
				const Output* prototype = getPrototype((*i).first);
				if (prototype != NULL)
				{
					if (prototype->isIdDefined(config))
					{
						vector<Parameter> parameter(prototype->getParameters());
						vector<Value> values(parameter.size());
						for (unsigned int k = 0; k < parameter.size(); k++)
						{
							if (config->valid(parameter[k].keyword) && !(parameter[k].keyword.empty()))
							{
								values[k].set((*config)[parameter[k].keyword]);
							}
							else if (parameter[k].defaultValue.valid())
							{
								values[k].set(parameter[k].defaultValue);
							}
							else
							{
								values[k] = parameter[k].value;
								values[k].clear();
							}
						}
						if (!prototype->checkParameters(values))
							prototype->adjustWithDefaultParameters(values, config);
						if (prototype->checkParameters(errMsg, values))
						{
							string err;
							res->adoptOutput(doMake(err, (*i).first, values));
							errMsg += err;
						}
					}
				}
			}
		}
		return res;
	}

	string OutputFactoryDetails::doPrint() const
	{
		string res;

		for (map<string, const Output*, ltstrNocase>::const_iterator i = myExemplars.begin(); i != myExemplars.end(); ++i)
		{
			res += (i == myExemplars.begin() ? "" : "\n") + i->first;
			vector<Parameter> parameter(i->second->getParameters());
			for (unsigned int k = 0; k < parameter.size(); k++)
			{
				if (!parameter[k].keyword.empty())
				{
					res += "\n" + Constant::PRINTINDENT + Constant::PRINTINDENT + getRightFill(parameter[k].keyword, Constant::PRINTMAXWIDTH);
				}
				res += (parameter[k].defaultValue.valid() ? parameter[k].defaultValue.getDefinitionTypeString() : parameter[k].value.getDefinitionTypeString());
				if (!parameter[k].text.empty())
					res += "\t # " + parameter[k].text;
			}
		}
		res += "\nAlias:";
		for (map<string, const Output*, ltstrNocase>::const_iterator j = myAliasExemplars.begin(); j != myAliasExemplars.end(); ++j)
			res += "\n" + j->first + " : " + j->second->getId() + " (" + j->second->getIdNoAlias() + ")";
		return res;
	}

	void OutputFactoryDetails::doRegisterHelpText() const
	{
		for (map<string, const Output*, ltstrNocase>::const_iterator i = myExemplars.begin(); i != myExemplars.end(); ++i)
		{
			HelpText helpText;
			i->second->getParameters(helpText.parameters);
			helpText.id = i->second->getIdNoAlias();
			helpText.text = i->second->getText();
			helpText.scope = i->second->getScope();
			if (i->second->addDoKeyword())
			{
				helpText.parameters.push_back(Parameter("do" + i->second->getId(), Value(true, Value::undefined), Text("flag to switch on/off the output")));
			}
			HelpTextFactory::registerExemplar(i->second->getId(), helpText);

			string txt = string("a parameter of " + i->second->getId());
			for (unsigned int j = 0; j < helpText.parameters.size(); j++)
			{
				helpText.text = txt;
				if (!equalNocase(i->second->getId(), helpText.parameters[j].keyword))
				{
					if (!helpText.parameters[j].text.empty())
					{
						helpText.text = txt + ", " + helpText.parameters[j].text;
					}
					HelpTextFactory::registerExemplar(helpText.parameters[j].keyword, helpText);
				}
			}

			HelpText alias;
			alias.text = "alias for \'" + i->second->getId() + "\'";
			alias.scope = i->second->getScope();
			for (map<string, const Output*, ltstrNocase>::const_iterator j = myAliasExemplars.begin(); j != myAliasExemplars.end(); ++j)
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
