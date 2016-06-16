#include "ForceFactoryDetails.h"
#include "stringutilities.h"
#include "Report.h"
#include "simpleTypes.h"
#include "Parallel.h"
#include "CompareForce.h"
#include "TimeForce.h"
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
	//_____________________________________________________ ForceFactoryDetails
	ForceFactoryDetails::ForceFactoryDetails(): FactoryBase<Force>(), myLastCompareForce(NULL)
	{
		//report << plain <<"[ForceFactoryDetails::ForceFactoryDetails()]"<< endr;
	}

	ForceFactoryDetails::~ForceFactoryDetails()
	{
		//report << plain <<"[ForceFactoryDetails::~ForceFactoryDetails()]"<< endr;
	}

	ForceFactoryDetails::ForceFactoryDetails(const ForceFactoryDetails&)
	{
	}

	ForceFactoryDetails& ForceFactoryDetails::operator=(const ForceFactoryDetails&)
	{
		return *this;
	}

	Force* ForceFactoryDetails::make(string& errMsg, const string& id, vector<Value> values)
	{
		errMsg = "";
		return TFactory::instance().doMake(errMsg, id, values);
	}

	Force* ForceFactoryDetails::doMake(string& errMsg, const string& idInput, vector<Value> values) const
	{
		// Just try .. and see if we can find it directly
		string id = normalizeString(idInput);
		const Force* prototype = getPrototype(id);
		if (prototype != NULL && values.size() != prototype->getParameterSize())
			prototype = NULL;
		bool excat = true;

		if (prototype == NULL)
		{
			// Lookup tables
			if (!myCache)
				updateCache();


			vector<string> splitId(splitForceString(id));

			// Find the force keyword(s) ...
			string keyword = splitId[0];
			map<string, ForceType, ltstrNocase>::const_iterator itr = myForceTypes.find(keyword);

			if (itr == myForceTypes.end())
			{
				// ... sort force keyword(s) and try again ...
				excat = false;
				vector<string> tmp(splitString(keyword));
				sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
				keyword = mergeString(tmp);
				if (myForceTypesSorted.find(keyword) != myForceTypesSorted.end())
				{
					keyword = myForceTypesSorted[keyword];
					itr = myForceTypes.find(keyword);
				}
			}

			string newId;
			if (itr != myForceTypes.end())
			{
				// ... retrieve all valid policies for the force keyword(s) from the string ...
				newId = keyword;
				for (unsigned int i = 1; i < splitId.size(); ++i)
				{
					if (myForceTypes[keyword].policy.find(splitId[i]) != myForceTypes[keyword].policy.end())
						newId += " " + (*myForceTypes[keyword].policy.find(splitId[i]));
				}
				excat = excat && equalBeginNocase(newId, id);
				prototype = getPrototype(newId);
			}

			if (prototype == NULL && itr != myForceTypes.end())
			{
				// ... sort policies and try again ...
				newId = sortForceString(newId);
				for (map<string, string, ltstrNocase>::const_iterator j = itr->second.policiesSorted.begin(); j != itr->second.policiesSorted.end(); ++j)
				{
					if (equalNocase(j->first, newId))
						newId = j->second;
				}
				excat = false;
				prototype = getPrototype(newId);
			}

			if (prototype == NULL && itr != myForceTypes.end())
			{
				// ... remove double policies and try again ...
				newId = uniqueForceString(newId);
				for (map<string, string, ltstrNocase>::const_iterator j = itr->second.policiesSorted.begin(); j != itr->second.policiesSorted.end(); ++j)
				{
					if (equalNocase(uniqueForceString(j->first), newId))
						newId = j->second;
				}
				prototype = getPrototype(newId);
				excat = false;
			}

			//  ... try with alias!
			if (prototype == NULL)
			{
				prototype = getPrototype(splitId[0]);
				if (prototype != NULL)
				{
					newId = prototype->getId();
				}
			}

			// ... ok final last call SK945 ... time and compare
			if (prototype == NULL && equalStartNocase(CompareForce::keyword, id))
			{
				Force* actualForce = doMake(errMsg, id.substr(CompareForce::keyword.size() + 1), values);
				if (actualForce == NULL)
					return NULL;
				CompareForce* compareForce = NULL;
				if (myLastCompareForce == NULL)
				{
					compareForce = actualForce->makeCompareForce(actualForce,NULL);
					myLastCompareForce = compareForce;
				}
				else
				{
					compareForce = actualForce->makeCompareForce(actualForce, myLastCompareForce);
					report << hint << "Comparing " << compareForce->getIdNumber() / 2 << " :\'" << myLastCompareForce->getForceObject()->getId()
						<< "\' with \'" << actualForce->getId() << "\'" << endr;
					myLastCompareForce = NULL;
					if (Parallel::isParallel())
						report << hint << "Comparing forces in parallel environment not supported." << endr;
				}
				return compareForce;
			}
			else if (prototype == NULL && equalStartNocase(TimeForce::keyword, idInput))
			{
				Force* actualForce = doMake(errMsg, id.substr(TimeForce::keyword.size() + 1), values);
				if (actualForce == NULL)
					return NULL;
				TimeForce* timeForce = actualForce->makeTimeForce(actualForce);
				report << hint << "Timing " << timeForce->getIdNumber() << " :\'" << actualForce->getId() << "\'" << endr;
				return timeForce;
			}

			string parametersString;
			if (prototype != NULL)
			{
				// ok, we found a match, retrieve the parameters to
				// be passed the make method ... now remove
				// all force keyword(s) and policies.
				vector<string> splitNewId(splitForceString(newId));
				set<string, ltstrNocase> unique;
				for (unsigned int i = 0; i < splitNewId.size(); ++i)
					unique.insert(splitNewId[i]);

				for (unsigned int i = 0; i < splitId.size(); ++i)
				{
					if (unique.find(splitId[i]) == unique.end())
						parametersString += (parametersString.empty() ? "" : " ") + splitId[i];
				}
				parametersString += " "; // some compiler need this

				vector<Parameter> parameters;
				prototype->getParameters(parameters);
				values.resize(parameters.size());

				stringstream ssp(parametersString);
				bool foundLast = false;

				for (unsigned int i = 0; i < parameters.size(); ++i)
				{
					// parse the the parameters one by one
					string strp;
					values[i] = parameters[i].value;
					values[i].clear();

					//report << i << ":"<<parameters[i].keyword<<","<<(long)ssp.tellg()<<",";
					bool found = false;
					bool retry = true;
					if (!(ssp))
					{
						//report <<"R0,";
						ssp.str(parametersString);
						ssp.clear();
						ssp.seekg(std::ios::beg);
					}
					while (ssp >> strp || retry)
					{
						//report << strp <<",";
						if (!(ssp) && retry)
						{
							//report <<"R1,";
							ssp.str(parametersString);
							ssp.clear();
							ssp.seekg(std::ios::beg);
							retry = false;
							strp = "";
							continue;
						}
						if (equalNocase(strp, parameters[i].keyword) && !parameters[i].keyword.empty())
						{
							//report <<"A,";
							ssp >> values[i];
							found = true;
							break;
						}
						else if (foundLast && parameters[i].keyword.empty())
						{
							//report <<"B,"<<strp<<",";
							ssp.seekg((-1) * static_cast<int>(strp.size()), std::ios::cur);
							ssp.clear();
							ssp >> values[i];
							found = true;
							break;
						}
					}
					foundLast = found;

					// If still undefined take default value is available
					if (!found && parameters[i].defaultValue.valid())
					{
						values[i].set(parameters[i].defaultValue);
					}
					//report << values[i].valid()<<","<<values[i].getString()<<endr;
				}

				id = newId;
			}
			else
			{
				if (itr != myForceTypes.end())
				{
					errMsg += " Could not find a complete match for \'" + idInput + "\' in " + Force::scope + "Factory of force \'" + keyword + "\'. Possible definitions are:";
					for (set<string, ltstrNocase>::const_iterator j = itr->second.policies.begin(); j != itr->second.policies.end(); ++j)
					{
						errMsg += string("\n") + itr->first + string(" ") + (*j);
					}
				}
				else
				{
					errMsg += " Could not find any match for \'" + idInput + "\' in " + Force::scope + "Factory. Possible forces are:";
					for (map<string, ForceType, ltstrNocase>::const_iterator i = myForceTypes.begin(); i != myForceTypes.end(); ++i)
					{
						errMsg += string("\n") + i->first;
					}
				}
			}
		}

		// Make
		Force* newObj = NULL;
		if (prototype != NULL)
		{
			// Check parameter list
			if (!prototype->checkParameters(errMsg, values))
			{
				return NULL;
			}
			newObj = prototype->make(errMsg, values);
			if (!errMsg.empty())
			{
				delete newObj;
				return NULL;
			}
			if (newObj != NULL)
			{
				newObj->setAlias(id);
			}
		}

		return newObj;
	}

	string ForceFactoryDetails::doPrint() const
	{
		if (!myCache)
			updateCache();
		string res;
		for (map<string, ForceType, ltstrNocase>::const_iterator i = myForceTypes.begin(); i != myForceTypes.end(); ++i)
		{
			res += (res.empty() ? "" : "\n") + i->first;
			for (set<string, ltstrNocase>::const_iterator j = i->second.policies.begin(); j != i->second.policies.end(); ++j)
			{
				if (!(*j).empty())
					res += "\n" + Constant::PRINTINDENT + (*j);
				vector<Parameter> parameter = getPrototype(i->first + ((*j).empty() ? "" : " " + (*j)))->getParameters();
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
		}
		res += "\nAlias:";
		for (map<string, const Force*, ltstrNocase>::const_iterator j = myAliasExemplars.begin(); j != myAliasExemplars.end(); ++j)
			res += "\n" + j->first + " : " + j->second->getId() + " (" + j->second->getIdNoAlias() + ")";
		return res;
	}


	void ForceFactoryDetails::updateCache() const
	{
		myForceTypes.clear();
		myForceTypesSorted.clear();

		set<PairString> removes;

		for (set<const Force*>::const_iterator i = myPointers.begin(); i != myPointers.end(); ++i)
		{
			string id = (*i)->getId();
			string idSorted = sortForceString(id);
			vector<string> idSplit = splitForceString(id);
			string keyword = idSplit[0];

			vector<string> tmp(splitString(keyword));
			sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
			string keywordSorted = mergeString(tmp);

			for (unsigned int j = 1; j < idSplit.size(); ++j)
				myForceTypes[keyword].policy.insert(idSplit[j]);

			if (!equalNocase(keywordSorted, keyword))
				myForceTypesSorted[keywordSorted] = keyword;

			tmp.resize(idSplit.size() - 1);
			copy(idSplit.begin() + 1, idSplit.end(), tmp.begin());
			myForceTypes[keyword].policies.insert(mergeString(tmp));

			if (myForceTypes[keyword].policiesSorted.find(idSorted) != myForceTypes[keyword].policiesSorted.end() &&
				!equalNocase(myForceTypes[keyword].policiesSorted[idSorted], id))
			{
				//report << hint <<"[ForceFactoryDetails::updateCache] Sorted force id matching."<<endr;
				removes.insert(PairString(keyword, idSorted));
			}
			else
			{
				myForceTypes[keyword].policiesSorted[idSorted] = id;
			}
		}

		// Remove all sorted id entries which match when sorted
		for (set<PairString>::const_iterator i = removes.begin(); i != removes.end(); ++i)
		{
			myForceTypes[i->first].policiesSorted.erase(myForceTypes[i->first].policiesSorted.find(i->second));
		}

		myCache = true;
	}

	vector<string> ForceFactoryDetails::splitForceString(const string& id) const
	{
		vector<string> res(1);
		vector<string> tmp(splitString(id));
		unsigned int i = 0;
		while (i < tmp.size() && tmp[i][0] != '-')
			res[0] += (res[0].empty() ? "" : " ") + tmp[i++];
		for (; i < tmp.size(); ++i)
		{
			if (tmp[i][0] == '-' && tmp[i].size() > 1 && !isdigit(tmp[i][1]))
				res.push_back(tmp[i]);
			else
				res[res.size() - 1] += " " + tmp[i];
		}
		return res;
	}

	vector<string> ForceFactoryDetails::splitForceStringSorted(const string& id) const
	{
		vector<string> res(splitForceString(id));

		vector<string> tmp(splitString(res[0]));
		sort(tmp.begin(), tmp.end(), ltstrNocaseOp);
		res[0] = mergeString(tmp);

		sort(res.begin() + 1, res.end(), ltstrNocaseOp);
		return res;
	}

	string ForceFactoryDetails::sortForceString(const string& id) const
	{
		return mergeString(splitForceStringSorted(id));
	}

	string ForceFactoryDetails::uniqueForceString(const string& id) const
	{
		vector<string> tmp(splitForceString(id));
		string res = tmp[0];
		set<string, ltstrNocase> unique;
		for (unsigned int i = 1; i < tmp.size(); ++i)
		{
			if (unique.find(tmp[i]) == unique.end())
				res += " " + tmp[i];
			unique.insert(tmp[i]);
		}
		return res;
	}

	void ForceFactoryDetails::doRegisterHelpText() const
	{
	}
}
