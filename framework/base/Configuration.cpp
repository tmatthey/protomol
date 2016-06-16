#include "Configuration.h"
#include "stringutilities.h"
#include "Report.h"

using std::vector;
using std::string;
using std::map;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//________________________________________________________ Configuration
	void Configuration::registerKeyword(const string& keyword, Value value)
	{
		if (end() != myValues.find(keyword) && myAliases.end() == myAliases.find(keyword))
		{
			if (!value.equalTypeAndConstraint(myValues[keyword]))
				report << hint << "[Configuration::registerKeyword] overwriting already registered keyword \'"
					<< keyword << "\' with diffent type/constraint." << endr;
			myValues[keyword] = value;
		}
		else if (myAliases.end() == myAliases.find(keyword))
		{
			myValues[keyword] = value;
		}
		else
		{
			report << hint << "[Configuration::registerKeyword] keyword \'" << keyword
				<< "\' already used for alias \'" << myAliases[keyword] << "\'." << endr;
		}
	}

	void Configuration::registerAliases(const string& keyword, const vector<string>& aliases)
	{
		for (unsigned int i = 0; i < aliases.size(); ++i)
		{
			if (!aliases[i].empty() && !keyword.empty())
			{
				if (myAliases.end() == myAliases.find(aliases[i]))
				{
					myAliases[aliases[i]] = keyword;
				}
				else
				{
					if (!equalNocase(myAliases[aliases[i]], keyword))
						report << hint << "[Configuration::registerAliases] alias \'"
							<< aliases[i] << "\' already used for \'" << myAliases[aliases[i]] << "\'." << endr;
				}
			}
		}
	}

	void Configuration::unregisterKeyword(const std::string& keyword)
	{
		iterator i = myValues.find(keyword);
		if (i != end())
			myValues.erase(i);
		for (AliasMapType::iterator i = myAliases.begin(); i != myAliases.end();)
		{
			if (equalNocase(i->second, keyword))
				myAliases.erase(i++);
			else
				++i;
		}

		TextMapType::iterator j = myTexts.find(keyword);
		if (j != myTexts.end())
			myTexts.erase(j);
	}

	bool Configuration::empty(const string& keyword) const
	{
		if (keyword == "")
			return myValues.empty();

		return (find(keyword) == end());
	}

	bool Configuration::defined(const string& keyword) const
	{
		const_iterator i = find(keyword);

		if (i != end())
			return i->second.defined();

		return false;
	}

	bool Configuration::valid(const string& keyword) const
	{
		const_iterator i = find(keyword);

		if (i != end())
			return i->second.valid();

		return false;
	}

	bool Configuration::set(const string& keyword, Value v)
	{
		iterator i = find(keyword);

		if (i != end())
		{
			i->second.set(v);
			return i->second.valid();
		}
		return false;
	}

	bool Configuration::set(const string& keyword, const string& v)
	{
		iterator i = find(keyword);

		if (i != end())
		{
			if (v == "")
			{
				i->second.init();
			}
			else
			{
				i->second.set(v);
			}
			return i->second.valid();
		}
		return false;
	}


	bool Configuration::set(const vector<vector<string>>& values)
	{
		const unsigned int count = values.size();
		bool ok = true;
		for (unsigned int i = 0; i < count; ++i)
		{
			if (values[i].empty())
				continue;
			iterator itr = find(values[i][0]);
			if (itr == end())
			{
				report << recoverable << "Keyword \'" << values[i][0] << "\' not recognized." << endr;
				continue;
			}
			string val("");
			for (unsigned int j = 1; j < values[i].size(); ++j)
			{
				val += (j == 1 ? "" : " ") + values[i][j];
			}
			if (val == "")
			{
				ok &= itr->second.init();
			}
			else if (itr->second.set(val))
			{
				ok = false;
			}
			if (!itr->second.valid())
				report << hint << "Value \'" << val << "\'  could not be assigned \'"
					<< values[i][0] << "\' " << itr->second.getDefinitionTypeString() << "." << endr;
		}
		return false;
	}

	bool Configuration::set(const string& keyword, const vector<vector<string>>& values)
	{
		const unsigned int count = values.size();
		for (unsigned int i = 0; i < count; ++i)
		{
			if (values[i].empty())
				continue;
			if (!equalNocase(values[i][0], keyword))
				continue;
			string val("");
			for (unsigned int j = 1; j < values[i].size(); ++j)
			{
				val += (j == 1 ? "" : " ") + values[i][j];
			}
			return set(keyword, val);
		}
		return false;
	}

	const Value& Configuration::operator[](const string& keyword) const
	{
		const_iterator i = find(keyword);

		if (i == end())
		{
			report << error << "Configuration index \'" << keyword << "\' of of range." << endr;
		}
		return i->second;
	}

	Value& Configuration::operator[](const string& keyword)
	{
		iterator i = find(keyword);

		if (i == end())
		{
			report << error << "Configuration index \'" << keyword << "\' of of range." << endr;
		}
		return i->second;
	}


	Value Configuration::get(const string& keyword) const
	{
		const_iterator i = find(keyword);

		if (i != end())
		{
			return i->second;
		}
		return Value();
	}


	vector<Value> Configuration::get(const vector<Parameter>& parameters) const
	{
		vector<Value> values(parameters.size());

		for (unsigned int i = 0; i < parameters.size(); ++i)
		{
			values[i] = get(parameters[i].keyword);
		}

		return values;
	}


	bool Configuration::validConfiguration() const
	{
		string tmp;
		return validConfiguration(tmp);
	}

	bool Configuration::validConfiguration(string& errMsg) const
	{
		errMsg = "";
		for (const_iterator i = begin(); i != end(); ++i)
		{
			if (!i->second.valid())
			{
				if (errMsg.empty())
					errMsg = "Undefined keyword(s):";
				errMsg += "\n" + getRightFill(i->first, Constant::PRINTMAXWIDTH) + i->second.getDefinitionTypeString();
			}
		}
		return (errMsg.empty());
	}

	string Configuration::print() const
	{
		string str;
		str += "Keywords:\n\n";
		for (const_iterator i = begin(); i != end(); ++i)
		{
			str += getRightFill(i->first, Constant::PRINTMAXWIDTH) + i->second.getDefinitionTypeString();
			if (myTexts.find(i->first) != myTexts.end())
				str += " \t # " + myTexts.find(i->first)->second;
			str += "\n";
		}
		str += "\nAliases:\n\n";
		for (AliasMapType::const_iterator i = myAliases.begin(); i != myAliases.end(); ++i)
			str += getRightFill(i->first, Constant::PRINTMAXWIDTH) + i->second + "\n";
		return str;
	}

	Configuration::iterator Configuration::find(const string& keyword)
	{
		iterator i = myValues.find(keyword);

		if (i == end())
		{
			AliasMapType::iterator j = myAliases.find(keyword);
			if (j != myAliases.end())
			{
				i = myValues.find(j->second);
			}
		}
		return i;
	}

	Configuration::const_iterator Configuration::find(const string& keyword) const
	{
		const_iterator i = myValues.find(keyword);

		if (i == end())
		{
			AliasMapType::const_iterator j = myAliases.find(keyword);
			if (j != myAliases.end())
			{
				i = myValues.find(j->second);
			}
		}
		return i;
	}

	bool Configuration::setText(const std::string& keyword, const string& text)
	{
		if (!empty(keyword))
		{
			myTexts[keyword] = text;
			return true;
		}
		return false;
	}

	string Configuration::getText(const std::string& keyword) const
	{
		TextMapType::const_iterator i = myTexts.find(keyword);
		if (i != myTexts.end())
		{
			return i->second;
		}
		else
		{
			return string();
		}
	}

	vector<string> Configuration::getAliases(const string& keyword) const
	{
		vector<string> res;
		for (AliasMapType::const_iterator i = myAliases.begin(); i != myAliases.end(); ++i)
		{
			if (equalNocase(i->second, keyword))
			{
				res.push_back(i->first);
			}
		}
		return res;
	}
}
