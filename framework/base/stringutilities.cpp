#include "stringutilities.h"
#include "mathutilities.h"

#include <algorithm>
#include <stdlib.h>
#include <errno.h>

#include "Report.h"

using std::find;
using std::vector;
using std::string;
using std::stringstream;
using std::transform;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_____________________________________________________________________ uppercase
	string uppercase(const string& word)
	{
		string tmp(word);
		transform(word.begin(), word.end(), tmp.begin(), toupper);
		return tmp;
	}

	//_____________________________________________________________________ lowercase
	string lowercase(const string& word)
	{
		string tmp(word);
		transform(word.begin(), word.end(), tmp.begin(), tolower);
		return tmp;
	}

	//_____________________________________________________________________ equal
	bool equal(const string& s1, const string& s2)
	{
		return (s1 == s2);
	}

	//_____________________________________________________________________ equalNocase
	bool equalNocase(const string& s1, const string& s2)
	{
		const string::size_type i1 = s1.size();
		if (i1 != s2.size())
			return false;
		for (string::size_type i = 0; i < i1; ++i)
			if (toupper(s1[i]) != toupper(s2[i]))
				return false;
		return true;
	}

	//_____________________________________________________________________ equalBegin
	bool equalBegin(const string& s1, const string& s2)
	{
		const string::size_type i1 = s1.size();
		const string::size_type i2 = s2.size();
		if (i1 < i2)
			return (s2.substr(0, i1) == s1);
		else
			return (s1.substr(0, i2) == s2);
	}

	//_____________________________________________________________________ equalBeginNocase
	bool equalBeginNocase(const string& s1, const string& s2)
	{
		return equalBegin(uppercase(s1), uppercase(s2));
	}


	//_____________________________________________________________________ equalStart
	bool equalStart(const string& s1, const string& s2)
	{
		const string::size_type i1 = s1.size();
		if (i1 <= s2.size())
			return (s2.substr(0, i1) == s1);
		else
			return false;
	}

	//_____________________________________________________________________ equalStartNocase
	bool equalStartNocase(const string& s1, const string& s2)
	{
		return equalStart(uppercase(s1), uppercase(s2));
	}

	//_____________________________________________________________________ equalEnd
	bool equalEnd(const string& s1, const string& s2)
	{
		string::size_type i1 = s1.size();
		string::size_type i2 = s2.size();
		if (i1 < i2)
			return (s2.substr(i2 - i1) == s1);
		else
			return (s1.substr(i1 - i2) == s2);
	}

	//_____________________________________________________________________ equalEndNocase
	bool equalEndNocase(const string& s1, const string& s2)
	{
		return equalEnd(uppercase(s1), uppercase(s2));
	}

	//_____________________________________________________________________ equalTerminate
	bool equalTerminate(const string& s1, const string& s2)
	{
		string::size_type i1 = s1.size();
		string::size_type i2 = s2.size();
		if (i1 <= i2)
			return (s2.substr(i2 - i1) == s1);
		else
			return false;
	}

	//_____________________________________________________________________ equalTerminateNocase
	bool equalTerminateNocase(const string& s1, const string& s2)
	{
		return equalTerminate(uppercase(s1), uppercase(s2));
	}


	//_____________________________________________________________________ toStringGeneric
	template <class T>
	inline string toStringGeneric(T x)
	{
		// http://www.bespecific.com/dialog/becodetalk/archive/980405/0058.html
		stringstream ss;
		ss << x;
		return string(ss.str());
	}

	//_____________________________________________________________________ toString
	string toString(Real x)
	{
		stringstream ss;
		ss.precision(sizeof(Real) > sizeof(float) ? 15 : 9);
		ss << x;
		return string(ss.str());
	}

	//_____________________________________________________________________ toString
	string toString(Real x, unsigned int n, unsigned int m)
	{
		stringstream ss;
		ss.setf(std::ios::showpoint | std::ios::fixed);
		ss.precision(m);
		ss.width(n + m + 1);
		ss << x;
		return string(ss.str());
	}


	//_____________________________________________________________________ toString
	string toString(bool x)
	{
		if (x)
			return "true";
		else
			return "false";
	}

	//_____________________________________________________________________ toString
	string toString(const Vector3D& c)
	{
		return string(toString(c.x) + " " + toString(c.y) + " " + toString(c.z));
	}

	//_____________________________________________________________________ toString
	string toString(const vector<Real>& v)
	{
		string res;
		for (unsigned int i = 0; i < v.size(); ++i)
			res += string(i > 0 ? " " : "") + toString(v[i]);
		return res;
	}

	//_____________________________________________________________________ isReal
	bool isReal(const string& word)
	{
		Real r = 0.0;
		return toReal(word, r);
	}

	//_____________________________________________________________________ toReal
	Real toReal(const string& word)
	{
		Real r = 0.0;
		toReal(word, r);
		return r;
	}

	//_____________________________________________________________________ toReal
	// http://www.dinkumware.com/htm_cpl/stdlib.html#strtod
	bool toReal(const string& word, Real& r)
	{
		char* endptr = NULL;
		double d = strtod(word.c_str(), &endptr);
		r = static_cast<Real>(d);
		return (!word.empty() && ((fabs(d) >= Constant::MINREAL && fabs(d) <= Constant::MAXREAL) || fabs(d) == 0.0) && errno != ERANGE && (endptr == NULL || isBlank(string(endptr))) && isPrintable(word));
	}

	//_____________________________________________________________________ isInt
	bool isInt(const string& word)
	{
		int i = 0;
		return toInt(word, i);
	}

	//_____________________________________________________________________ toInt
	int toInt(const string& word)
	{
		int i = 0;
		toInt(word, i);
		return i;
	}

	//_____________________________________________________________________ toInt
	// http://www.dinkumware.com/htm_cpl/stdlib.html#strtol
	bool toInt(const string& word, int& i)
	{
		char* endptr = NULL;
		long l = strtol(word.c_str(), &endptr, 10);
		i = static_cast<int>(l);
		if (!word.empty() && static_cast<long>(i) == l && errno != ERANGE && (endptr == NULL || isBlank(string(endptr))) && isPrintable(word))
			return true;

		Real r;
		if (toReal(word, r))
		{
			i = static_cast<unsigned int>(r);
			return (static_cast<Real>(i) == r);
		}
		return false;
	}

	//_____________________________________________________________________ isUInt
	bool isUInt(const string& word)
	{
		unsigned int i = 0;
		return toUInt(word, i);
	}

	//_____________________________________________________________________ toUInt
	unsigned int toUInt(const string& word)
	{
		unsigned int i = 0;
		toUInt(word, i);
		return i;
	}

	//_____________________________________________________________________ toUInt
	// http://www.dinkumware.com/htm_cpl/stdlib.html#strtol
	bool toUInt(const string& word, unsigned int& i)
	{
		char* endptr = NULL;
		unsigned long l = strtoul(word.c_str(), &endptr, 10);
		i = static_cast<unsigned int>(l);
		if (!word.empty() && static_cast<unsigned long>(i) == l && errno != ERANGE && (endptr == NULL || isBlank(string(endptr))) && isPrintable(word))
			return true;
		Real r;
		if (toReal(word, r))
		{
			i = static_cast<unsigned int>(r);
			return (static_cast<Real>(i) == r);
		}
		return false;
	}

	//_____________________________________________________________________ isBool
	bool isBool(const string& word)
	{
		bool b = false;
		return toBool(word, b);
	}

	//_____________________________________________________________________ toBool
	bool toBool(const string& word)
	{
		bool b = false;
		toBool(word, b);
		return b;
	}

	//_____________________________________________________________________ toBool
	bool toBool(const string& word, bool& b)
	{
		string s = removeBeginEndBlanks(word);
		if (equalNocase(s, "true") || equalNocase(s, "yes") || equalNocase(s, "on") || equalNocase(s, "1"))
		{
			b = true;
			return true;
		}
		else if (equalNocase(s, "false") || equalNocase(s, "no") || equalNocase(s, "off") || equalNocase(s, "0"))
		{
			b = false;
			return true;
		}
		else
		{
			return false;
		}
	}

	//_____________________________________________________________________ isVector3D
	bool isVector3D(const string& word)
	{
		Vector3D c(0.0, 0.0, 0.0);
		return toVector3D(word, c);
	}

	//_____________________________________________________________________ toVector3D
	Vector3D toVector3D(const string& word)
	{
		Vector3D c(0.0, 0.0, 0.0);
		toVector3D(word, c);
		return c;
	}

	//_____________________________________________________________________ toVector3D
	bool toVector3D(const string& word, Vector3D& c)
	{
		string s = removeBeginEndBlanks(word);
		stringstream ss(s);
		string x, y, z;
		ss >> x >> y >> z;
		bool bx, by, bz;
		bx = toReal(x, c.x);
		by = toReal(y, c.y);
		bz = toReal(z, c.z);
		return (ss.eof() && bx && by && bz);
	}

	//_____________________________________________________________________ isVector
	bool isVector(const string& word)
	{
		vector<Real> v;
		return toVector(word, v);
	}

	//_____________________________________________________________________ toVector
	vector<Real> toVector(const string& word)
	{
		vector<Real> v;
		toVector(word, v);
		return v;
	}

	//_____________________________________________________________________ toVector
	bool toVector(const string& word, vector<Real>& v)
	{
		string s = removeBeginEndBlanks(word);
		stringstream is(s);
		v.clear();
		string str;
		is >> str;
		if (isReal(str))
		{
			v.push_back(toReal(str));
			while (is >> str)
			{
				if (!isReal(str))
					break;
				v.push_back(toReal(str));
			}
		}
		else if (str.size() > 2 && str[0] == '-' && str[1] == '-' && isUInt(str.substr(2)))
		{
			unsigned int n = toUInt(str.substr(2));
			for (unsigned int i = 0; i < n; i++)
			{
				if (!(is >> str))
					return false;
				if (!isReal(str))
				{
					is.seekg((-1) * static_cast<int>(str.size()), std::ios::cur);
					is.clear();
					return false;
				}
				v.push_back(toReal(str));
			}
			return true;
		}
		is.seekg((-1) * static_cast<int>(str.size()), std::ios::cur);
		is.clear();
		return true;
	}

	//_____________________________________________________________________ isBlank
	bool isBlank(const string& word)
	{
		return (word.begin() == std::find_if(word.begin(), word.end(), ProtoMol::isblankchar));
	}

	//_____________________________________________________________________ isblankchar
	bool isblankchar(char c)
	{
		return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
	}

	//_____________________________________________________________________ isPrintable
	bool isPrintable(const string& word)
	{
		return (word.begin() == std::find_if(word.begin(), word.end(), isprintablechar));
	}

	//_____________________________________________________________________ isprintablechar
	bool isprintablechar(char c)
	{
		return (isblankchar(c) || isprint(c));
	}

	//_____________________________________________________________________ getBegin
	string getBegin(const string& s, string::size_type n)
	{
		if (s.size() <= n)
			return s;
		return s.substr(0, n);
	}

	//_____________________________________________________________________ getEnd
	string getEnd(const string& s, string::size_type n)
	{
		const string::size_type i = s.size();
		if (i <= n)
			return s;
		return s.substr(i - n);
	}

	//_____________________________________________________________________ getRightFill
	string getRightFill(const string& s, string::size_type n)
	{
		const string::size_type i = s.size();
		if (i < n)
			return s + string(n - i, ' ');
		return s.substr(0, n);
	}

	//_____________________________________________________________________ getLeftFill
	string getLeftFill(const string& s, string::size_type n)
	{
		const string::size_type i = s.size();
		if (i < n)
			return string(n - i, ' ') + s;
		return s.substr(0, n);
	}

	//_____________________________________________________________________ removeBeginEndBlanks
	string removeBeginEndBlanks(const string& s)
	{
		string::size_type a = s.find_first_not_of(" \t\n\r");
		if (a == string::npos)
			return "";
		return string(&s[a], &s[s.find_last_not_of(" \t\n\r") + 1]);
	}

	//_____________________________________________________________________ ltstrNocase
	bool ltstrNocase::operator()(const string& s1, const string& s2) const
	{
		return strcmp(uppercase(s1).c_str(), uppercase(s2).c_str()) < 0;
	}

	//_____________________________________________________________________ ltstrNocaseOp
	bool ltstrNocaseOp(const string& s1, const string& s2)
	{
		return strcmp(uppercase(s1).c_str(), uppercase(s2).c_str()) < 0;
	}

	//_____________________________________________________________________ equalWildcard
	int equalWildcard(const string& wildcard, const string& name)
	{
		// Match with no wildcards
		if (wildcard == name)
		{
			return 2;
		}

		// Return if no wildcards found
		if (find(wildcard.begin(), wildcard.end(), '*') == wildcard.end() &&
			find(wildcard.begin(), wildcard.end(), '%') == wildcard.end() &&
			find(wildcard.begin(), wildcard.end(), '#') == wildcard.end() &&
			find(wildcard.begin(), wildcard.end(), '+') == wildcard.end())
		{
			return 0;
		}

		// Move to first wildcard
		unsigned int pos = 0;
		for (unsigned int i = 0; i < wildcard.size(); i++)
		{
			pos = i;
			if (wildcard[i] == '*' ||
				wildcard[i] == '%' ||
				wildcard[i] == '#' ||
				wildcard[i] == '+')
				break;
			if (!(i < name.size()))
			{
				return 0;
			}
			if (wildcard[i] != name[i])
			{
				return 0;
			}
		}

		if (pos + 1 == wildcard.size())
		{
			// Test if last wildcard
			if (wildcard[pos] == '*')
			{
				return 1;
			}
			else if (wildcard[pos] == '#')
			{
				for (unsigned int i = pos; i < name.size(); i++)
				{
					if (!isdigit(name[i]))
					{
						return 0;
					}
				}
				return 1;
			}
			else if (wildcard[pos] == '%' && pos + 1 == name.size())
			{
				return 1;
			}
			else if (wildcard[pos] == '+' && pos + 1 == name.size() && isdigit(name[pos]))
			{
				return 1;
			}
		}
		else
		{
			// Recursive test if inside wildcard
			if (wildcard[pos] == '*')
			{
				int ok = 0;
				for (unsigned int i = pos; i <= name.size(); i++)
				{
					if (equalWildcard(string(wildcard.begin() + pos + 1, wildcard.end()), string(name.begin() + i, name.end())) > 0)
						ok = 1;
				}
				return ok;
			}
			else if (wildcard[pos] == '%')
			{
				if (pos < name.size() &&
					equalWildcard(string(wildcard.begin() + pos + 1, wildcard.end()), string(name.begin() + pos + 1, name.end())) > 0)
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
			else if (wildcard[pos] == '#')
			{
				int ok = 0;
				for (unsigned int i = pos; i <= name.size(); i++)
				{
					if (equalWildcard(string(wildcard.begin() + pos + 1, wildcard.end()), string(name.begin() + i, name.end())) > 0)
						ok = 1;
					if (i < name.size() && !isdigit(name[i]))
						break;
				}
				return ok;
			}
			else if (wildcard[pos] == '+')
			{
				if (pos < name.size() && isdigit(name[pos]) &&
					equalWildcard(string(wildcard.begin() + pos + 1, wildcard.end()), string(name.begin() + pos + 1, name.end())) > 0)
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
		}
		return 0;
	}

	//_____________________________________________________________________ splitString
	vector<string> splitString(const string& id)
	{
		stringstream ss(id);
		vector<string> res;
		string str;
		while (ss >> str)
		{
			if (!str.empty())
				res.push_back(str);
		}
		return res;
	}

	//_____________________________________________________________________ mergeString
	string mergeString(const vector<string>& id)
	{
		string res;
		for (unsigned int i = 0; i < id.size(); i++)
		{
			res += (i > 0 ? " " : "") + id[i];
		}
		return res;
	}

	//_____________________________________________________________________ normalizeString
	string normalizeString(const string& word)
	{
		stringstream ss(word);
		string res, str;
		while (ss >> str)
		{
			if (!str.empty())
				res += (res.empty() ? "" : " ") + str;
		}
		return res;
	}

	//_____________________________________________________________________ headString
	string headString(const string& word)
	{
		stringstream ss(word);
		string str;
		ss >> str;
		return str;
	}

	//_____________________________________________________________________ tailString
	string tailString(const string& word)
	{
		stringstream ss(word);
		string res, str;
		ss >> str;
		while (ss >> str)
		{
			if (!str.empty())
				res += (res.empty() ? "" : " ") + str;
		}
		return res;
	}
}
