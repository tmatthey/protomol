#include "ValueType.h"
using namespace ProtoMol::Report;
using std::string;
using std::istream;

namespace ProtoMol
{
	//__________________________________________________________ValueType
	const string ValueEnum::str[static_cast<int>(LAST) - static_cast<int>(FIRST)] = {
		// Order is essential, must be in relation to Enum ordering
		string("undefined"), // Returned when no enum matches
		string("string"),
		string("int"),
		string("uint"),
		string("real"),
		string("boolean"),
		string("coordinates"),
		string("vector"),
		string("integrator"),
		string("force")
	};
}
