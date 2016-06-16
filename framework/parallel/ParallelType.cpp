#include "ParallelType.h"
using std::string;

namespace ProtoMol
{
	//________________________________________________________ ParallelType

	const string ParallelEnum::str[static_cast<int>(LAST) - static_cast<int>(FIRST)] = {
		// Order is essential, must be in relation to Enum
		string("undefined"), // Returned when no enum matches
		string("static"),
		string("dynamic"),
		string("masterSlave")
	};
}
