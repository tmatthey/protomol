#include "BSplineType.h"
using std::string;

namespace ProtoMol
{
	//________________________________________________________ BSplineType

	const string BSplineEnum::str[static_cast<int>(LAST) - static_cast<int>(FIRST)] = {
		// Order is essential, must be in relation to Enum
		string("undefined"), // Returned when no enum matches
		string("short"),
		string("long"),
	};
}
