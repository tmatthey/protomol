#include "SystemForce.h"
#include "ForceGroup.h"
#include "SystemCompareForce.h"
#include "SystemTimeForce.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ SystemForce

	void SystemForce::addToForceGroup(ForceGroup* forceGroup)
	{
		forceGroup->addSystemForce(this);
	}

	CompareForce* SystemForce::makeCompareForce(Force* actualForce, CompareForce* compareForce) const
	{
		return new SystemCompareForce(actualForce, compareForce);
	}

	TimeForce* SystemForce::makeTimeForce(Force* actualForce) const
	{
		return new SystemTimeForce(actualForce);
	}
}
