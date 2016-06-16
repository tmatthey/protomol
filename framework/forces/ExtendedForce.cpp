#include "ExtendedForce.h"
#include "ForceGroup.h"
#include "ExtendedCompareForce.h"
#include "ExtendedTimeForce.h"

namespace ProtoMol
{
	//_________________________________________________________________ ExtendedForce

	void ExtendedForce::addToForceGroup(ForceGroup* forceGroup)
	{
		forceGroup->addExtendedForce(this);
	}

	CompareForce* ExtendedForce::makeCompareForce(Force* actualForce, CompareForce* compareForce) const
	{
		return (new ExtendedCompareForce(actualForce, compareForce));
	}

	TimeForce* ExtendedForce::makeTimeForce(Force* actualForce) const
	{
		return (new ExtendedTimeForce(actualForce));
	}
}
