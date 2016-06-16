#include "CoulombEwaldRealForce.h"
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ CoulombEwaldRealForce

	const string CoulombEwaldRealForce::keyword("CoulombEwaldReal");

	CoulombEwaldRealForce::CoulombEwaldRealForce(): myAlpha(-1.0)
	{
	}

	CoulombEwaldRealForce::CoulombEwaldRealForce(Real a): myAlpha(a),
	                                                      myAlphaSquared(a * a),
	                                                      my2AlphaPI(2.0 * a / sqrt(Constant::M_PI))
	{
	}

	void CoulombEwaldRealForce::getParameters(std::vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-alpha", Value(myAlpha, ConstraintValueType::Positive()), Text("Ewald splitting")));
	}

	CoulombEwaldRealForce CoulombEwaldRealForce::make(std::string&, const std::vector<Value>& values)
	{
		return CoulombEwaldRealForce(values[0]);
	}
}
