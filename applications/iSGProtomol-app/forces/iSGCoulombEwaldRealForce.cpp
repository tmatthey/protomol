#include "iSGCoulombEwaldRealForce.h"
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {
  //_________________________________________________________________ iSGCoulombEwaldRealForce

  const string iSGCoulombEwaldRealForce::keyword("iSGCoulombEwaldReal");

  iSGCoulombEwaldRealForce::iSGCoulombEwaldRealForce():myAlpha(-1.0){}

  iSGCoulombEwaldRealForce::iSGCoulombEwaldRealForce(Real a):myAlpha(a),
								 myAlphaSquared(a*a),
								 my2AlphaPI(2.0*a/sqrt(M_PI)){}

  void iSGCoulombEwaldRealForce::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-alpha",Value(myAlpha,ConstraintValueType::Positive()),Text("Ewald splitting")));
  }

  iSGCoulombEwaldRealForce iSGCoulombEwaldRealForce::make(std::string& , const std::vector<Value>& values) {
    return iSGCoulombEwaldRealForce(values[0]);
  }

}
