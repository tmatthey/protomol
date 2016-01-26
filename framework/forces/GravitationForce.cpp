#include "GravitationForce.h"
using std::string;
namespace ProtoMol {
  //_________________________________________________________________ GravitationForce
  const string GravitationForce::keyword("Gravitation");

  void GravitationForce::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-G",Value(myG)));
  }

  GravitationForce GravitationForce::make(std::string& ,const std::vector<Value>& values) {
    return GravitationForce(values[0]);
  }
}
