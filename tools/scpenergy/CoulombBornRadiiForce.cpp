#include "CoulombBornRadiiForce.h"
using std::string;
namespace ProtoMol {
  //_________________________________________________________________ CoulombForce

  const string CoulombBornRadiiForce::keyword("CoulombBornRadii");
  const string CoulombBornRadiiForce::C1::keyword("C1");
  const string CoulombBornRadiiForce::C2::keyword("C2");
  const string CoulombBornRadiiForce::C3::keyword("C3");
  const string CoulombBornRadiiForce::C4::keyword("C4");

    // Default constructor
    CoulombBornRadiiForce::CoulombBornRadiiForce(): sw(1){}
    // Constructor with parameters
    CoulombBornRadiiForce::CoulombBornRadiiForce(int s) : sw(s){}


    void CoulombBornRadiiForce::getParameters(std::vector<Parameter>& parameters) const{
      parameters.push_back(Parameter("-bornswitch",Value(sw,ConstraintValueType::NotNegative()),1,Text("Born switch")));
    }

    CoulombBornRadiiForce CoulombBornRadiiForce::make(std::string&, const std::vector<Value>& values) {
      return CoulombBornRadiiForce(values[0]);
    }

}
