#include "GenericTopology.h"
using std::string;
using std::vector;
namespace ProtoMol {
  //_________________________________________________________________ GenericTopology

  const string  GenericTopology::scope("Topology");
  const string  GenericTopology::keyword("Topology");

  GenericTopology::GenericTopology():Makeable(),
				     exclude(ExclusionType::ONE4MODIFIED), 
				     coulombScalingFactor(1.0), 
				     time(0.0), 
				     min(Vector3D(Constant::MAXREAL,Constant::MAXREAL,Constant::MAXREAL)),
				     max(Vector3D(-Constant::MINREAL,-Constant::MINREAL,-Constant::MINREAL)),
				     minimalMolecularDistances(false),
  doSCPISM(0){}

  GenericTopology::GenericTopology(Real c, const ExclusionType& e):Makeable(),
								   exclude(e), 
								   coulombScalingFactor(c), 
								   time(0.0), 
								   min(Vector3D(Constant::MAXREAL,Constant::MAXREAL,Constant::MAXREAL)),
								   max(Vector3D(-Constant::MINREAL,-Constant::MINREAL,-Constant::MINREAL)),
								   minimalMolecularDistances(false),
  doSCPISM(0){}



  
  GenericTopology* GenericTopology::make(string& errMsg, const vector<Value>& values)const{
    errMsg = "";
    if(!checkParameters(errMsg,values))
      return NULL;
    return adjustAlias(doMake(errMsg,values));
  }
}
