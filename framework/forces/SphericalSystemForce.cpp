#include "SphericalSystemForce.h"
#include "GenericTopology.h"
#include "Vector3DBlock.h"
#include "Parallel.h"
#include "ScalarStructure.h"
#include "mathutilities.h"
#include "pmconstants.h"
using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {
  //_________________________________________________________________ SphericalSystemForce
  const string SphericalSystemForce::keyword("Spherical");

  //static const Real SI_FORCE_FACTOR = 1e-3*Constant::SI::KCAL/Constant::SI::LENGTH_AA;
  //static const Real SI_ENERGY_FACTOR = 1e-3*Constant::SI::KCAL;

  SphericalSystemForce::SphericalSystemForce():SystemForce(),
					       myCenter(Vector3D(0.0,0.0,0.0)),
					       myRadius(0.0),
					       myK(0.0),
					       myJ(0),
					       // myU(0.0),
					       myRadius2(0.0)
  {
  }

  SphericalSystemForce::SphericalSystemForce(Vector3D center, Real radius, Real k, int j /*, Real u*/):SystemForce(),
												       myCenter(center),
												       myRadius(radius),
												       myK(k),
												       myJ(j),
												       // myU(u),
												       myRadius2(radius*radius)
  {
  }

  void SphericalSystemForce::doEvaluate(const GenericTopology* /*topo*/,
					const Vector3DBlock* positions, 
					Vector3DBlock* forces,
					ScalarStructure* energies, int from, int to){
    Real e = 0.0;
    Real c = myJ*myK;
    for(int i=from;i<to;i++){
      Vector3D rc = (*positions)[i]-myCenter;
      Real rcNorm = rc.normSquared();
      if(rcNorm > myRadius2){
	rcNorm = sqrt(rcNorm);
	Real a = rcNorm-myRadius;
	(*forces)[i] -= rc*(c*power(a,myJ-1)/rcNorm);
	e += power(a,myJ);
      }
    }
    (*energies)[ScalarStructure::OTHER] += e*myK;

    //     Real e = 0.0;
    //     Real c = myJ*myK*SI_FORCE_FACTOR;
    //     for(int i=from;i<to;i++){
    //       Vector3D rc = (*positions)[i]-myCenter;
    //       Real rcNorm = rc.normSquared();
    //       if(rcNorm > myRadius2){
    // 	rcNorm = sqrt(rcNorm);
    // 	Real m = topo->atoms[i].scaledMmass;
    // 	Real a = (rcNorm-myRadius)/myU;
    // 	(*forces)[i] -= rc*(c*power(a,myJ-1)/rcNorm)*m;
    // 	e += power(a,myJ)*m;
    //       }
    //     }
    //     (*energies)[ScalarStructure::OTHER] += e*myK*SI_ENERGY_FACTOR;
  }

  void SphericalSystemForce::parallelEvaluate(const GenericTopology* topo,
					      const Vector3DBlock* positions, 
					      Vector3DBlock* forces,
					      ScalarStructure* energies){
    unsigned int n = positions->size();
    unsigned int count = numberOfBlocks(topo,positions);
  
    for(unsigned int i = 0;i<count;i++){
      if(Parallel::next()){
	int to = (n*(i+1))/count;
	if(to > (int)n)
	  to = n;
	int from = (n*i)/count;
	doEvaluate(topo,positions,forces,energies,from,to);
      }
    }
  }

  void SphericalSystemForce::getParameters(vector<Parameter>& parameters) const {
    parameters.push_back(Parameter("-center",Value(myCenter)));
    parameters.push_back(Parameter("-radius",Value(myRadius,ConstraintValueType::Positive())));
    parameters.push_back(Parameter("-k",Value(myK,ConstraintValueType::NotZero()),Text("external force constant")));
    parameters.push_back(Parameter("-j",Value(myJ,ConstraintValueType::Positive()),Text("exponent for the potential function")));
    //    parameters.push_back(Parameter("-u",Value(myU,ConstraintValueType::Positive()),1.0,Text("unit factor")));
  }

  Force* SphericalSystemForce::doMake(string& errMsg, vector<Value> values) const {
    int j;

    values[3].get(j);
    if(!values[3].valid() || j < 2 || (j % 2) != 0){
      errMsg += keyword + " force: j (="+values[3].getString()+") > 1 and even.";
      return NULL;
    }
    return new SphericalSystemForce(values[0],values[1],values[2],j/*,values[4]*/);
  }

  unsigned int SphericalSystemForce::numberOfBlocks(const GenericTopology* , 
						    const Vector3DBlock* pos){    
    return std::min(Parallel::getAvailableNum(),static_cast<int>(pos->size()));
  }
}
