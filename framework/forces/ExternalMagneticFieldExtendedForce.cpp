#include "ExternalMagneticFieldExtendedForce.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"
#include "Parallel.h"
#include "ScalarStructure.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {
  //_________________________________________________________________ ExternalMagneticFieldExtendedForce


 
  const string ExternalMagneticFieldExtendedForce::keyword("ExternalMagneticField");
  static const Real SI_FORCE_FACTOR = Constant::SI::ELECTRON_CHARGE*Constant::SI::AVOGADRO*Constant::SI::KCAL/Constant::SI::LENGTH_AA*(1e8/sqrt(4184.0));

  ExternalMagneticFieldExtendedForce::ExternalMagneticFieldExtendedForce():ExtendedForce(),myB(Vector3D(0.0,0.0,0.0)){}

 
  ExternalMagneticFieldExtendedForce::ExternalMagneticFieldExtendedForce(const Vector3D& b ):ExtendedForce(),myB(b){}

 
  void ExternalMagneticFieldExtendedForce::getParameters(vector<Parameter>& parameters) const {
    parameters.push_back(Parameter("-B",Value(myB),Text("External magnetic field [T]")));
  }

 
  Force* ExternalMagneticFieldExtendedForce::doMake(string&, vector<Value> values) const {
    return new ExternalMagneticFieldExtendedForce(values[0]);
  }
 
  void ExternalMagneticFieldExtendedForce::evaluate(const GenericTopology* topo,
						    const Vector3DBlock* positions, 
						    const Vector3DBlock* velocities,
						    Vector3DBlock* forces,
						    ScalarStructure* energies){
    doEvaluate(topo,positions,velocities,forces,energies,0,positions->size());
  }

 
  void ExternalMagneticFieldExtendedForce::parallelEvaluate(const GenericTopology* topo,
							    const Vector3DBlock* positions, 
							    const Vector3DBlock* velocities,
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
	doEvaluate(topo,positions,velocities,forces,energies,from,to);
      }
    }
  }

 
  void ExternalMagneticFieldExtendedForce::doEvaluate(const GenericTopology* topo,
						      const Vector3DBlock* positions, 
						      const Vector3DBlock* velocities,
						      Vector3DBlock* forces,
						      ScalarStructure* /*energies*/,
						      int from, int to){
    Vector3D b(myB*SI_FORCE_FACTOR);
    Real e = 0.0;
    for(int i=from;i<to;i++){
      Vector3D f((*velocities)[i].cross(b)*topo->atomTypes[topo->atoms[i].type].charge);
      //report << (*forces)[i] <<",";
      (*forces)[i] += f;
      //report << b <<","<< SI_FORCE_FACTOR<<",";
      //report << (*velocities)[i].cross(b)*topo->atomTypes[topo->atoms[i].type].charge <<endr;
      e += f*(*positions)[i];
    }
    //(*energies)[ScalarStructure::OTHER] += e;
    
  }

 
  unsigned int ExternalMagneticFieldExtendedForce::numberOfBlocks(const GenericTopology*, 
								  const Vector3DBlock* pos){    
    return std::min(static_cast<unsigned int>(Parallel::getAvailableNum()),static_cast<unsigned int>(pos->size()));
  }

}  
