#include "MTSIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
using std::vector;
using std::string;

namespace ProtoMol {
  //_________________________________________________________________ MTSIntegrator

  MTSIntegrator::MTSIntegrator():StandardIntegrator(),myNextIntegrator(NULL),myCycleLength(0){}

  MTSIntegrator::MTSIntegrator (int cycles,
				ForceGroup *overloadedForces,
				StandardIntegrator *nextIntegrator)
    : StandardIntegrator(overloadedForces),
      myNextIntegrator(nextIntegrator),
      myCycleLength(cycles) {
  
    // The plan for this constructor:
    //  We make sure that we get a non-zero pointer to ForceGroup object.
    //  In some cases (for test purpose) we may not have any forces to evaluate, 
    //  but the force evaluation is still call ...
  
    myNextIntegrator->myPreviousIntegrator = this;
  }


  MTSIntegrator::~MTSIntegrator(){
    delete myNextIntegrator;
  }


  void MTSIntegrator::doDriftOrNextIntegrator() {
    preDriftOrNextModify();
    myNextIntegrator->run(myCycleLength);
    postDriftOrNextModify();
  }


  void MTSIntegrator::initialize(GenericTopology *topo,
				 Vector3DBlock *positions,
				 Vector3DBlock *velocities, 
				 ScalarStructure *energies){
    myNextIntegrator->initialize(topo,positions,velocities,energies);
    StandardIntegrator::initialize(topo,positions,velocities,energies);
    //Report::report <<"[MTSIntegrator::initialize]"<<Report::endr;
  }

  void MTSIntegrator::getParameters(vector<Parameter>& parameters) const {
    parameters.push_back(Parameter("cyclelength",Value(myCycleLength,ConstraintValueType::Positive())));
  }

  MTSIntegrator*  MTSIntegrator::make(string& errMsg, const vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    errMsg = "";
    if(!checkParameters(errMsg,values))
      return NULL;
    return adjustAlias(doMake(errMsg,values,fg,nextIntegrator));
  }

}
