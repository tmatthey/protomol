#include "Output.h"
#include "OutputCache.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "inputValueDefinitions.h"
#include "ScalarStructure.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

namespace ProtoMol {
  //________________________________________________________ Output
  const string  Output::scope("Output");

  Output::Output():Makeable(),
		   myFirstStep(0),
		   myLastStep(0),
		   myNextStep(0),
		   myFirst(true),
		   myOutputFreq(0),
		   myConfig(NULL), 
		   myTopology(NULL), 
		   myIntegrator(NULL),
		   myEnergies(NULL),
		   myPositions(NULL),
		   myVelocities(NULL),
		   myCache(NULL){}

    Output::Output(int freq):Makeable(),
			     myFirstStep(0),
			     myLastStep(0),
			     myNextStep(0),
			     myFirst(true),
			     myOutputFreq(freq),
			     myConfig(NULL), 
			     myTopology(NULL), 
			     myIntegrator(NULL),
			     myEnergies(NULL),
			     myPositions(NULL),
			     myVelocities(NULL),
			     myCache(NULL){}

  void Output::initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
			  const Vector3DBlock* pos, const Vector3DBlock*  vel, const ScalarStructure* energies){
    myFirst      = true;
    myConfig     = config; 
    myTopology   = topo; 
    myIntegrator = integrator;
    myEnergies   = energies;  
    myPositions  = pos;       
    myVelocities = vel;       
    if(config->valid(InputFirststep::keyword)){
      myNextStep   = (*config)[InputFirststep::keyword];
      myFirstStep  = (*config)[InputFirststep::keyword];
      myLastStep   = myFirstStep;
    }
    if(config->valid(InputNumsteps::keyword)){
      myLastStep  = myLastStep + (*config)[InputNumsteps::keyword].operator int();
    }
    doInitialize();
  }

  void Output::updateNextStep(int step) {
    int n = (step - myNextStep)/myOutputFreq;
    myNextStep += std::max(n,1)*myOutputFreq;
  }

  void Output::run(int step){
    if(step >= myNextStep){
      int n = (step - myNextStep)/myOutputFreq;
      myNextStep += std::max(n,1)*myOutputFreq;
      if(myEnergies->output())
	doRun(step);      
    }
    myFirst = false;
  }

  void Output::finalize(int step){
    if(myEnergies->output())
      doRun(step);
    doFinalize(step);
  }

  Output* Output::make(string& errMsg, const vector<Value>& values) const{
    errMsg = "";
    if(!checkParameters(errMsg,values))
      return NULL;
    return adjustAlias(doMake(errMsg,values));
  }


  void Output::setCache(const OutputCache* cache){
    myCache = cache;
  }

  bool Output::isIdDefined(const Configuration* config) const{
    if(!addDoKeyword())
      return config->valid(getId());

    string str("do"+getId());
    if(!config->valid(getId()) || config->empty(str))
      return false;

    if(!config->valid(str))
      return true;
    return (bool)(*config)[str];
  }

}
