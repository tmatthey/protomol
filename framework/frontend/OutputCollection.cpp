#include "OutputCollection.h"
#include "Output.h"
#include "OutputCache.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "Integrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "topologyutilities.h"

namespace ProtoMol {
  //________________________________________________________ OutputCollection
  OutputCollection::OutputCollection():myCache(new OutputCache()),
				       myConfig(NULL), 
				       myTopology(NULL), 
				       myIntegrator(NULL),
				       myEnergies(NULL),
				       myPositions(NULL),
				       myVelocities(NULL){}

  OutputCollection::~OutputCollection(){
    for(iterator  i=begin();i!=end();++i)
      delete (*i);    
    delete myCache;
  }

  void OutputCollection::initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
				    const Vector3DBlock* pos, const Vector3DBlock*  vel, const ScalarStructure* energies){
    myConfig     = config; 
    myTopology   = topo; 
    myIntegrator = integrator;
    myEnergies   = energies;  
    myPositions  = pos;       
    myVelocities = vel;       
    myCache->initialize(config, integrator, topo, pos, vel, energies);
    for(iterator  i=begin();i!=end();++i)
      (*i)->initialize(config, integrator, topo, pos, vel, energies);
  }
  
  void OutputCollection::run(int step){
    myCache->uncache();
    for(iterator  i=begin();i!=end();++i)
      (*i)->run(step);
  }

  void OutputCollection::updateNextStep(int step){
    myCache->uncache();
    for(iterator  i=begin();i!=end();++i)
      (*i)->updateNextStep(step);
  }

  void OutputCollection::restoreState(CheckpointInputStream& is) {
    for (iterator i = begin(); i != end(); i++)
      (*i)->restoreState(is);
  }

  void OutputCollection::finalize(int step){
    myCache->uncache();
    for(iterator  i=begin();i!=end();++i)
      (*i)->finalize(step);    
  }
  
  void OutputCollection::adoptOutput(Output* output){
    if(output != NULL){
      output->setCache(myCache);
      myOutputList.push_back(output);	
    }
  }

  int  OutputCollection::getNext() const{
    int next = Constant::MAX_INT;
    for(const_iterator  i=begin();i!=end();++i)
      next = std::min((*i)->getNext(),next);
    return next;
  }

}
