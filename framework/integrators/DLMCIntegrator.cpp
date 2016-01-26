#include "DLMCIntegrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "ForceGroup.h"
#include "Force.h"

using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
  //____________________________________________________________ DLMCIntegrator

  const string DLMCIntegrator::keyword("DihedralLiftMC");

  const int DLMCIntegrator::myNumParameters(2);

  //  ---------------------------------------------------------------------  //

  DLMCIntegrator::DLMCIntegrator()    
    : MCIntegrator() {
    hdf = 0;
    lift = false; // The system will run normal MD first    
  }

  //  ---------------------------------------------------------------------  //

  DLMCIntegrator::DLMCIntegrator(int cycles,
			       Real initialTemperature,
			       ForceGroup *overloadedForces,
			       StandardIntegrator *nextIntegrator)
    
    : MCIntegrator(cycles, initialTemperature, overloadedForces,
            nextIntegrator) {
    hdf = 0;
    lift = false; // The system will run normal MD first
  }
  
  //  ---------------------------------------------------------------------  //

  void DLMCIntegrator::initialize(GenericTopology *topo,
				 Vector3DBlock *positions,
				 Vector3DBlock *velocities,
				 ScalarStructure *energies){
    
    MCIntegrator::initialize(topo,positions,velocities,energies);
    
    std::vector<Force*> myForcePtrs = next()->getForceGroup()->getForces();
    //std::cout << "myForce Ptrs size: " << myForcePtrs.size();
    std::vector<Force*>::iterator forceItr = myForcePtrs.begin();
    while( (forceItr != myForcePtrs.end()) && (hdf == 0))
    {
      std::cout << (*forceItr)->getId() << std::endl;
      if(equalNocase((*forceItr)->getId(),"HarmDihedral"))
        hdf = (*forceItr);
      forceItr++;
    }
    if (hdf == 0)
    {
      std::cout << "A HarmDihedralSystemForce pointer was not found" << std::endl;      
      report << debug(1) << "A HarmDihedralSystemForce pointer was not found" << endr;

    }
    else
    {
      std::cout << "hdf points to a force of type: ";
      std::cout << hdf->getId() << std::endl;
      std::cout.flush();
      report << debug(1) << hdf->getId() << ":" << hdf->getParameters()[1].value << endr;
    }
  }

  //  ---------------------------------------------------------------------  //

  MTSIntegrator*  DLMCIntegrator::doMake(string& , const vector<Value>& values, 
					ForceGroup* fg, StandardIntegrator *nextIntegrator)const{

      return new DLMCIntegrator(values[0],values[1], fg,nextIntegrator);

  }

  //  ---------------------------------------------------------------------  //
  
  void DLMCIntegrator::perturbSystem() {

    Real angle = randomNumber()*M_PI;
    if (randomNumber() < 0.50)
      angle = angle * -1;
    
    std::string err;    
    if (lift)
    {
      std::cout << "Trial angle for lift: " << angle << std::endl;      
      hdf->setParameter(err,0,0);
      hdf->setParameter(err,1,angle);      
      report << debug(1) << hdf->getId() << ":" << hdf->getParameters()[0].value
             << ":" << hdf->getParameters()[1].value << endr;      
      lift = !(lift);
    }
    else
    {
      hdf->setParameter(err,0,-1);
      report << debug(1) << hdf->getId() << ":" << hdf->getParameters()[0].value
             << ":" << hdf->getParameters()[1].value << endr;      
      lift = !(lift);
    }
    
    //randomVelocity(getInitialTemperature(), myTopo,myVelocities);
    buildMolecularMomentum(myVelocities,myTopo);
  }

  //  ---------------------------------------------------------------------  //
  
  void DLMCIntegrator::saveValues(){
    MCIntegrator::saveValues();
  }

  //  ---------------------------------------------------------------------  //
  
  void DLMCIntegrator::restoreValues(){
    MCIntegrator::restoreValues();   
  }

}
