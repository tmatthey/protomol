#include "HapticSystemForce.h"
#include "IMDElf.h"
#include "topologyutilities.h"
#include "Vector3DBlock.h"
//#include "ScalarStructure.h"
#include "GenericTopology.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ HapticSystemForce
  const string HapticSystemForce::keyword("Haptic");

  HapticSystemForce::HapticSystemForce():SystemForce(),
					 myHapticForces(NULL),
					 myIMDElf(NULL),
					 myStepInc(-1),
					 myStepCounter(-1),
					 myPort(-1),
					 myCached(false){}

  HapticSystemForce::HapticSystemForce(int port, int trate, int timeout, 
				       int stepInc, int waitIMD):SystemForce(),
								 myHapticForces(new Vector3DBlock()),
								 myIMDElf(new IMDElf()),
								 myStepInc(0),
								 myStepCounter(stepInc),
								 myPort(port),
								 myCached(false)
  {
    myIMDElf->trate    = trate;
    myIMDElf->timeout  = timeout;
    myIMDElf->wait_imd = waitIMD;

    if (myIMDElf->listen(port))
      report << recoverable <<"[HapticSystemForce::HapticSystemForce] IMD listen failed!"<<endr;
    if (myIMDElf->hookup())
      report << recoverable <<"[HapticSystemForce::HapticSystemForce] IMD connection attempt failed!"<<endr;
    else if (waitIMD)    
      report << plain << "IMD2 connection ESTABLISHED." << endr;

  }

//   HapticSystemForce::HapticSystemForce(const HapticSystemForce& haptic){
//     myIMDElf->trate    = haptic.myIMDElf->trate;
//     myIMDElf->timeout  = haptic.myIMDElf->timeout;
//     myIMDElf->wait_imd = haptic.myIMDElf->wait_imd;
//     myStepCounter  = 0;
//     myStepInc      = haptic.myStepInc;
//     myHapticForces = new Vector3DBlock();
//   }

  HapticSystemForce::~HapticSystemForce(){
    if(myHapticForces != NULL)
      delete myHapticForces;
    if(myIMDElf != NULL)
      delete myIMDElf;
  }

  void HapticSystemForce::setForce(Vector3DBlock* forces){
    if(myHapticForces != forces && myHapticForces != NULL)
      delete myHapticForces;
    myHapticForces = forces;
  }


  void HapticSystemForce::getParameters(vector<Parameter>& parameters) const {
    parameters.push_back(Parameter("-port",Value(myPort,ConstraintValueType::NotNegative()),2000));
    parameters.push_back(Parameter("-trate",Value(myIMDElf!=NULL?myIMDElf->trate:-1,ConstraintValueType::NotNegative()),1));
    parameters.push_back(Parameter("-timeout",Value(myIMDElf!=NULL?myIMDElf->timeout:-1,ConstraintValueType::NotNegative()),1000));
    parameters.push_back(Parameter("-step_inc", Value(myStepInc,ConstraintValueType::NotNegative()),1));
    parameters.push_back(Parameter("-wait",Value(myIMDElf!=NULL?myIMDElf->wait_imd:-1,ConstraintValueType::NotNegative()),0));
  }

  Force* HapticSystemForce::doMake(string& , vector<Value> values) const {
    return (new HapticSystemForce(values[0],values[1],values[2],values[3],values[4]));
  }

  void HapticSystemForce::evaluate (const GenericTopology* topo,
				    const Vector3DBlock* positions,
				    Vector3DBlock* forces,
				    ScalarStructure* ) {
    // Set the correct size of the haptic forces coordinate block
    myHapticForces->zero(positions->size());

    //Check for new connections if client is detached
    if (!myIMDElf->client_alive())
      myIMDElf->hookup();
  
    // UPDATE theData.theTopo.hapticForces HERE!!!
    if (myIMDElf->client_alive()) {
      if (myIMDElf->get_haptics(*myHapticForces))
	report << recoverable << "[HapticSystemForce::evaluate] Updating haptic forces for current time-step." << endr;
    }

    //Center the coordinates to fit 
    //within a VMD window.
    Vector3DBlock data = *positions;
    if(!myCached){
      myCenter = centerOfMass(positions,topo);
      myCached = true;
      data = *positions;
    }
    else {
      Vector3D delta = myCenter-centerOfMass(positions,topo);
      for(unsigned int i=0;i<positions->size();i++)
	data[i] = (*positions)[i]+delta;
    }
    // Attempt to send imd coordinates
    if (myIMDElf->client_alive())
      if ((myStepCounter % myIMDElf->trate) == 0)
	if (myIMDElf->send_coords(data) < 0)
	  report << recoverable << "[HapticSystemForce::evaluate] Sending coordinates to client for time-step: " << myStepCounter<<"." << endr;
  
    forces->intoAdd(*myHapticForces);
    myStepCounter += myStepInc;

    //  energies.hapticEnergy = myHapticEnergy;
  }
}
