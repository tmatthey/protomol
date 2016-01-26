#include "oSGregisterOutputExemplars.h"
#include "OutputEnergies.h"
#include "OutputOSGProperties.h"
#include "OutputFactory.h"  
#include "OutputFinalXYZBinPos.h"
#include "OutputFinalXYZBinVel.h"
#include "OutputFinalXYZPos.h"
#include "OutputFinalXYZVel.h"
#include "OutputMomentum.h"
#include "OutputPDBFramePos.h"
#include "OutputScreen.h"

#include "OutputOSGDCDTrajectory.h"
#include "OutputXYZTrajectoryForce.h"
#include "OutputOSGXYZTrajectoryPos.h"
#include "OutputXYZTrajectoryVel.h"

#include "OutputFinalXSC.h"
#include "OutputFinalPSF.h"

namespace ProtoMol {

  void oSGregisterOutputExemplars(){
    OutputFactory::registerExemplar(new OutputEnergies());
    OutputFactory::registerExemplar(new OutputFinalXYZBinPos());
    OutputFactory::registerExemplar(new OutputFinalXYZBinVel());
    OutputFactory::registerExemplar(new OutputFinalXYZPos());
    OutputFactory::registerExemplar(new OutputFinalXYZVel());
    OutputFactory::registerExemplar(new OutputMomentum());
    OutputFactory::registerExemplar(new OutputPDBFramePos());
    OutputFactory::registerExemplar(new OutputScreen());
    
    OutputFactory::registerExemplar(new OutputOSGDCDTrajectory());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryForce());
    OutputFactory::registerExemplar(new OutputOSGXYZTrajectoryPos());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryVel());
    
    OutputFactory::registerExemplar(new OutputOSGProperties());
    OutputFactory::registerExemplar(new OutputFinalXSC());
    OutputFactory::registerExemplar(new OutputFinalPSF());
  }
}
