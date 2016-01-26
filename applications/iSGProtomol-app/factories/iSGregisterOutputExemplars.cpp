#include "iSGregisterOutputExemplars.h"
#include "OutputDCDTrajectory.h"
#include "OutputEnergies.h"
#include "OutputISGProperties.h"
#include "OutputFactory.h"  
#include "OutputFinalXYZBinPos.h"
#include "OutputFinalXYZBinVel.h"
#include "OutputFinalXYZPos.h"
#include "OutputFinalXYZVel.h"
#include "OutputMomentum.h"
#include "OutputPDBFramePos.h"
#include "OutputScreen.h"
#include "OutputXYZTrajectoryForce.h"
#include "OutputXYZTrajectoryPos.h"
#include "OutputXYZTrajectoryVel.h"
#include "OutputFinalXSC.h"
#include "OutputFinalPSF.h"

namespace ProtoMol {

  void iSGregisterOutputExemplars(){
    OutputFactory::registerExemplar(new OutputDCDTrajectory());
    OutputFactory::registerExemplar(new OutputEnergies());
    OutputFactory::registerExemplar(new OutputFinalXYZBinPos());
    OutputFactory::registerExemplar(new OutputFinalXYZBinVel());
    OutputFactory::registerExemplar(new OutputFinalXYZPos());
    OutputFactory::registerExemplar(new OutputFinalXYZVel());
    OutputFactory::registerExemplar(new OutputMomentum());
    OutputFactory::registerExemplar(new OutputPDBFramePos());
    OutputFactory::registerExemplar(new OutputScreen());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryForce());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryPos());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryVel());
    OutputFactory::registerExemplar(new OutputISGProperties());
    OutputFactory::registerExemplar(new OutputFinalXSC());
    OutputFactory::registerExemplar(new OutputFinalPSF());
  }
}
