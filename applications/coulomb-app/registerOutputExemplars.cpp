#include "registerOutputExemplars.h"
#include "OutputDCDTrajectory.h"
#include "OutputFactory.h"
//#include "OutputFinalXYZBinPos.h"
//#include "OutputFinalXYZBinVel.h"
#include "OutputFinalXYZPos.h"
#include "OutputFinalXYZVel.h"
#include "OutputPaulTrap.h"
#include "OutputXYZTrajectoryPos.h"

namespace ProtoMol {

  void registerOutputExemplars(){
    OutputFactory::registerExemplar(new OutputDCDTrajectory());
//    OutputFactory::registerExemplar(new OutputFinalXYZBinPos());
//    OutputFactory::registerExemplar(new OutputFinalXYZBinVel());
    OutputFactory::registerExemplar(new OutputFinalXYZPos());
    OutputFactory::registerExemplar(new OutputFinalXYZVel());
    OutputFactory::registerExemplar(new OutputPaulTrap());
    OutputFactory::registerExemplar(new OutputXYZTrajectoryPos());
  }
}
