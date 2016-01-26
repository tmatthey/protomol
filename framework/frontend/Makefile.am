# frontend

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libfrontend_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# frontend Sources
P_FRONTEND_H =\
HelpTextFactory.h\
Factory.h\
ForceFactory.h\
ForceFactoryDetails.h\
InputPosVel.h\
InputPosVelType.h\
InputValue.h\
IntegratorFactory.h\
IntegratorFactoryDetails.h\
Output.h\
OutputBSDL.h\
OutputCache.h\
OutputCollection.h\
OutputDCDTrajectory.h\
OutputDCDTrajectoryNoWater.h\
OutputDCDTrajectoryVel.h\
OutputDiffusion.h\
OutputDihedrals.h\
OutputEnergies.h\
OutputFactory.h\
OutputFactoryDetails.h\
OutputFile.h\
OutputFinalXYZBinPos.h\
OutputFinalXYZBinVel.h\
OutputFinalXYZBinVelRev.h\
OutputFinalXYZPos.h\
OutputFinalPDBPos.h\
OutputFinalXYZVel.h\
OutputMomentum.h\
OutputPaulTrap.h\
OutputPDBFramePos.h\
OutputREMExchangeRate.h\
OutputREMHistory.h\
OutputScreen.h\
OutputState.h\
OutputTemperatures.h\
OutputXYZTrajectoryForce.h\
OutputXYZTrajectoryPos.h\
OutputXYZTrajectoryVel.h\
OutputXYZBinTrajectoryPos.h\
OutputXYZBinTrajectoryVel.h\
TopologyFactory.h\
TopologyFactoryDetails.h\
buildTopology.h\
inputValueDefinitions.h\
parseCommandLine.h

# frontend Includes
P_FRONTEND_SRC =\
HelpTextFactory.cpp\
ForceFactoryDetails.cpp\
InputPosVel.cpp\
InputPosVelType.cpp\
IntegratorFactoryDetails.cpp\
Output.cpp\
OutputBSDL.cpp\
OutputCache.cpp\
OutputCollection.cpp\
OutputDCDTrajectory.cpp\
OutputDCDTrajectoryNoWater.cpp\
OutputDCDTrajectoryVel.cpp\
OutputDiffusion.cpp\
OutputDihedrals.cpp\
OutputEnergies.cpp\
OutputFactoryDetails.cpp\
OutputFile.cpp\
OutputFinalXYZBinPos.cpp\
OutputFinalXYZBinVel.cpp\
OutputFinalXYZBinVelRev.cpp\
OutputFinalXYZPos.cpp\
OutputFinalPDBPos.cpp\
OutputFinalXYZVel.cpp\
OutputMomentum.cpp\
OutputPaulTrap.cpp\
OutputPDBFramePos.cpp\
OutputREMExchangeRate.cpp\
OutputREMHistory.cpp\
OutputScreen.cpp\
OutputState.cpp\
OutputTemperatures.cpp\
OutputXYZTrajectoryForce.cpp\
OutputXYZTrajectoryPos.cpp\
OutputXYZTrajectoryVel.cpp\
OutputXYZBinTrajectoryPos.cpp\
OutputXYZBinTrajectoryVel.cpp\
TopologyFactoryDetails.cpp\
buildTopology.cpp\
inputValueDefinitions.cpp\
parseCommandLine.cpp

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libfrontend.a

libfrontend_a_SOURCES = $(P_FRONTEND_SRC) $(P_FRONTEND_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND)# -I$(P_FACTORIES)

SOURCEFILES = $(P_FRONTEND_SRC)

include $(top_srcdir)/Makefile.depend
