# integrators

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libintegrators_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories


# integrators Sources
P_INTEGRATORS_SRC =\
BBKIntegrator.cpp\
BSplineMOLLYIntegrator.cpp\
CGMinimizerIntegrator.cpp\
DLMCIntegrator.cpp\
DMDLeapfrogIntegrator.cpp\
DihedralHMCIntegrator.cpp\
EquilibriumMOLLYIntegrator.cpp\
HMCIntegrator.cpp\
HessianInt.cpp \
ImpulseIntegrator.cpp\
Integrator.cpp\
IntegratorDefinition.cpp\
LangevinImpulseIntegrator.cpp\
LangevinLeapfrogIntegrator.cpp\
LeapfrogIntegrator.cpp\
LeapfrogTruncatedShadow.cpp\
MCIntegrator.cpp\
MOLLYIntegrator.cpp\
MTSIntegrator.cpp\
ModifierMetaRattle.cpp\
ModifierMetaRattleShake.cpp\
ModifierMetaShake.cpp\
ModifierNPTRattle.cpp\
ModifierNPTShake.cpp\
ModifierNVTRattle.cpp\
ModifierNVTShake.cpp\
ModifierRattle.cpp\
ModifierShadow.cpp\
ModifierShake.cpp\
NPTVerletIntegrator.cpp\
NVTVerletIntegrator.cpp\
NormModeInt.cpp\
NormModeMin.cpp\
NormModeSmplMin.cpp\
NormModeDiag.cpp\
NormModeVisual.cpp\
NoseNVTLeapfrogIntegrator.cpp\
RMTIntegrator.cpp\
BerendsenIntegrator.cpp\
PLeapfrogIntegrator.cpp\
PaulTrapIntegrator.cpp\
STSIntegrator.cpp\
ShadowHMCIntegrator.cpp\
StandardIntegrator.cpp\
Hessian.cpp\
S2HMCIntegrator.cpp\
NormalModeUtilities.cpp\
NormalModeLangevin.cpp\
NormalModeMinimizer.cpp\
NormalModeDiagonalize.cpp\
NormalModeMori.cpp\
NormalModeRelax.cpp\
NormalModeBrownian.cpp\
NumericalDifferentiation.cpp

# integrators Includes
P_INTEGRATORS_H =\
BBKIntegrator.h\
BSplineMOLLYIntegrator.h\
CGMinimizerIntegrator.h\
DLMCIntegrator.h\
DMDLeapfrogIntegrator.h\
DihedralHMCIntegrator.h\
EquilibriumMOLLYIntegrator.h\
HMCIntegrator.h\
HessianInt.h\
ImpulseIntegrator.h\
Integrator.h\
IntegratorDefinition.h\
LangevinImpulseIntegrator.h\
LangevinLeapfrogIntegrator.h\
LeapfrogIntegrator.h\
LeapfrogTruncatedShadow.h\
MCIntegrator.h\
MOLLYIntegrator.h\
MTSIntegrator.h\
Modifier.h\
ModifierAveraging.h\
ModifierFriction.h\
ModifierIncrementTimestep.h\
ModifierMetaRattle.h\
ModifierMetaRattleShake.h\
ModifierMetaShake.h\
ModifierMollification.h\
ModifierNPTRattle.h\
ModifierNPTShake.h\
ModifierNVTRattle.h\
ModifierNVTShake.h\
ModifierPostForceBarostat.h\
ModifierPostForceThermostat.h\
ModifierPreForceBarostat.h\
ModifierPreForceThermostat.h\
ModifierRattle.h\
ModifierRemoveAngularMomentum.h\
ModifierRemoveLinearMomentum.h\
ModifierShadow.h\
ModifierShake.h\
ModifierBerendsen.h\
NPTVerletIntegrator.h\
NVTVerletIntegrator.h\
NonStandardIntegrator.h\
NormModeInt.h\
NormModeMin.h\
NormModeSmplMin.h\
NormModeDiag.h\
NormModeVisual.h\
NoseNPTLeapfrogIntegrator.h\
NoseNVTLeapfrogIntegrator.h\
RMTIntegrator.h\
BerendsenIntegrator.h\
PLeapfrogIntegrator.h\
PaulTrapIntegrator.h\
SimTKlapack.h\
STSIntegrator.h\
ShadowHMCIntegrator.h\
StandardIntegrator.h\
Hessian.h\
S2HMCIntegrator.h\
NormalModeUtilities.h\
NormalModeLangevin.h\
NormalModeMinimizer.h\
NormalModeDiagonalize.h\
NormalModeMori.h\
NormalModeRelax.h\
NormalModeBrownian.h\
NumericalDifferentiation.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libintegrators.a

libintegrators_a_SOURCES = $(P_INTEGRATORS_SRC) $(P_INTEGRATORS_H)

INCLUDES = -I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) \
           -I$(P_IMD) -I$(P_INTEGRATORS) @LAPACK_CFLAGS@

SOURCEFILES = $(P_INTEGRATORS_SRC)

include $(top_srcdir)/Makefile.depend
