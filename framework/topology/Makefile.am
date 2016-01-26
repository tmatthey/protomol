# topology

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libtopology_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# topology Sources
P_TOPOLOGY_SRC =\
ArrayCellListStructure.cpp\
C1SwitchingFunctionBase.cpp\
C1SwitchingFunction.cpp\
C2SwitchingFunctionBase.cpp\
C2SwitchingFunction.cpp\
CnSwitchingFunctionBase.cpp\
CnSwitchingFunction.cpp\
CmpCnCnSwitchingFunctionBase.cpp\
CmpCnCnSwitchingFunction.cpp\
CoulombSCPISMParameterTable.cpp\
buildCellLists.cpp\
ComplementSwitchingFunctionBase.cpp\
CoulombSCPISMParameterTable.cpp\
CubicCellManager.cpp\
CutoffSwitchingFunctionBase.cpp\
CutoffSwitchingFunction.cpp\
ExclusionTable.cpp\
ExclusionType.cpp\
GenericTopology.cpp\
LennardJonesParameterTable.cpp\
PeriodicBoundaryConditions.cpp\
RangeSwitchingFunctionBase.cpp\
ShiftSwitchingFunctionBase.cpp\
ShiftSwitchingFunction.cpp\
topologyutilities.cpp\
UniversalSwitchingFunctionBase.cpp\
VacuumBoundaryConditions.cpp \
AngleInfo.cpp

# topology Includes
P_TOPOLOGY_H =\
Angle.h\
ArrayCellListStructure.h\
Bond.h\
buildCellLists.h\
C1SwitchingFunctionBase.h\
C1SwitchingFunction.h\
C2SwitchingFunctionBase.h\
C2SwitchingFunction.h\
CnSwitchingFunctionBase.h\
CnSwitchingFunction.h\
CmpCnCnSwitchingFunctionBase.h\
CmpCnCnSwitchingFunction.h\
CoulombSCPISMParameters.h\
CoulombSCPISMParameterTable.h\
CellListEnumerator.h\
CellListEnumerator_periodicBoundaries.h\
CellListEnumerator_standard.h\
ComplementSwitchingFunctionBase.h\
ComplementSwitchingFunction.h\
CoulombSCPISMParameterTable.h\
CoulombSCPISMParameters.h\
CubicCellLocation.h\
CubicCellManager.h\
CutoffSwitchingFunctionBase.h\
CutoffSwitchingFunction.h\
ExclusionTable.h\
ExclusionType.h\
GenericTopology.h\
LennardJonesParameters.h\
LennardJonesParameterTable.h\
BankLennardJonesParameterTable.h\
PeriodicBoundaryConditions.h\
RangeSwitchingFunctionBase.h\
RangeSwitchingFunction.h\
SemiGenericTopology.h\
ShiftSwitchingFunctionBase.h\
ShiftSwitchingFunction.h\
Topology.h\
topologyutilities.h\
Torsion.h\
UniversalSwitchingFunctionBase.h\
UniversalSwitchingFunction.h\
VacuumBoundaryConditions.h \
AngleInfo.h


# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libtopology.a

libtopology_a_SOURCES = $(P_TOPOLOGY_SRC) $(P_TOPOLOGY_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY)# -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_TOPOLOGY_SRC)

include $(top_srcdir)/Makefile.depend
