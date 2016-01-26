CXX=g++  # CHANGE TO icpc if necessary
CFILE = main.cpp \
scpenergy.cpp \
AngleInfo.cpp \
ArrayCellListStructure.cpp \
buildExclusionTable.cpp \
C1SwitchingFunctionBase.cpp \
C1SwitchingFunction.cpp \
C2SwitchingFunctionBase.cpp \
C2SwitchingFunction.cpp \
CmpCnCnSwitchingFunctionBase.cpp \
CmpCnCnSwitchingFunction.cpp \
CnSwitchingFunctionBase.cpp \
CnSwitchingFunction.cpp \
CompareForce.cpp \
ConstraintValueType.cpp \
CoulombBornRadiiForce.cpp \
CoulombSCPISMForce.cpp \
CoulombSCPISMParameterTable.cpp \
CubicCellManager.cpp \
CutoffSwitchingFunctionBase.cpp \
CutoffSwitchingFunction.cpp \
ExclusionTable.cpp \
ExclusionType.cpp \
ExtendedCompareForce.cpp \
ExtendedForce.cpp \
ExtendedTimeForce.cpp \
Force.cpp \
ForceGroup.cpp \
GenericTopology.cpp \
LennardJonesParameterTable.cpp \
Makeable.cpp \
Matrix3by3.o \
MollyForce.cpp \
NonbondedCutoffForceBase.cpp \
Parallel.cpp \
Parameter.cpp \
Report.cpp \
ScalarStructure.cpp \
ShiftSwitchingFunctionBase.cpp \
ShiftSwitchingFunction.cpp \
SystemCompareForce.cpp \
SystemForce.cpp \
SystemTimeForce.cpp \
TimeForce.cpp \
Timer.cpp \
TimerStatistic.cpp \
VacuumBoundaryConditions.cpp \
Value.cpp \
ValueType.cpp \
Vector3DBlock.cpp \
buildCellLists.cpp \
mathutilities.cpp \
pmconstants.cpp \
stringutilities.cpp \
systemutilities.cpp \
topologyutilities.cpp
OFILES=${CFILE:.cpp=.o}
CPPFLAGS = -g

scpenergy: ${OFILES}
	$(CXX) $(OFILES) -o scpenergy 

%.o:%.cpp
	${CXX} -c ${CPPFLAGS} ${INCLUDE} $<

clean:
	rm -f *.o
	rm -f scpenergy core
