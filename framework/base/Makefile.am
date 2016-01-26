# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libbase_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# base Sources
P_BASE_SRC =\
BSpline.cpp\
BSplineType.cpp\
Configuration.cpp\
ConstraintValueType.cpp\
FFTComplex.cpp\
Lagrange.cpp\
Makeable.cpp\
Matrix3by3.cpp\
Hermite.cpp\
iSGPAR.cpp\
PAR.cpp\
PDB.cpp\
PPM.cpp\
PGM.cpp\
PSF.cpp\
Parameter.cpp\
ReducedHessAngle.cpp\
reducedHessBond.cpp\
Report.cpp\
ScalarStructure.cpp\
StateRestore.cpp\
Timer.cpp\
TimerStatistic.cpp\
Value.cpp\
ValueType.cpp\
Vector3DBlock.cpp\
XYZ.cpp\
XSC.cpp\
pmconstants.cpp\
fft.c\
mathutilities.cpp\
stringutilities.cpp\
systemutilities.cpp

# base Includes
P_BASE_H =\
AbstractEnumType.h\
Array.h\
Atom.h\
AtomType.h\
BSpline.h\
BSplineType.h\
Configuration.h\
ConstraintValueType.h\
EigenvectorInfo.h\
FFTComplex.h\
Grid.h\
Lagrange.h\
LookUpTableBase.h\
Makeable.h\
MakeableDefinition.h\
Matrix3by3.h\
Molecule.h\
MultiGrid.h\
Hermite.h\
PAR.h\
PDB.h\
PGM.h\
PPM.h\
PSF.h\
Parameter.h\
Proxy.h\
Real.h\
ReducedHessAngle.h\
Report.h\
ScalarStructure.h\
Stack.h \
StateRestore.h \
Timer.h\
TimerStatistic.h\
Value.h\
ValueType.h\
Vector.h\
Vector3D.h\
Vector3DBlock.h\
XYZ.h\
XSC.h\
iSGPAR.h\
mathutilities.h\
pmconstants.h\
protomol.h\
reducedHessBond.h\
simpleTypes.h\
stringutilities.h\
systemutilities.h\
typeSelection.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libbase.a

libbase_a_SOURCES = $(P_BASE_SRC) $(P_BASE_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES) @BOOST_CFLAGS@

SOURCEFILES = $(P_BASE_SRC)

include $(top_srcdir)/Makefile.depend
