# forces

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libforces_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# forces Sources
P_FORCES_SRC =\
AngleSystemForceBase.cpp\
BondSystemForceBase.cpp\
CompareForce.cpp\
CoulombBornRadiiForce.cpp\
CoulombEwaldRealForce.cpp\
CoulombEwaldRealTableForceBase.cpp\
CoulombForce.cpp\
CoulombForceDiElec.cpp\
CoulombMultiGridDirectForceBase.cpp\
CoulombMultiGridDirectTableForceBase.cpp\
CoulombSCPISMForce.cpp\
CoulombTableForceBase.cpp\
DihedralSystemForceBase.cpp\
ElectricFieldSystemForceBase.cpp\
ExtendedCompareForce.cpp\
ExtendedForce.cpp\
ExtendedTimeForce.cpp\
ExternalGravitationSystemForce.cpp\
ExternalMagneticFieldExtendedForce.cpp\
Force.cpp\
ForceGroup.cpp\
FrictionExtendedForce.cpp\
GravitationForce.cpp\
HarmDihedralSystemForceBase.cpp\
ImproperSystemForceBase.cpp\
LennardJonesForce.cpp\
LennardJonesTableForceBase.cpp\
MagneticDipoleForce.cpp\
MagneticDipoleMirrorSystemForceBase.cpp\
MetaForce.cpp\
MollyForce.cpp\
NonbondedCutoffForceBase.cpp\
NonbondedFullEwaldSystemForceBase.cpp\
NonbondedFullSystemForceBase.cpp\
NonbondedMultiGridSystemForceBase.cpp\
NonbondedPMEwaldSystemForceBase.cpp\
NonbondedSimpleFullSystemForceBase.cpp\
PaulTrapExtendedForceBase.cpp\
ReducedHessCoulomb.cpp\
ReducedHessCoulombBornRadii.cpp\
ReducedHessCoulombDiElec.cpp\
ReducedHessCoulombSCPISM.cpp\
ReducedHessLennardJones.cpp\
SphericalRestraintSystemForce.cpp\
SphericalSystemForce.cpp\
SystemCompareForce.cpp\
SystemForce.cpp\
SystemTimeForce.cpp\
TimeForce.cpp\
WrapperMetaForce.cpp\
HessDihedral.cpp

# forces Includes
P_FORCES_H =\
AngleSystemForce.h\
AngleSystemForceBase.h\
BondSystemForce.h\
BondSystemForceBase.h\
CompareForce.h\
CoulombBornRadiiForce.h\
CoulombEwaldRealForce.h\
CoulombEwaldRealTableForceBase.h\
CoulombEwaldRealTableForce.h\
CoulombForce.h\
CoulombForceDiElec.h\
CoulombMultiGridDirectForceBase.h\
CoulombMultiGridDirectForce.h\
CoulombMultiGridDirectTableForceBase.h\
CoulombMultiGridDirectTableForce.h\
CoulombSCPISMForce.h\
CoulombTableForceBase.h\
CoulombTableForce.h\
DihedralSystemForce.h\
DihedralSystemForceBase.h\
ElectricFieldSystemForce.h\
ElectricFieldSystemForceBase.h\
ExtendedCompareForce.h\
ExtendedForce.h\
ExtendedTimeForce.h\
ExternalGravitationSystemForce.h\
ExternalMagneticFieldExtendedForce.h\
Force.h\
ForceGroup.h\
FrictionExtendedForce.h\
GravitationForce.h\
HarmDihedralSystemForce.h\
HarmDihedralSystemForceBase.h\
ImproperSystemForce.h\
ImproperSystemForceBase.h\
LennardJonesForce.h\
LennardJonesTableForce.h\
LennardJonesTableForceBase.h\
MTorsionSystemForce.h\
MagneticDipoleForce.h\
MagneticDipoleMirrorSystemForce.h\
MagneticDipoleMirrorSystemForceBase.h\
MetaForce.h\
MollyForce.h\
NonbondedCutoffBornForce.h\
NonbondedCutoffForce.h\
NonbondedCutoffForceBase.h\
NonbondedCutoffMollyForce.h\
NonbondedCutoffSystemForce.h\
NonbondedFullBornForce.h\
NonbondedFullEwaldSystemForce.h\
NonbondedFullEwaldSystemForceBase.h\
NonbondedFullSystemForce.h\
NonbondedFullSystemForceBase.h\
NonbondedMultiGridSystemForce.h\
NonbondedMultiGridSystemForceBase.h\
NonbondedPMEwaldSystemForce.h\
NonbondedPMEwaldSystemForceBase.h\
NonbondedSimpleFullBornForce.h\
NonbondedSimpleFullSystemForce.h\
NonbondedSimpleFullSystemForceBase.h\
OneAtomPair.h\
OneAtomPairFull.h\
OneAtomPairTwo.h\
OneAtomPairTwoFull.h\
PaulTrapExtendedForce.h\
PaulTrapExtendedForceBase.h\
ReducedHessCoulomb.h\
ReducedHessCoulombBornRadii.h\
ReducedHessCoulombDiElec.h\
ReducedHessCoulombSCPISM.h\
ReducedHessLennardJones.h\
SphericalRestraintSystemForce.h\
SphericalSystemForce.h\
SystemCompareForce.h\
SystemForce.h\
SystemTimeForce.h\
TimeForce.h\
WrapperMetaForce.h\
evaluateBorn.h\
oneAtomContraints.h\
HessDihedral.h


# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libforces.a

libforces_a_SOURCES = $(P_FORCES_SRC) $(P_FORCES_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES)# -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_FORCES_SRC)

include $(top_srcdir)/Makefile.depend
