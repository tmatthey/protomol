# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libISGforces_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories
P_ISGBASE     = $(top_srcdir)/applications/iSGProtomol-app/base
P_ISGIO       = $(top_srcdir)/applications/iSGProtomol-app/io
P_ISGFORCES   = $(top_srcdir)/applications/iSGProtomol-app/forces
P_ISGFRONTEND = $(top_srcdir)/applications/iSGProtomol-app/frontend
P_ISGTOPOLOGY = $(top_srcdir)/applications/iSGProtomol-app/topology

# base Sources
P_ISGFORCES_SRC =\
iSGAngleSystemForceBase.cpp\
iSGBondSystemForceBase.cpp\
iSGDihedralSystemForceBase.cpp\
iSGImproperSystemForceBase.cpp\
iSGLennardJonesForce.cpp\
iSGCoulombForce.cpp\
iSGNonbondedFullEwaldSystemForceBase.cpp\
iSGNonbondedPMEwaldSystemForceBase.cpp\
iSGIdealGasCoulombForce.cpp\
iSGIdealGasLennardJonesForce.cpp\
iSGLennardJonesTableForce.cpp\
iSGCoulombEwaldRealForce.cpp\
iSGCoulombEwaldRealTableForce.cpp

# base Includes
P_ISGFORCES_H =\
iSGAngleSystemForceBase.h\
iSGAngleSystemForce.h\
iSGBondSystemForceBase.h\
iSGBondSystemForce.h\
iSGDihedralSystemForceBase.h\
iSGDihedralSystemForce.h\
iSGImproperSystemForceBase.h\
iSGImproperSystemForce.h\
iSGMTorsionSystemForce.h\
iSGLennardJonesForce.h\
iSGOneAtomPairFull.h\
iSGOneAtomPair.h\
iSGCoulombForce.h\
iSGNonbondedFullEwaldSystemForceBase.h\
iSGNonbondedFullEwaldSystemForce.h\
iSGNonbondedPMEwaldSystemForceBase.h\
iSGNonbondedPMEwaldSystemForce.h\
iSGIdealGasCoulombForce.h\
iSGIdealGasLennardJonesForce.h\
iSGOneAtomPairFullTwo.h\
iSGOneAtomPairTwo.h\
iSGLennardJonesTableForce.h\
iSGCoulombEwaldRealForce.h\
iSGCoulombEwaldRealTableForce.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libISGforces.a

libISGforces_a_SOURCES = $(P_ISGFORCES_SRC) $(P_ISGFORCES_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_ISGFORCES) -I$(P_BASE) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_ISGTOPOLOGY) -I$(P_ISGBASE) -I$(P_PARALLEL) -I$(P_FORCES) -I$(P_IO)

SOURCEFILES = $(P_ISGFORCES_SRC)

include $(top_srcdir)/Makefile.depend
