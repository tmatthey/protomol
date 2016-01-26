# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libOSGforces_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories
P_OSGBASE     = $(top_srcdir)/applications/oSGProtomol-app/base
P_OSGIO       = $(top_srcdir)/applications/oSGProtomol-app/io
P_OSGFORCES   = $(top_srcdir)/applications/oSGProtomol-app/forces
P_OSGFRONTEND = $(top_srcdir)/applications/oSGProtomol-app/frontend
P_OSGTOPOLOGY = $(top_srcdir)/applications/oSGProtomol-app/topology

# base Sources
P_OSGFORCES_SRC =\
oSGLennardJonesForce.cpp
#iSGCoulombForce.cpp\
#iSGNonbondedFullEwaldSystemForceBase.cpp\
#iSGNonbondedPMEwaldSystemForceBase.cpp\
#iSGCoulombEwaldRealForce.cpp\
#iSGCoulombEwaldRealTableForce.cpp

# base Includes
P_OSGFORCES_H =\
oSGOneAtomPairFull.h\
oSGOneAtomPair.h\
oSGLennardJonesForce.h\
oSGOneAtomPairFullTwo.h\
oSGOneAtomPairTwo.h
#iSGCoulombForce.h\
#iSGNonbondedFullEwaldSystemForceBase.h\
#iSGNonbondedFullEwaldSystemForce.h\
#iSGNonbondedPMEwaldSystemForceBase.h\
#iSGNonbondedPMEwaldSystemForce.h\
#iSGCoulombEwaldRealForce.h\
#iSGCoulombEwaldRealTableForce.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libOSGforces.a

libOSGforces_a_SOURCES = $(P_OSGFORCES_SRC) $(P_OSGFORCES_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_OSGFORCES) -I$(P_BASE) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_OSGTOPOLOGY) -I$(P_OSGBASE) -I$(P_PARALLEL) -I$(P_FORCES) -I$(P_IO)

SOURCEFILES = $(P_OSGFORCES_SRC)

include $(top_srcdir)/Makefile.depend
