# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libISGintegrators_a_AR = $(AR2)

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
P_ISGFACTORIES = $(top_srcdir)/applications/iSGProtomol-app/factories

# base Sources
P_ISGINTEGRATORS_SRC =\
iSGIntegrator.cpp

# base Includes
P_ISGINTEGRATORS_H =\
iSGIntegrator.h\
iSGModifyForces.h
#iSGModifierPreForceChemostat.h\
#iSGModifierPostForceChemostat.h
#iSGCheckForCompletion.h\
#iSGUpdateTopology.h\
#iSGPickNewMolecule.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libISGintegrators.a

libISGintegrators_a_SOURCES = $(P_ISGINTEGRATORS_SRC) $(P_ISGINTEGRATORS_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_ISGFRONTEND) -I$(P_BASE) -I$(P_FRONTEND) -I$(P_ISGBASE) -I$(P_TOPOLOGY) -I$(P_FRONTEND) -I$(P_INTEGRATORS) -I$(P_FORCES) -I$(P_PARALLEL) -I$(P_ISGIO) -I$(P_IO)

SOURCEFILES = $(P_ISGINTEGRATORS_SRC)

include $(top_srcdir)/Makefile.depend
