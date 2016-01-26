# factories

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libOSGfactories_a_AR = $(AR2)

P_BASE           = $(top_srcdir)/framework/base
P_IO             = $(top_srcdir)/framework/io
P_PARALLEL       = $(top_srcdir)/framework/parallel
P_TOPOLOGY       = $(top_srcdir)/framework/topology
P_FORCES         = $(top_srcdir)/framework/forces
P_IMD            = $(top_srcdir)/framework/imd
P_INTEGRATORS    = $(top_srcdir)/framework/integrators
P_FRONTEND       = $(top_srcdir)/framework/frontend
P_FACTORIES      = $(top_srcdir)/framework/factories
P_OSGBASE        = $(top_srcdir)/applications/oSGProtomol-app/base
P_OSGIO          = $(top_srcdir)/applications/oSGProtomol-app/io
P_OSGFORCES      = $(top_srcdir)/applications/oSGProtomol-app/forces
P_OSGFACTORIES   = $(top_srcdir)/applications/oSGProtomol-app/factories
P_OSGFRONTEND    = $(top_srcdir)/applications/oSGProtomol-app/frontend
P_OSGINTEGRATORS = $(top_srcdir)/applications/oSGProtomol-app/integrators
P_OSGPROTOMOL    = $(top_srcdir)/applications/oSGProtomol-app

# factories Sources
P_OSGFACTORIES_SRC =\
oSGregisterOutputExemplars.cpp\
oSGregisterIntegratorExemplars.cpp\
oSGregisterForceExemplars.cpp\
registerForceExemplarsBonded.cpp\
oSGregisterForceExemplarsCutoff.cpp\
oSGregisterForceExemplarsFull.cpp\
oSGregisterForceExemplarsSimpleFull.cpp

# factories Includes
P_OSGFACTORIES_H =\
oSGregisterOutputExemplars.h\
oSGregisterIntegratorExemplars.h\
oSGregisterForceExemplars.h\
registerForceExemplarsBonded.h\
oSGregisterForceExemplarsCutoff.h\
oSGregisterForceExemplarsFull.h\
oSGregisterForceExemplarsSimpleFull.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libOSGfactories.a

libOSGfactories_a_SOURCES = $(P_OSGFACTORIES_SRC) $(P_OSGFACTORIES_H)

INCLUDES = -I$(P_OSGIO) -I$(P_OSGBASE) -I$(P_IO) -I$(P_BASE) -I$(P_OSGFORCES) -I$(P_TOPOLOGY) -I$(P_FRONTEND) -I$(P_OSGFRONTEND) -I$(P_INTEGRATORS) -I$(P_OSGPROTOMOL) -I$(P_OSGINTEGRATORS) -I$(P_OSGFACTORIES) -I$(P_PARALLEL) -I$(P_FORCES) -I$(P_FACTORIES) @FFT_CFLAGS@

SOURCEFILES = $(P_OSGFACTORIES_SRC)

include $(top_srcdir)/Makefile.depend
