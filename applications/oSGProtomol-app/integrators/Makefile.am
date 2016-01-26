# integrators

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libOSGintegrators_a_AR = $(AR2)

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
P_OSGFRONTEND = $(top_srcdir)/applications/osGProtomol-app/frontend
P_OSGFACTORIES = $(top_srcdir)/applications/oSGProtomol-app/factories

# integrators Sources
P_OSGINTEGRATORS_SRC =\
oSGIntegrator.cpp\
muVTIntegrator.cpp\
NfPTIntegrator.cpp

# integrators Includes
P_OSGINTEGRATORS_H =\
oSGModifierPreForceChemostat.h\
oSGModifierPostForceChemostat.h\
oSGIntegrator.h\
muVTIntegrator.h\
NfPTIntegrator.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libOSGintegrators.a

libOSGintegrators_a_SOURCES = $(P_OSGINTEGRATORS_SRC) $(P_OSGINTEGRATORS_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_OSGFRONTEND) -I$(P_BASE) -I$(P_FRONTEND) -I$(P_OSGBASE) -I$(P_TOPOLOGY) -I$(P_FRONTEND) -I$(P_INTEGRATORS) -I$(P_FORCES) -I$(P_PARALLEL) -I$(P_OSGIO) -I$(P_IO)

SOURCEFILES = $(P_OSGINTEGRATORS_SRC)

include $(top_srcdir)/Makefile.depend
