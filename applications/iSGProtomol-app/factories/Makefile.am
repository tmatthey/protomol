# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libISGfactories_a_AR = $(AR2)

P_BASE           = $(top_srcdir)/framework/base
P_IO             = $(top_srcdir)/framework/io
P_PARALLEL       = $(top_srcdir)/framework/parallel
P_TOPOLOGY       = $(top_srcdir)/framework/topology
P_FORCES         = $(top_srcdir)/framework/forces
P_IMD            = $(top_srcdir)/framework/imd
P_INTEGRATORS    = $(top_srcdir)/framework/integrators
P_FRONTEND       = $(top_srcdir)/framework/frontend
P_FACTORIES      = $(top_srcdir)/framework/factories
P_ISGBASE        = $(top_srcdir)/applications/iSGProtomol-app/base
P_ISGIO          = $(top_srcdir)/applications/iSGProtomol-app/io
P_ISGFORCES      = $(top_srcdir)/applications/iSGProtomol-app/forces
P_ISGFRONTEND    = $(top_srcdir)/applications/iSGProtomol-app/frontend
P_ISGTOPOLOGY    = $(top_srcdir)/applications/iSGProtomol-app/topology
P_ISGFACTORIES   = $(top_srcdir)/applications/iSGProtomol-app/factories
P_ISGINTEGRATORS = $(top_srcdir)/applications/iSGProtomol-app/integrators
P_ISGPROTOMOL    = $(top_srcdir)/applications/iSGProtomol-app

# base Sources
P_ISGFACTORIES_SRC =\
iSGregisterOutputExemplars.cpp\
iSGregisterIntegratorExemplars.cpp\
iSGregisterForceExemplars.cpp\
iSGregisterForceExemplarsBonded.cpp\
iSGregisterForceExemplarsCutoff.cpp\
iSGregisterForceExemplarsFull.cpp\
iSGregisterForceExemplarsSimpleFull.cpp\
iSGregisterForceExemplarsFastElectrostatic.cpp\
iSGregisterForceExemplarsIdealGas.cpp

# base Includes
P_ISGFACTORIES_H =\
iSGregisterOutputExemplars.h\
iSGregisterForceExemplars.h\
iSGregisterIntegratorExemplars.h\
iSGregisterForceExemplars.h\
iSGregisterForceExemplarsBonded.h\
iSGregisterForceExemplarsCutoff.h\
iSGregisterForceExemplarsFull.h\
iSGregisterForceExemplarsSimpleFull.h\
iSGregisterForceExemplarsFastElectrostatic.h\
iSGregisterForceExemplarsIdealGas.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libISGfactories.a

libISGfactories_a_SOURCES = $(P_ISGFACTORIES_SRC) $(P_ISGFACTORIES_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_ISGIO) -I$(P_ISGBASE) -I$(P_IO) -I$(P_BASE) -I$(P_ISGTOPOLOGY) -I$(P_ISGFORCES) -I$(P_TOPOLOGY) -I$(P_ISGFRONTEND) -I$(P_FRONTEND) -I$(P_INTEGRATORS) -I$(P_ISGPROTOMOL) -I$(P_ISGINTEGRATORS) -I$(P_ISGFACTORIES) -I$(P_PARALLEL) -I$(P_FORCES) -I$(P_FACTORIES)

SOURCEFILES = $(P_ISGFACTORIES_SRC)

include $(top_srcdir)/Makefile.depend
