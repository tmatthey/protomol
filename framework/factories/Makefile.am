# factories

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libfactories_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# factories Sources
P_FACTORIES_SRC =\
registerForceExemplars.cpp\
registerForceExemplarsBonded.cpp\
registerForceExemplarsCutoff.cpp\
registerForceExemplarsFastElectrostatic.cpp\
registerForceExemplarsFull.cpp\
registerForceExemplarsOther.cpp\
registerForceExemplarsSimpleFull.cpp\
registerIntegratorExemplars.cpp\
registerOutputExemplars.cpp\
registerTopologyExemplars.cpp

# factories Includes
P_FACTORIES_H =\
registerForceExemplars.h\
registerForceExemplarsBonded.h\
registerForceExemplarsCutoff.h\
registerForceExemplarsFastElectrostatic.h\
registerForceExemplarsFull.h\
registerForceExemplarsOther.h\
registerForceExemplarsSimpleFull.h\
registerIntegratorExemplars.h\
registerOutputExemplars.h\
registerTopologyExemplars.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libfactories.a

libfactories_a_SOURCES = $(P_FACTORIES_SRC) $(P_FACTORIES_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_FACTORIES_SRC)

include $(top_srcdir)/Makefile.depend
