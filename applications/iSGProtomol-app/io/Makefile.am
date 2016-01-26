# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libISGio_a_AR = $(AR2)

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
P_ISGIO_SRC =\
iSGPARReader.cpp\
TRANSReader.cpp\
XSCReader.cpp\
XSCWriter.cpp

# base Includes
P_ISGIO_H =\
iSGPARReader.h\
TRANSReader.h\
XSCReader.h\
XSCWriter.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libISGio.a

libISGio_a_SOURCES = $(P_ISGIO_SRC) $(P_ISGIO_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_ISGIO) -I$(P_ISGBASE) -I$(P_IO) -I$(P_BASE) -I$(P_ISGTOPOLOGY) #-I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(=_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_ISGIO_SRC)

include $(top_srcdir)/Makefile.depend
