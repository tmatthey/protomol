# imd

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libimd_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# imd Sources
P_IMD_SRC = IMDElf.cpp imd.cpp vmdsock.cpp HapticSystemForce.cpp

# imd Includes
P_IMD_H = IMDElf.h imd.h vmdsock.h HapticSystemForce.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libimd.a

libimd_a_SOURCES = $(P_IMD_SRC) $(P_IMD_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD)# -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_IMD_SRC)

include $(top_srcdir)/Makefile.depend
