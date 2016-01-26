# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

bin_PROGRAMS = coulomb

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

P_COULOMB_SRC = \
coulomb.cpp\
registerForceExemplars.cpp\
registerIntegratorExemplars.cpp\
registerOutputExemplars.cpp\
registerTopologyExemplars.cpp

P_COULOMB_H =\
registerForceExemplars.h\
registerIntegratorExemplars.h\
registerOutputExemplars.h\
registerTopologyExemplars.h

## Coulomb
coulomb_SOURCES = $(P_COULOMB_SRC) $(P_COULOMB_H)
coulomb_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a $(P_PARALLEL)/libparallel.a $(P_TOPOLOGY)/libtopology.a $(P_FORCES)/libforces.a $(P_INTEGRATORS)/libintegrators.a $(P_FRONTEND)/libfrontend.a 
INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_INTEGRATORS) -I$(P_FRONTEND)

coulomb_LDADD =-L$(P_BASE) -L$(P_IO) -L$(P_PARALLEL) -L$(P_TOPOLOGY) -L$(P_FORCES) -L$(P_INTEGRATORS) -L$(P_FRONTEND) -lfrontend -lintegrators -lforces -ltopology -lparallel -lio -lbase @SOCKET_LIBS@ @FFT_LIBS@
SOURCEFILES		= $(coulomb_SOURCES)


include $(top_srcdir)/Makefile.depend
