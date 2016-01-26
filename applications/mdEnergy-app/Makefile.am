# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

bin_PROGRAMS = mdEnergy

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

P_MDENERGY_SRC = \
mdEnergy.cpp

P_MDENERGY_H = 

## mdEnergy
mdEnergy_SOURCES = $(P_MDENERGY_SRC) $(P_MDENERGY_H)
mdEnergy_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a $(P_PARALLEL)/libparallel.a $(P_TOPOLOGY)/libtopology.a $(P_FORCES)/libforces.a $(P_IMD)/libimd.a $(P_INTEGRATORS)/libintegrators.a $(P_FRONTEND)/libfrontend.a $(P_FACTORIES)/libfactories.a
INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

mdEnergy_LDADD = -L$(P_BASE) -L$(P_IO) -L$(P_PARALLEL) -L$(P_TOPOLOGY) -L$(P_FORCES) -L$(P_IMD) -L$(P_INTEGRATORS) -L$(P_FRONTEND) -L$(P_FACTORIES) -lfactories -lfrontend -lintegrators -limd -lforces -ltopology -lparallel -lio -lbase @SOCKET_LIBS@ @FFT_LIBS@ @PNG_LIBS@
SOURCEFILES		= $(mdEnergy_SOURCES)


include $(top_srcdir)/Makefile.depend
