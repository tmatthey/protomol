# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

bin_PROGRAMS = eigset

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

P_EIGSET_SRC = \
eigset.cpp

P_EIGSET_H =

## eigset
eigset_SOURCES = $(P_EIGSET_SRC) $(P_EIGSET_H)
eigset_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a 
INCLUDES =-I$(P_BASE) -I$(P_IO)

eigset_LDADD =-L$(P_BASE) -L$(P_IO) -lio -lbase @SOCKET_LIBS@ @FFT_LIBS@
SOURCEFILES		= $(eigset_SOURCES)


include $(top_srcdir)/Makefile.depend
