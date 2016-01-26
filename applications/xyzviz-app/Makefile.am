# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

bin_PROGRAMS = xyzviz

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

P_XYZVIZ_SRC = \
openglutilities.cpp\
xyzviz.cpp

P_XYZVIZ_H = \
openglutilities.h

## xyzviz
xyzviz_SOURCES = $(P_XYZVIZ_SRC) $(P_XYZVIZ_H)
xyzviz_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a 
INCLUDES =-I$(P_BASE) -I$(P_IO) @GLUT_CFLAGS@

xyzviz_LDADD =-L$(P_BASE) -L$(P_IO) @SOCKET_LIBS@ @FFT_LIBS@ @GLUT_LIBS@ -lio -lbase 
SOURCEFILES		= $(xyzviz_SOURCES)


include $(top_srcdir)/Makefile.depend
