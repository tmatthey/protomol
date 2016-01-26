# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

SUBDIRS = base io factories forces frontend integrators

bin_PROGRAMS = oSGProtomol

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
P_OSGFRONTEND    = $(top_srcdir)/applications/oSGProtomol-app/frontend
P_OSGFACTORIES   = $(top_srcdir)/applications/oSGProtomol-app/factories
P_OSGINTEGRATORS = $(top_srcdir)/applications/oSGProtomol-app/integrators


## oSGProtomol
oSGProtomol_SOURCES = oSGProtomol.cpp

oSGProtomol_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a $(P_PARALLEL)/libparallel.a $(P_TOPOLOGY)/libtopology.a $(P_FORCES)/libforces.a $(P_IMD)/libimd.a $(P_INTEGRATORS)/libintegrators.a $(P_FRONTEND)/libfrontend.a $(P_FACTORIES)/registerTopologyExemplars.o $(P_OSGBASE)/libOSGbase.a $(P_OSGIO)/libOSGio.a $(P_OSGINTEGRATORS)/libOSGintegrators.a
INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)  -I$(P_OSGBASE) -I$(P_OSGIO) -I$(P_OSGFORCES) -I$(P_OSGFRONTEND) -I$(P_OSGFACTORIES) -I$(P_OSGINTEGRATORS)

oSGProtomol_LDADD =-L$(P_BASE) -L$(P_IO) -L$(P_PARALLEL) -L$(P_TOPOLOGY) -L$(P_FORCES) -L$(P_IMD) -L$(P_INTEGRATORS) -L$(P_FRONTEND) -L$(P_OSGBASE) -L$(P_OSGIO) -L$(P_OSGFORCES) -L$(P_OSGFRONTEND) -L$(P_OSGFACTORIES) -L$(P_OSGINTEGRATORS) -lOSGintegrators -lOSGfactories -lOSGfrontend -lOSGforces -lOSGio -lOSGbase $(P_FACTORIES)/registerTopologyExemplars.o -lfrontend -lintegrators -limd -lforces -ltopology -lparallel -lio -lbase @SOCKET_LIBS@ @FFT_LIBS@ @PNG_LIBS@
SOURCEFILES		= $(oSGProtomol_SOURCES)


include $(top_srcdir)/Makefile.depend
