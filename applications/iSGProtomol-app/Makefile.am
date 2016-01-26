# -*- makefile -*-
#
# See the LICENSE file in the top level directory for license
# information

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

SUBDIRS =  base io forces frontend factories integrators

bin_PROGRAMS = iSGProtomol

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
P_ISGFACTORIES   = $(top_srcdir)/applications/iSGProtomol-app/factories
P_ISGINTEGRATORS = $(top_srcdir)/applications/iSGProtomol-app/integrators


## ISGProtomol
iSGProtomol_SOURCES = iSGProtomol.cpp
iSGProtomol_DEPENDENCIES = $(P_BASE)/libbase.a $(P_IO)/libio.a $(P_PARALLEL)/libparallel.a $(P_TOPOLOGY)/libtopology.a $(P_FORCES)/libforces.a $(P_IMD)/libimd.a $(P_INTEGRATORS)/libintegrators.a $(P_FRONTEND)/libfrontend.a $(P_FACTORIES)/registerTopologyExemplars.o $(P_ISGBASE)/libISGbase.a $(P_ISGIO)/libISGio.a $(P_ISGINTEGRATORS)/libISGintegrators.a
INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES) -I$(P_ISGBASE) -I$(P_ISGIO) -I$(P_ISGFORCES) -I$(P_ISGFRONTEND) -I$(P_ISGFACTORIES) -I$(P_ISGINTEGRATORS)

iSGProtomol_LDADD =-L$(P_BASE) -L$(P_IO) -L$(P_PARALLEL) -L$(P_TOPOLOGY) -L$(P_FORCES) -L$(P_IMD) -L$(P_INTEGRATORS) -L$(P_FRONTEND) -L$(P_ISGBASE) -L$(P_ISGIO) -L$(P_ISGFORCES) -L$(P_ISGFRONTEND) -L$(P_ISGFACTORIES) -L$(P_ISGINTEGRATORS) -lISGintegrators -lISGfactories -lISGfrontend -lISGforces -lISGio -lISGbase $(P_FACTORIES)/registerTopologyExemplars.o -lfrontend -lintegrators -limd -lforces -ltopology -lparallel -lio -lbase @SOCKET_LIBS@ @FFT_LIBS@ @PNG_LIBS@
SOURCEFILES		= $(iSGProtomol_SOURCES)


include $(top_srcdir)/Makefile.depend
