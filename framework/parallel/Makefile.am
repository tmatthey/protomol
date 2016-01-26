# parallel

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libparallel_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# parallel Sources
P_PARALLEL_SRC =\
Parallel.cpp\
ParallelType.cpp

# parallel Includes
P_PARALLEL_H =\
Parallel.h\
ParallelType.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libparallel.a

libparallel_a_SOURCES = $(P_PARALLEL_SRC) $(P_PARALLEL_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL)#-I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES)

SOURCEFILES = $(P_PARALLEL_SRC)

include $(top_srcdir)/Makefile.depend
