# base

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libOSGfrontend_a_AR = $(AR2)

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
P_OSGINTEGRATORS = $(top_srcdir)/applications/oSGProtomol-app/integrators
P_OSGPROTOMOL    = $(top_srcdir)/applications/oSGProtomol-app

# base Sources
P_OSGFRONTEND_SRC =\
OutputOSGProperties.cpp\
OutputFinalXSC.cpp\
OutputFinalPSF.cpp\
buildOSGTopology.cpp\
OutputOSGDCDTrajectory.cpp\
OutputOSGXYZTrajectoryPos.cpp

# base Includes
P_OSGFRONTEND_H =\
OutputOSGProperties.h\
OutputFinalXSC.h\
OutputFinalPSF.h\
buildOSGTopology.h\
OutputOSGDCDTrajectory.h\
OutputOSGXYZTrajectoryPos.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libOSGfrontend.a

libOSGfrontend_a_SOURCES = $(P_OSGFRONTEND_SRC) $(P_OSGFRONTEND_H)

INCLUDES = @FFT_CFLAGS@ -I$(P_OSGFRONTEND) -I$(P_OSGBASE) -I$(P_OSGIO) -I$(P_OSGINTEGRATORS) -I$(P_BASE) -I$(P_IO) -I$(P_FRONTEND) -I$(P_TOPOLOGY) -I$(P_INTEGRATORS) -I$(P_FORCES) -I$(P_OSGPROTOMOL)

SOURCEFILES = $(P_OSGFRONTEND_SRC)

include $(top_srcdir)/Makefile.depend
