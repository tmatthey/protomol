# io

AUTOMAKE_OPTIONS	= foreign dist-tarZ no-dependencies

#Change AR Flag
AR = $(AR1)
libio_a_AR = $(AR2) 

P_BASE        = $(top_srcdir)/framework/base
P_IO          = $(top_srcdir)/framework/io
P_PARALLEL    = $(top_srcdir)/framework/parallel
P_TOPOLOGY    = $(top_srcdir)/framework/topology
P_FORCES      = $(top_srcdir)/framework/forces
P_IMD         = $(top_srcdir)/framework/imd
P_INTEGRATORS = $(top_srcdir)/framework/integrators
P_FRONTEND    = $(top_srcdir)/framework/frontend
P_FACTORIES   = $(top_srcdir)/framework/factories

# io Sources
P_IO_SRC =\
CheckpointInputStream.cpp\
CheckpointOutputStream.cpp\
ConfigurationReader.cpp\
DCDTrajectoryReader.cpp\
DCDTrajectoryWriter.cpp\
File.cpp\
AMBERReader.cpp\
EigenvectorReader.cpp\
EigenvectorTextReader.cpp\
PARReader.cpp\
PARWriter.cpp\
PDBReader.cpp\
CRDReader.cpp\
PDBWriter.cpp\
PNGWriter.cpp\
PPMReader.cpp\
PPMWriter.cpp\
PGMReader.cpp\
PGMWriter.cpp\
PSFReader.cpp\
PSFWriter.cpp\
Reader.cpp\
REMConfigurationReader.cpp\
TempFileReader.cpp\
Writer.cpp\
XFigWriter.cpp\
XYZBinReader.cpp\
XYZBinWriter.cpp\
XYZBinRevWriter.cpp\
XYZReader.cpp\
XYZTrajectoryReader.cpp\
XYZTrajectoryWriter.cpp\
XYZWriter.cpp

# io Includes
P_IO_H =\
CheckpointInputStream.h\
CheckpointOutputStream.h\
ConfigurationReader.h\
DCDTrajectoryReader.h\
DCDTrajectoryWriter.h\
EigenvectorReader.h\
EigenvectorTextReader.h\
File.h\
AMBERReader.h\
PARReader.h\
PARWriter.h\
PDBReader.h\
CRDReader.h\
PDBWriter.h\
PNGWriter.h\
PPMReader.h\
PPMWriter.h\
PGMReader.h\
PGMWriter.h\
PSFReader.h\
PSFWriter.h\
Reader.h\
REMConfigurationReader.h\
TempFileReader.h\
Writer.h\
XFigWriter.h\
XYZBinReader.h\
XYZBinRevWriter.h\
XYZBinWriter.h\
XYZReader.h\
XYZTrajectoryReader.h\
XYZTrajectoryWriter.h\
XYZWriter.h

# Temporary library to get subdirectory source files, blech!
noinst_LIBRARIES = libio.a

libio_a_SOURCES = $(P_IO_SRC) $(P_IO_H)

INCLUDES =-I$(P_BASE) -I$(P_IO) -I$(P_PARALLEL) -I$(P_TOPOLOGY) -I$(P_FORCES) -I$(P_IMD) -I$(P_INTEGRATORS) -I$(P_FRONTEND) -I$(P_FACTORIES) @BOOST_CFLAGS@

SOURCEFILES = $(P_IO_SRC)

include $(top_srcdir)/Makefile.depend
