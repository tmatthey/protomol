/*  -*- c++ -*-  */
#ifndef PARSECOMMANDLINE_H
#define PARSECOMMANDLINE_H

#include <vector>
#include <string>

namespace ProtoMol
{
	class Configuration;
	class GenericTopology;
	//________________________________________________________ parseCommandLine
	std::vector<std::vector<std::string>> parseCommandLine(int argc, char** argv,
	                                                       const Configuration* config = NULL,
	                                                       void (*registerForceExemplarsFunction)(const GenericTopology*) = NULL);

	// Parses the command line arguments and returns 2d vector string.
	// The parsing expects keyword with Constant::ParseCommandLine::prefix prefix.
	// NB: We pass the function (registerForceExemplarsFunction) to
	// decouple registerForceExemplars from parseCommandLine, thus
	// we do not need to link registerForceExemplars if called with NULL
}
#endif
