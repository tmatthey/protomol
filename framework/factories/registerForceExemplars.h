/* -*- c++ -*- */
#ifndef REGISTERFORCEEXEMPLARS_H
#define REGISTERFORCEEXEMPLARS_H


namespace ProtoMol
{
	class GenericTopology;

	/// registers all force prototypes by calling all registerForceExemplars* to be recognized by the parser
	void registerForceExemplars(const GenericTopology* topo);
}
#endif /* REGISTERFORCEEXEMPLARS_H */
