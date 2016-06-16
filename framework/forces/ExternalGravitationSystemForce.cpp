#include "ExternalGravitationSystemForce.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"
#include "Parallel.h"
#include "pmconstants.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ ExternalGravitationSystemForce


	static const Real SI_FACTOR = 1e-3 * Constant::SI::KCAL / Constant::SI::LENGTH_AA;


	const string ExternalGravitationSystemForce::keyword("ExternalGravitation");

	ExternalGravitationSystemForce::ExternalGravitationSystemForce(): SystemForce(), myG(Vector3D(0.0, 0.0, 0.0))
	{
	}

	ExternalGravitationSystemForce::ExternalGravitationSystemForce(const Vector3D& g): SystemForce(), myG(g)
	{
	}


	void ExternalGravitationSystemForce::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-g", Value(myG), Vector3D(0.0, 0.0, -9.81), Text("external parallel gravitation [m/s^2]")));
	}


	Force* ExternalGravitationSystemForce::doMake(string&, vector<Value> values) const
	{
		return new ExternalGravitationSystemForce(values[0]);
	}

	void ExternalGravitationSystemForce::evaluate(const GenericTopology* topo,
	                                              const Vector3DBlock* positions,
	                                              Vector3DBlock* forces,
	                                              ScalarStructure* energies)
	{
		doEvaluate(topo, positions, forces, energies, 0, positions->size());
	}


	void ExternalGravitationSystemForce::parallelEvaluate(const GenericTopology* topo,
	                                                      const Vector3DBlock* positions,
	                                                      Vector3DBlock* forces,
	                                                      ScalarStructure* energies)
	{
		unsigned int n = positions->size();
		unsigned int count = numberOfBlocks(topo, positions);

		for (unsigned int i = 0; i < count; i++)
		{
			if (Parallel::next())
			{
				int to = (n * (i + 1)) / count;
				if (to > (int)n)
					to = n;
				int from = (n * i) / count;
				doEvaluate(topo, positions, forces, energies, from, to);
			}
		}
	}


	void ExternalGravitationSystemForce::doEvaluate(const GenericTopology* topo,
	                                                const Vector3DBlock* /*positions*/,
	                                                Vector3DBlock* forces,
	                                                ScalarStructure*,
	                                                int from, int to)
	{
		Vector3D g(myG * SI_FACTOR);
		for (int i = from; i < to; i++)
		{
			(*forces)[i] += g * topo->atoms[i].scaledMass;
		}
	}


	unsigned int ExternalGravitationSystemForce::numberOfBlocks(const GenericTopology*,
	                                                            const Vector3DBlock* pos)
	{
		return std::min(static_cast<unsigned int>(Parallel::getAvailableNum()), static_cast<unsigned int>(pos->size()));
	}
}
