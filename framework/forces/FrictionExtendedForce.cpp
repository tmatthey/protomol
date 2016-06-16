#include "FrictionExtendedForce.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"
#include "Parallel.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ FrictionExtendedForce


	const string FrictionExtendedForce::keyword("Friction");

	FrictionExtendedForce::FrictionExtendedForce(): ExtendedForce(), myF(0.0), myRnd(Vector3D(0.0, 0.0, 0.0))
	{
	}


	FrictionExtendedForce::FrictionExtendedForce(Real f, const Vector3D& random): ExtendedForce(), myF(f), myRnd(random)
	{
	}


	void FrictionExtendedForce::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-k", Value(myF), 0.0));
		parameters.push_back(Parameter("-rnd", Value(myRnd), Vector3D(0.0, 0.0, 0.0)));
	}


	Force* FrictionExtendedForce::doMake(string&, vector<Value> values) const
	{
		return new FrictionExtendedForce(values[0], values[1]);
	}

	void FrictionExtendedForce::evaluate(const GenericTopology* topo,
	                                     const Vector3DBlock* positions,
	                                     const Vector3DBlock* velocities,
	                                     Vector3DBlock* forces,
	                                     ScalarStructure* energies)
	{
		doEvaluate(topo, positions, velocities, forces, energies, 0, positions->size());
	}


	void FrictionExtendedForce::parallelEvaluate(const GenericTopology* topo,
	                                             const Vector3DBlock* positions,
	                                             const Vector3DBlock* velocities,
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
				doEvaluate(topo, positions, velocities, forces, energies, from, to);
			}
		}
	}


	void FrictionExtendedForce::doEvaluate(const GenericTopology*,
	                                       const Vector3DBlock* /*positions*/,
	                                       const Vector3DBlock* velocities,
	                                       Vector3DBlock* forces,
	                                       ScalarStructure*,
	                                       int from, int to)
	{
		for (int i = from; i < to; i++)
		{
			(*forces)[i].x += (*velocities)[i].x * myF + (randomNumber() - .5) * myRnd.x;
			(*forces)[i].y += (*velocities)[i].y * myF + (randomNumber() - .5) * myRnd.y;
			(*forces)[i].z += (*velocities)[i].z * myF + (randomNumber() - .5) * myRnd.z;
		}
	}


	unsigned int FrictionExtendedForce::numberOfBlocks(const GenericTopology*,
	                                                   const Vector3DBlock* pos)
	{
		return std::min(static_cast<unsigned int>(Parallel::getAvailableNum()), static_cast<unsigned int>(pos->size()));
	}
}
