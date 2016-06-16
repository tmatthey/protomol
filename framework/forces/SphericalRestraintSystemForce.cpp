#include "SphericalRestraintSystemForce.h"
#include "GenericTopology.h"
#include "Vector3DBlock.h"
#include "Parallel.h"
#include "ScalarStructure.h"
#include "mathutilities.h"
#include "pmconstants.h"
using std::vector;
using std::string;
using std::cout;
using std::endl;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ SphericalRestraintSystemForce
	const string SphericalRestraintSystemForce::keyword("SphericalRestraint");


	SphericalRestraintSystemForce::SphericalRestraintSystemForce(): SystemForce(),
	                                                                //myCenter(Vector3D(0.0,0.0,0.0)),
	                                                                myK(0.0)
	{
	}

	SphericalRestraintSystemForce::SphericalRestraintSystemForce(Real k): SystemForce(),
	                                                                      //      myCenter(center),
	                                                                      myK(k)
	{
	}

	void SphericalRestraintSystemForce::doEvaluate(const GenericTopology* topo,
	                                               const Vector3DBlock* positions,
	                                               Vector3DBlock* forces,
	                                               ScalarStructure* energies, int from, int to)
	{
		Real e = 0.0;
		Real c = 2 * myK;
		for (int i = from; i < to; i++)
		{
			Vector3D rc = (*positions)[i];
			Real rc2 = rc.normSquared();
			(*forces)[i] -= rc * c;
			e += rc2;
		}
		(*energies)[ScalarStructure::OTHER] += e * myK;
	}

	void SphericalRestraintSystemForce::parallelEvaluate(const GenericTopology* topo,
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

	void SphericalRestraintSystemForce::getParameters(vector<Parameter>& parameters) const
	{
		//    parameters.push_back(Parameter("-center",Value(myCenter)));
		parameters.push_back(Parameter("-k", Value(myK, ConstraintValueType::NotZero()), Text("external force constant")));
	}

	Force* SphericalRestraintSystemForce::doMake(string& /*errMsg*/, vector<Value> values) const
	{
		return new SphericalRestraintSystemForce(values[0]);
	}

	unsigned int SphericalRestraintSystemForce::numberOfBlocks(const GenericTopology*,
	                                                           const Vector3DBlock* pos)
	{
		return std::min(Parallel::getAvailableNum(), static_cast<int>(pos->size()));
	}
}
