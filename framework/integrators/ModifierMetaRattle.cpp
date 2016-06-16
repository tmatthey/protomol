#include "ModifierMetaRattle.h"
#include "Topology.h"

namespace ProtoMol
{
	//__________________________________________________ ModifierMetaRattle
	ModifierMetaRattle::ModifierMetaRattle(Real eps, int maxIter, int order): ModifierMetaRattleShake(eps, maxIter, order)
	{
	}


	Real ModifierMetaRattle::calcError() const
	{
		// the error for the RATTLE algorithm is defined as fabs( [v1-v2] * [r1-r2] ) < tolerance
		// which is the constraint imposed upon the velocities by RATTLE.  It is this
		// constraint upon the velocities that allows us to compute the
		// multipliers (lambdas) at time t + delta_t

		Real error = 0;
		for (unsigned int i = 0; i < myListOfConstraints->size(); i++)
		{
			int a1 = (*myListOfConstraints)[i].atom1;
			int a2 = (*myListOfConstraints)[i].atom2;
			Vector3D vab = (*myVelocities)[a1] - (*myVelocities)[a2];
			Vector3D rab = (*myPositions)[a1] - (*myPositions)[a2];
			Real err = fabs(rab * vab);
			error += err;
		}
		return error /= myListOfConstraints->size();
	}
}
