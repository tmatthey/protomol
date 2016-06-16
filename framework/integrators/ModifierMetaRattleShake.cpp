#include "ModifierMetaRattleShake.h"
#include "Topology.h"
#include "topologyutilities.h"

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//__________________________________________________ ModifierMetaRattleShake
	ModifierMetaRattleShake::ModifierMetaRattleShake(Real eps, int maxIter, int order): Modifier(order),
	                                                                                    myEpsilon(eps),
	                                                                                    myMaxIter(maxIter),
	                                                                                    myListOfConstraints(NULL)
	{
	}

	void ModifierMetaRattleShake::doInitialize()
	{
		myLastPositions = (*myPositions);

		// ... maybe it's a second inittialize, add the old back.
		myTopology->degreesOfFreedom += myTopology->bondRattleShakeConstraints.size();

		buildRattleShakeBondConstraintList(myTopology, myTopology->bondRattleShakeConstraints);
		// this list contains bonded pairs, and UB-bonded pairs excluding
		// (heavy atom)-H pairs and (heavy)-(heavy) pairs

		// subtract the # of constraints from the # of degrees of freedom of the system
		// This is needed so that we get the correct temperature
		myTopology->degreesOfFreedom -= myTopology->bondRattleShakeConstraints.size();

		myListOfConstraints = &(myTopology->bondRattleShakeConstraints);
	}
}
