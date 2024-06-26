#include "ModifierNPTRattle.h"
#include "Integrator.h"
#include "Topology.h"
#include "ScalarStructure.h"
#include "pmconstants.h"

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//__________________________________________________ ModifierNPTRattleDetails
	ModifierNPTRattleDetails::ModifierNPTRattleDetails(Real eps, int maxIter, int order): ModifierMetaRattle(eps, maxIter, order)
	{
	}

	void ModifierNPTRattleDetails::doExecute()
	{
		// estimate the current error in all velocity constraints
		Real error = calcError();

		// delta_t
		Real dt = getTimestep() / Constant::TIMEFACTOR;

		// multiplicative constants for the box volume velocity and thermostat velocity
		const Real pVol_term = 0.5 * getTimestep() * getEpsilonVel();
		const Real pEta_term = 0.5 * getTimestep() * getEtaVel();
		const Real OneOverN_term = 1.0 + 1.0 / myTopology->atoms.size();

		int iter = 0;
		while (error > myEpsilon)
		{
			for (unsigned int i = 0; i < myListOfConstraints->size(); i++)
			{
				// find the ID#s of the two atoms in the current constraint
				int a1 = (*myListOfConstraints)[i].atom1;
				int a2 = (*myListOfConstraints)[i].atom2;

				// get the ID# of the molecule that this bond belongs to
				int Mol = myTopology->atoms[a1].molecule;

				// reciprocal atomic masses
				Real rM1 = 1 / myTopology->atoms[a1].scaledMass;
				Real rM2 = 1 / myTopology->atoms[a2].scaledMass;
				Real m1_over_M = 1 / (rM1 * myTopology->molecules[Mol].mass);
				Real m2_over_M = 1 / (rM2 * myTopology->molecules[Mol].mass);

				// multiplicative constants due to change in box volume
				// (for constant pressure simulations -- see Equation 56 of G. Kalibaeva,
				//  M. Ferrario, and G. Ciccotti, "Constant pressure-constant temperature molecular
				//  dynamics: a correct constrained NPT ensemble using the molecular virial", Mol. Phys.
				//  101(6), 2003, p. 765-778)
				Real Exp2_atom1 = (1 + pVol_term * OneOverN_term * m1_over_M) * exp(-pEta_term - pVol_term * OneOverN_term * m1_over_M);
				Real Exp2_atom2 = (1 + pVol_term * OneOverN_term * m2_over_M) * exp(-pEta_term - pVol_term * OneOverN_term * m2_over_M);

				// now lets compute the lambdas.
				// compute the current bond vector
				Vector3D rab = (*myPositions)[a1] - (*myPositions)[a2];
				Real rabsq = rab.normSquared();

				// compute the current velocity vector
				Vector3D vab = (*myVelocities)[a1] - (*myVelocities)[a2];

				// dot product of distance and velocity vectors
				Real rvab = rab * vab;

				// compute the change in lambda
				Real gab = -rvab / (dt * (rM1 * Exp2_atom1 + rM2 * Exp2_atom2) * rabsq);
				Vector3D dp = rab * gab;

				// move the velocities based upon the multiplier
				(*myVelocities)[a1] += dp * dt * rM1 * Exp2_atom1;
				(*myVelocities)[a2] -= dp * dt * rM2 * Exp2_atom2;

				// the constraint adds a force to each atom since their positions
				// had to be changed.  This constraint force therefore contributes
				// to the atomic virial.  Note that the molecular virial is independent of
				// any intramolecular constraint forces.
				if (myEnergies->virial()) myEnergies->addVirial(dp * 2, rab);
			}

			// compute the error in all the velocity constraints after this RATTLE iteration
			error = calcError();
			iter ++;
			if (iter > myMaxIter)
			{
				report << warning << "Rattle maxIter = " << myMaxIter
					<< " reached, but still not converged ... error is " << error << endr;
				break;
			}
		}
	}
}
