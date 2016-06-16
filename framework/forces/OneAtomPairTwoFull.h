/* -*- c++ -*- */
#ifndef ONEATOMPAIRTWOFULL_H
#define ONEATOMPAIRTWOFULL_H

#include "Topology.h"
#include "Parameter.h"
#include "oneAtomContraints.h"

namespace ProtoMol
{
	//_____________________________________________________________ OneAtomPairTwoFull

	template <class TBoundaryConditions,
	          class TSwitchingFunctionFirst,
	          class TNonbondedForceFirst,
	          class TSwitchingFunctionSecond,
	          class TNonbondedForceSecond,
	          class TConstraint=NoConstraint>
	class OneAtomPairTwoFull
	{
		// Computes the interaction for a given force between to atoms.
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Typedef
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		typedef TBoundaryConditions BoundaryConditions;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OneAtomPairTwoFull()
		{
		}

		OneAtomPairTwoFull(TNonbondedForceFirst f1,
		                   TNonbondedForceSecond f2,
		                   TSwitchingFunctionFirst sF1,
		                   TSwitchingFunctionSecond sF2):
			switchingFunctionFirst(sF1),
			nonbondedForceFunctionFirst(f1),
			switchingFunctionSecond(sF2),
			nonbondedForceFunctionSecond(f2)
		{
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class OneAtomPairTwoFull
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void initialize(const SemiGenericTopology<TBoundaryConditions>* topo,
		                const Vector3DBlock* pos,
		                Vector3DBlock* f,
		                ScalarStructure* e,
		                const std::vector<Vector3D>* l)
		{
			realTopo = topo;
			positions = pos;
			forces = f;
			energies = e;
			lattice = l;
		}

		void doOneAtomPair(const int i, const int j);

		// Computes the force and energy for atom i and j.

		void getParameters(std::vector<Parameter>& parameters) const
		{
			nonbondedForceFunctionFirst.getParameters(parameters);
			nonbondedForceFunctionSecond.getParameters(parameters);
			switchingFunctionFirst.getParameters(parameters);
			switchingFunctionSecond.getParameters(parameters);
		}

		static unsigned int getParameterSize()
		{
			return
				TNonbondedForceFirst::getParameterSize() +
				TNonbondedForceSecond::getParameterSize() +
				TSwitchingFunctionFirst::getParameterSize() +
				TSwitchingFunctionSecond::getParameterSize();
		}

		static OneAtomPairTwoFull make(std::string& errMsg, std::vector<Value> values)
		{
			unsigned int l1 = TNonbondedForceFirst::getParameterSize();
			unsigned int l2 = TNonbondedForceSecond::getParameterSize() + l1;
			unsigned int l3 = TSwitchingFunctionFirst::getParameterSize() + l2;
			return OneAtomPairTwoFull(TNonbondedForceFirst::make(errMsg, std::vector<Value>(values.begin(), values.begin() + l1)),
			                          TNonbondedForceSecond::make(errMsg, std::vector<Value>(values.begin() + l1, values.begin() + l2)),
			                          TSwitchingFunctionFirst::make(errMsg, std::vector<Value>(values.begin() + l2, values.begin() + l3)),
			                          TSwitchingFunctionSecond::make(errMsg, std::vector<Value>(values.begin() + l3, values.end())));
		}

		static std::string getId()
		{
			return TConstraint::getPrefixId() + TNonbondedForceFirst::getId() + TConstraint::getPostfixId() + " " +
				TConstraint::getPrefixId() + TNonbondedForceSecond::getId() + TConstraint::getPostfixId() + std::string((!(TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE)) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunctionFirst::getId() + " -switchingFunction " + TSwitchingFunctionSecond::getId()));
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		const SemiGenericTopology<TBoundaryConditions>* realTopo;
		const Vector3DBlock* positions;
		Vector3DBlock* forces;
		ScalarStructure* energies;
		TSwitchingFunctionFirst switchingFunctionFirst;
		TNonbondedForceFirst nonbondedForceFunctionFirst;
		TSwitchingFunctionSecond switchingFunctionSecond;
		TNonbondedForceSecond nonbondedForceFunctionSecond;
		const std::vector<Vector3D>* lattice;
	};

	//______________________________________________________________________ INLINES

	template <class TBoundaryConditions,
	          class TSwitchingFunctionFirst,
	          class TNonbondedForceFirst,
	          class TSwitchingFunctionSecond,
	          class TNonbondedForceSecond,
	          class TConstraint>
	inline void OneAtomPairTwoFull<TBoundaryConditions,
	                               TSwitchingFunctionFirst,
	                               TNonbondedForceFirst,
	                               TSwitchingFunctionSecond,
	                               TNonbondedForceSecond,
	                               TConstraint>::doOneAtomPair(const int i, const int j)
	{
		if (TConstraint::PRE_CHECK)
			if (!TConstraint::check(realTopo, i, j))
				return;

		// Do we have something to do?
		bool same = (i == j);
		if (same && lattice->empty())
			return;

		Vector3D diffMinimal(realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j]));
		if (!same)
		{
			Real distSquared = diffMinimal.normSquared();
			// Do switching function rough test, if necessary.
			if ((TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE) && !switchingFunctionFirst.roughTest(distSquared)
				&& !switchingFunctionSecond.roughTest(distSquared))
				return;

			// Check for an exclusion.
			ExclusionClass excl = realTopo->exclusions.check(i, j);
			if (excl != EXCLUSION_FULL)
			{
				// Calculate the force and energy.
				Real rawEnergy1, rawForce1, rawEnergy2, rawForce2;
				Real rDistSquared = 1.0 / distSquared;
				nonbondedForceFunctionFirst(rawEnergy1, rawForce1, distSquared, rDistSquared, diffMinimal, realTopo, i, j, excl);
				nonbondedForceFunctionSecond(rawEnergy2, rawForce2, distSquared, rDistSquared, diffMinimal, realTopo, i, j, excl);

				// Calculate the switched force and energy.
				Real energy1, force1, energy2, force2;
				if (TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE)
				{
					Real switchingValue1, switchingDeriv1;
					Real switchingValue2, switchingDeriv2;
					switchingFunctionFirst(switchingValue1, switchingDeriv1, distSquared);
					switchingFunctionSecond(switchingValue2, switchingDeriv2, distSquared);
					energy1 = rawEnergy1 * switchingValue1;
					energy2 = rawEnergy2 * switchingValue2;
					// This has a - sign because the force is the negative of the 
					// derivative of the energy (divided by the distance between the atoms).
					force1 = rawForce1 * switchingValue1 - rawEnergy1 * switchingDeriv1;
					force2 = rawForce2 * switchingValue2 - rawEnergy2 * switchingDeriv2;
				}
				else
				{
					energy1 = rawEnergy1;
					energy2 = rawEnergy2;
					force1 = rawForce1;
					force2 = rawForce2;
				}
				// Add this energy into the total system energy.
				nonbondedForceFunctionFirst.accumulateEnergy(energies, energy1);
				nonbondedForceFunctionSecond.accumulateEnergy(energies, energy2);
				// Add this force into the atom forces.
				Vector3D fij = -diffMinimal * (force1 + force2);
				(*forces)[i] += fij;
				(*forces)[j] -= fij;

				// compute the vector between molecular centers of mass
				Vector3D molDiff = realTopo->boundaryConditions.minimalDifference(realTopo->molecules[realTopo->atoms[i].molecule].position,
				                                                                  realTopo->molecules[realTopo->atoms[j].molecule].position);

				// Add to the atomic and molecular virials
				energies->addVirial(fij, -diffMinimal, -molDiff);

				if (TConstraint::POST_CHECK)
					TConstraint::check(realTopo, i, j, diffMinimal, energy1 + energy2, fij);
			}
		}

		for (unsigned int k = 0; k < lattice->size(); k++)
		{
			Vector3D diff(diffMinimal + (*lattice)[k]);
			// Get atom distance.
			Real distSquared = diff.normSquared();

			// Do switching function rough test, if necessary.
			if ((TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE) && !switchingFunctionFirst.roughTest(distSquared)
				&& !switchingFunctionSecond.roughTest(distSquared))
				continue;
			// Calculate the force and energy.
			Real rawEnergy1, rawForce1, rawEnergy2, rawForce2;
			Real rDistSquared = 1.0 / distSquared;
			nonbondedForceFunctionFirst(rawEnergy1, rawForce1, distSquared, rDistSquared, diff, realTopo, i, j, EXCLUSION_NONE);
			nonbondedForceFunctionSecond(rawEnergy2, rawForce2, distSquared, rDistSquared, diff, realTopo, i, j, EXCLUSION_NONE);

			// Calculate the switched force and energy.
			Real energy1, force1, energy2, force2;
			if (TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE)
			{
				Real switchingValue1, switchingDeriv1;
				Real switchingValue2, switchingDeriv2;
				switchingFunctionFirst(switchingValue1, switchingDeriv1, distSquared);
				switchingFunctionSecond(switchingValue2, switchingDeriv2, distSquared);
				energy1 = rawEnergy1 * switchingValue1;
				energy2 = rawEnergy2 * switchingValue2;
				// This has a - sign because the force is the negative of the 
				// derivative of the energy (divided by the distance between the atoms).
				force1 = rawForce1 * switchingValue1 - rawEnergy1 * switchingDeriv1;
				force2 = rawForce2 * switchingValue2 - rawEnergy2 * switchingDeriv2;
			}
			else
			{
				energy1 = rawEnergy1;
				energy2 = rawEnergy2;
				force1 = rawForce1;
				force2 = rawForce2;
			}
			// Correct the energy by factor 1/2 when same atom since
			// there is only one pair (i,j) with i==j, where 
			// there are two pairs with same contribution with i !=j
			if (same)
			{
				energy1 /= 2;
				energy2 /= 2;
			}
			else
			{
				// Add this force into the atom forces.
				Vector3D fij = -diff * (force1 + force2);
				(*forces)[i] += fij;
				(*forces)[j] -= fij;

				// compute the vector between molecular centers of mass
				Vector3D molDiff = realTopo->boundaryConditions.minimalDifference(realTopo->molecules[realTopo->atoms[i].molecule].position,
				                                                                  realTopo->molecules[realTopo->atoms[j].molecule].position) +
					(*lattice)[k];

				// Add to the atomic and molecular virials
				energies->addVirial(fij, -diff, -molDiff);
			}

			// Add this energy into the total system energy.
			nonbondedForceFunctionFirst.accumulateEnergy(energies, energy1);
			nonbondedForceFunctionSecond.accumulateEnergy(energies, energy2);

			if (TConstraint::POST_CHECK)
				TConstraint::check(realTopo, i, j, diff, energy1 + energy2, -diff * (force1 + force2));
		}
		// End of force computation.
	}
}
#endif /* ONEATOMPAIRTWO_H */
