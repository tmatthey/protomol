/* -*- c++ -*- */
#ifndef ONEATOMPAIRTWO_H
#define ONEATOMPAIRTWO_H


#include "Topology.h"
#include "Parameter.h"
#include "oneAtomContraints.h"

namespace ProtoMol
{
	//_____________________________________________________________ OneAtomPairTwo
	/**
	 * Computes the interaction for a given force between two atoms with the
	 * template arguments defining the boundary conditions, two switching functions,
	 * two potentials and optional constraint.
	 */
	template <class TBoundaryConditions,
	          class TSwitchingFunctionFirst,
	          class TNonbondedForceFirst,
	          class TSwitchingFunctionSecond,
	          class TNonbondedForceSecond,
	          class TConstraint=NoConstraint>
	class OneAtomPairTwo
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Typedef & sub classes
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		typedef TBoundaryConditions BoundaryConditions;
		// Make the boundary conditions visible

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OneAtomPairTwo()
		{
		}

		OneAtomPairTwo(TNonbondedForceFirst f1,
		               TSwitchingFunctionFirst sF1,
		               TNonbondedForceSecond f2,
		               TSwitchingFunctionSecond sF2):
			switchingFunctionFirst(sF1),
			nonbondedForceFunctionFirst(f1),
			switchingFunctionSecond(sF2),
			nonbondedForceFunctionSecond(f2), mySquaredCutoff(std::max(Cutoff<TNonbondedForceFirst::CUTOFF>::cutoff(sF1, f1), Cutoff<TNonbondedForceSecond::CUTOFF>::cutoff(sF2, f2)))
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class OneAtomPairTwo
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void initialize(const SemiGenericTopology<TBoundaryConditions>* topo,
		                const Vector3DBlock* pos, Vector3DBlock* f, ScalarStructure* e)
		{
			realTopo = topo;
			positions = pos;
			forces = f;
			energies = e;
		}

		void doOneAtomPair(const int i, const int j);

		// Computes the force and energy for atom i and j.

		void getParameters(std::vector<Parameter>& parameters) const
		{
			nonbondedForceFunctionFirst.getParameters(parameters);
			switchingFunctionFirst.getParameters(parameters);
			nonbondedForceFunctionSecond.getParameters(parameters);
			switchingFunctionSecond.getParameters(parameters);
		}

		static unsigned int getParameterSize()
		{
			return
				TNonbondedForceFirst::getParameterSize() +
				TSwitchingFunctionFirst::getParameterSize() +
				TNonbondedForceSecond::getParameterSize() +
				TSwitchingFunctionSecond::getParameterSize();
		}

		static OneAtomPairTwo make(std::string& errMsg, std::vector<Value> values)
		{
			unsigned int l1 = TNonbondedForceFirst::getParameterSize();
			unsigned int l2 = TSwitchingFunctionFirst::getParameterSize() + l1;
			unsigned int l3 = TNonbondedForceSecond::getParameterSize() + l2;
			return OneAtomPairTwo(TNonbondedForceFirst::make(errMsg, std::vector<Value>(values.begin(), values.begin() + l1)),
			                      TSwitchingFunctionFirst::make(errMsg, std::vector<Value>(values.begin() + l1, values.begin() + l2)),
			                      TNonbondedForceSecond::make(errMsg, std::vector<Value>(values.begin() + l2, values.begin() + l3)),
			                      TSwitchingFunctionSecond::make(errMsg, std::vector<Value>(values.begin() + l3, values.end())));
		}

		static std::string getId()
		{
			return
				TConstraint::getPrefixId() + headString(TNonbondedForceFirst::getId()) + TConstraint::getPostfixId() + " " +
				TConstraint::getPrefixId() + headString(TNonbondedForceSecond::getId()) + TConstraint::getPostfixId() +
				(tailString(TNonbondedForceFirst::getId()).empty() ? "" : " ") + tailString(TNonbondedForceFirst::getId()) +
				(tailString(TNonbondedForceSecond::getId()).empty() ? "" : " ") + tailString(TNonbondedForceSecond::getId()) +
				std::string((!TSwitchingFunctionFirst::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunctionFirst::getId())) +
				std::string((!TSwitchingFunctionSecond::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunctionSecond::getId()));
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
		Real mySquaredCutoff;
	};

	//______________________________________________________________________ INLINES

	template <class TBoundaryConditions,
	          class TSwitchingFunctionFirst,
	          class TNonbondedForceFirst,
	          class TSwitchingFunctionSecond,
	          class TNonbondedForceSecond,
	          class TConstraint>
	inline void OneAtomPairTwo<TBoundaryConditions,
	                           TSwitchingFunctionFirst,
	                           TNonbondedForceFirst,
	                           TSwitchingFunctionSecond,
	                           TNonbondedForceSecond,
	                           TConstraint>::doOneAtomPair(const int i, const int j)
	{
		if (TConstraint::PRE_CHECK)
			if (!TConstraint::check(realTopo, i, j))
				return;

		// Get atom distance.
		Real distSquared;
		Vector3D diff(realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j], distSquared));
		// Do switching function rough test, if necessary.
		if ((TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE || TNonbondedForceFirst::CUTOFF || TNonbondedForceSecond::CUTOFF) && distSquared > mySquaredCutoff)
			return;

		// Check for an exclusion.
		int mi = realTopo->atoms[i].molecule;
		int mj = realTopo->atoms[j].molecule;
		bool same = (mi == mj);
		ExclusionClass excl = (same ? realTopo->exclusions.check(i, j) : EXCLUSION_NONE);
		if (excl == EXCLUSION_FULL)
			return;

		// Calculate the force and energy.
		Real rDistSquared = ((TNonbondedForceFirst::DIST_R2 || TNonbondedForceSecond::DIST_R2) ? 1.0 / distSquared : 1.0);
		Real energy1, force1, energy2, force2;
		nonbondedForceFunctionFirst(energy1, force1, distSquared, rDistSquared, diff, realTopo, i, j, excl);
		nonbondedForceFunctionSecond(energy2, force2, distSquared, rDistSquared, diff, realTopo, i, j, excl);
		//Report::report << "\t"<<i << "\t"<<j<<Report::endr;
		// Calculate the switched force and energy.
		if (TSwitchingFunctionFirst::MODIFY || TSwitchingFunctionSecond::MODIFY)
		{
			Real switchingValue, switchingDeriv;

			switchingFunctionFirst(switchingValue, switchingDeriv, distSquared);
			force1 = force1 * switchingValue - energy1 * switchingDeriv;
			energy1 = energy1 * switchingValue;

			switchingFunctionSecond(switchingValue, switchingDeriv, distSquared);
			force2 = force2 * switchingValue - energy2 * switchingDeriv;
			energy2 = energy2 * switchingValue;
		}

		// Add this energy into the total system energy.
		nonbondedForceFunctionFirst.accumulateEnergy(energies, energy1);
		nonbondedForceFunctionSecond.accumulateEnergy(energies, energy2);
		// Add this force into the atom forces.
		Vector3D fij(diff * (force1 + force2));
		(*forces)[i] -= fij;
		(*forces)[j] += fij;

		// compute the vector between molecular centers of mass
		if (!same && energies->molecularVirial())
		{
			// Add to the atomic and molecular virials
			energies->addVirial(fij, diff, realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
			                                                                              realTopo->molecules[mj].position));
		}
		else if (energies->virial())
		{
			energies->addVirial(fij, diff);
		}
		// End of force computation.
		if (TConstraint::POST_CHECK)
			TConstraint::check(realTopo, i, j, diff, energy1 + energy2, fij);
	}
}

#endif /* ONEATOMPAIRTWO_H */
