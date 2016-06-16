/* -*- c++ -*- */
#ifndef ELECTRICFIELDFORCE_H
#define ELECTRICFIELDFORCE_H

#include "SystemForce.h"
#include "ElectricFieldSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"
#include "pmconstants.h"
#include "Parameter.h"
#include "oneAtomContraints.h"

namespace ProtoMol
{
	//_________________________________________________________________ ElectricFieldSystemForce
	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint=NoConstraint>
	class ElectricFieldSystemForce : public SystemForce, private ElectricFieldSystemForceBase
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ElectricFieldSystemForce();
		ElectricFieldSystemForce(TSwitchingFunction sF, Vector3D origin, Real e, Real a);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class ElectricFieldSystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		void doEvaluate(const GenericTopology* topo,
		                const Vector3DBlock* positions,
		                Vector3DBlock* forces,
		                ScalarStructure* energies, int from, int to);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class SystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* position,
		                      Vector3DBlock* forces,
		                      ScalarStructure* energies);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              Vector3DBlock* forces,
		                              ScalarStructure* energies);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Force
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getKeyword() const
		{
			return keyword;
		}

		virtual unsigned int numberOfBlocks(const GenericTopology* topo,
		                                    const Vector3DBlock* pos);
	private:
		virtual Force* doMake(std::string&, std::vector<Value>) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return TConstraint::getPrefixId() + keyword + TConstraint::getPostfixId() + std::string((!TSwitchingFunction::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunction::getId()));
		}

		virtual unsigned int getParameterSize() const
		{
			return 3 + TSwitchingFunction::getParameterSize();
		}

		virtual void getParameters(std::vector<Parameter>&) const;


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		TSwitchingFunction switchingFunction;
		Vector3D myOrigin;
		Real myE;
		Real myA;
		Real myC;
		Real mySquaredCutoff;
	};

	//______________________________________________________________________ INLINES
	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::ElectricFieldSystemForce(): SystemForce(), myOrigin(Vector3D(0.0, 0.0, 0.0)), myE(0.0), myA(0.0), myC(0.0)
	{
	}

	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::ElectricFieldSystemForce(TSwitchingFunction sF, Vector3D origin, Real e, Real a): SystemForce(), switchingFunction(sF), myOrigin(origin), myE(e), myA(a), myC(Constant::SQRTCOULOMBCONSTANT * e), mySquaredCutoff(sF.cutoffSquared())
	{
	}


	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline void ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::getParameters(std::vector<Parameter>& parameters) const
	{
		switchingFunction.getParameters(parameters);
		parameters.push_back(Parameter("-origin", Value(myOrigin)));
		parameters.push_back(Parameter("-e", Value(myE), Text("Electric field source, [e]")));
		parameters.push_back(Parameter("-a", Value(myA, ConstraintValueType::Positive()), Text("Scaling factor")));
	}

	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline Force* ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::doMake(std::string& errMsg, std::vector<Value> values) const
	{
		unsigned int n = TSwitchingFunction::getParameterSize();
		return new ElectricFieldSystemForce(TSwitchingFunction::make(errMsg, std::vector<Value>(values.begin(), values.begin() + n)),
		                                    values[0 + n], values[1 + n], values[2 + n]);
	}

	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline void ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::evaluate(const GenericTopology* topo,
	                                                                                                     const Vector3DBlock* positions,
	                                                                                                     Vector3DBlock* forces,
	                                                                                                     ScalarStructure* energies)
	{
		doEvaluate(topo, positions, forces, energies, 0, (int)positions->size());
	}

	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline void ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::parallelEvaluate(const GenericTopology* topo,
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
				if (to > static_cast<int>(n))
					to = n;
				int from = (n * i) / count;
				doEvaluate(topo, positions, forces, energies, from, to);
			}
		}
	}

	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline void ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::doEvaluate(const GenericTopology* topo,
	                                                                                                       const Vector3DBlock* positions,
	                                                                                                       Vector3DBlock* forces,
	                                                                                                       ScalarStructure* energies, int from, int to)
	{
		const TBoundaryConditions& boundary =
			(dynamic_cast<const SemiGenericTopology<TBoundaryConditions>&>(*topo)).boundaryConditions;

		Real e = 0.0;

		for (int i = from; i < to; i++)
		{
			if (!TConstraint::check(topo, i))
				continue;
			Vector3D pos(boundary.minimalPosition((*positions)[i] - myOrigin));
			Real distSquared = pos.normSquared();

			if (TSwitchingFunction::USE && distSquared > mySquaredCutoff)
				continue;

			Real d = sqrt(distSquared);
			Real r = d * myA + 1;

			// Calculate the force and energy.
			Real energy = myC * topo->atoms[i].scaledCharge / r;
			Real force = -energy * myA / r * (d > Constant::EPSILON ? 1.0 / d : 0.0);

			// Calculate the switched force and energy.
			if (TSwitchingFunction::USE)
			{
				Real switchingValue, switchingDeriv;
				switchingFunction(switchingValue, switchingDeriv, distSquared);
				// This has a - sign because the force is the negative of the 
				// derivative of the energy (divided by the distance between the atoms).
				force = force * switchingValue - energy * switchingDeriv;
				energy = energy * switchingValue;
			}
			e += energy;
			// Add this force into the atom forces.
			Vector3D fi(pos * force);
			(*forces)[i] -= fi;

			// End of force computation.
		}
		(*energies)[ScalarStructure::OTHER] += e;
	}


	template <class TBoundaryConditions,
	          class TSwitchingFunction,
	          class TConstraint>
	inline unsigned int ElectricFieldSystemForce<TBoundaryConditions, TSwitchingFunction, TConstraint>::numberOfBlocks(const GenericTopology*,
	                                                                                                                   const Vector3DBlock* pos)
	{
		return std::min(Parallel::getAvailableNum(), static_cast<int>(pos->size()));
	}
}
#endif /* ELECTRICFIELDFORCE_H */
