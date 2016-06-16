/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFFORCE_H
#define NONBONDEDCUTOFFFORCE_H

#include "Force.h"
#include "Parallel.h"
#include "NonbondedCutoffForceBase.h"

namespace ProtoMol
{
	//_________________________________________________________________ NonbondedCutoffForce

	template <class TCellManager, class TOneAtomPair, class TForce, class TImplForce>
	class NonbondedCutoffForce: public TForce, private NonbondedCutoffForceBase
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Typedef
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		typedef typename TOneAtomPair::BoundaryConditions BoundaryConditions;
		typedef Topology<BoundaryConditions, TCellManager> RealTopologyType;
		typedef typename RealTopologyType::Enumerator EnumeratorType;
		typedef typename RealTopologyType::Enumerator::CellPair CellPairType;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NonbondedCutoffForce(): myCutoff(0.0)
		{
		}

		NonbondedCutoffForce(Real cutoff, TOneAtomPair oneAtomPair) :
			TForce(), myCutoff(cutoff), myOneAtomPair(oneAtomPair)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NonbondedCutoffForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		void doEvaluate(const GenericTopology* topo, unsigned int n);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Force
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual unsigned int numberOfBlocks(const GenericTopology* topo,
		                                    const Vector3DBlock* positions)
		{
			const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
			realTopo->updateCellLists(positions);
			return Parallel::getNumberOfPackages(realTopo->cellLists.size());
		}

		virtual std::string getKeyword() const
		{
			return keyword;
		}

	private:
		virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const;
		virtual void getParameters(std::vector<Parameter>& parameters) const;

		virtual unsigned int getParameterSize() const
		{
			return 1 + TOneAtomPair::getParameterSize();
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
	protected:
		Real myCutoff;
		TOneAtomPair myOneAtomPair;
		EnumeratorType enumerator;
	};


	template <class TCellManager, class TOneAtomPair, class TForce, class TImplForce>
	void NonbondedCutoffForce<TCellManager, TOneAtomPair, TForce, TImplForce>::doEvaluate(const GenericTopology* topo, unsigned int n)
	{
		CellPairType thisPair;
		unsigned int count = 0;
		for (; !enumerator.done(); enumerator.next())
		{
			enumerator.get(thisPair);
			bool notSameCell = enumerator.notSameCell();

			if (!notSameCell)
			{
				count++;
				if (count > n)
					break;
			}
			for (int i = thisPair.first; i != -1; i = topo->atoms[i].cellListNext)
			{
				for (int j = (notSameCell ? thisPair.second : topo->atoms[i].cellListNext); j != -1; j = topo->atoms[j].cellListNext)
				{
					myOneAtomPair.doOneAtomPair(i, j);
				}
			}
		}
	}


	template <class TCellManager, class TOneAtomPair, class TForce, class TImplForce>
	void NonbondedCutoffForce<TCellManager, TOneAtomPair, TForce, TImplForce>::getParameters(std::vector<Parameter>& parameters) const
	{
		myOneAtomPair.getParameters(parameters);
		parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("algorithm cutoff")));
	}

	template <class TCellManager, class TOneAtomPair, class TForce, class TImplForce>
	Force* NonbondedCutoffForce<TCellManager, TOneAtomPair, TForce, TImplForce>::doMake(std::string& errMsg, std::vector<Value> values) const
	{
		return (new TImplForce(values[values.size() - 1], TOneAtomPair::make(errMsg, std::vector<Value>(values.begin(), values.end() - 1))));
	}

	template <class TCellManager, class TOneAtomPair, class TForce, class TImplForce>
	std::string NonbondedCutoffForce<TCellManager, TOneAtomPair, TForce, TImplForce>::getIdNoAlias() const
	{
		return (TOneAtomPair::getId() + " -algorithm " + keyword);
	}
}
#endif /* NONBONDEDCUTOFFFORCE_H */
