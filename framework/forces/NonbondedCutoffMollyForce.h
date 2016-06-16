/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFMOLLYFORCE_H
#define NONBONDEDCUTOFFMOLLYFORCE_H

#include "MollyForce.h"
#include "NonbondedCutoffForce.h"

namespace ProtoMol
{
	//_________________________________________________________________ NonbondedCutoffMollyForce

	template <class TCellManager, class TOneAtomPair>
	class NonbondedCutoffMollyForce: public NonbondedCutoffForce<TCellManager, TOneAtomPair, MollyForce, NonbondedCutoffMollyForce<TCellManager, TOneAtomPair>>
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
		NonbondedCutoffMollyForce(): NonbondedCutoffForce<TCellManager, TOneAtomPair, MollyForce, NonbondedCutoffMollyForce>()
		{
		}

		NonbondedCutoffMollyForce(Real cutoff, TOneAtomPair oneAtomPair) : NonbondedCutoffForce<TCellManager, TOneAtomPair, MollyForce, NonbondedCutoffMollyForce>(cutoff, oneAtomPair)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class MollyForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* positions,
		                      std::vector<ReducedHessAngle>* angleFilter);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              std::vector<ReducedHessAngle>* angleFilter);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
	private:
	};


	template <class TCellManager, class TOneAtomPair>
	void NonbondedCutoffMollyForce<TCellManager, TOneAtomPair>::evaluate(const GenericTopology* topo,
	                                                                     const Vector3DBlock* positions,
	                                                                     std::vector<ReducedHessAngle>* angleFilter)
	{
		const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);

		this->myOneAtomPair.initialize(realTopo, positions, angleFilter);
		realTopo->updateCellLists(positions);
		this->enumerator.initialize(realTopo, this->myCutoff);
		this->doEvaluate(topo, realTopo->cellLists.size());
	}

	template <class TCellManager, class TOneAtomPair>
	void NonbondedCutoffMollyForce<TCellManager, TOneAtomPair>::parallelEvaluate(const GenericTopology* topo,
	                                                                             const Vector3DBlock* positions,
	                                                                             std::vector<ReducedHessAngle>* angleFilter)
	{
		const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);

		this->myOneAtomPair.initialize(realTopo, positions, angleFilter);
		realTopo->updateCellLists(positions);
		this->enumerator.initialize(realTopo, this->myCutoff);

		unsigned int n = realTopo->cellLists.size();
		unsigned int count = numberOfBlocks(realTopo, positions);

		for (unsigned int i = 0; i < count; i++)
		{
			unsigned int l = (n * (i + 1)) / count - (n * i) / count;
			if (Parallel::next())
				this->doEvaluate(topo, l);
			else
				this->enumerator.nextNewPair(l);
		}
	}
}
#endif /* NONBONDEDCUTOFFMOLLYFORCE_H */
