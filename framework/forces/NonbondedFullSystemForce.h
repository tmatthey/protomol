/* -*- c++ -*- */
#ifndef NONBONDEDFULLSYSTEMFORCE_H
#define NONBONDEDFULLSYSTEMFORCE_H

#include "SystemForce.h"
#include "Parallel.h"
#include "NonbondedFullSystemForceBase.h"
#include "simpleTypes.h"

//#include "evaluateBorn.h"
namespace ProtoMol
{
	//_________________________________________________________________ NonbondedFullSystemForce

	template <class TOneAtomPair>
	class NonbondedFullSystemForce: public SystemForce, private NonbondedFullSystemForceBase
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Typedef
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		typedef SemiGenericTopology<typename TOneAtomPair::BoundaryConditions> RealTopologyType;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NonbondedFullSystemForce(): SystemForce(), myCutoff(0.0), myBlockSize(0)
		{
		};

		NonbondedFullSystemForce(Real cutoff, TOneAtomPair oneAtomPair, unsigned int blockSize = defaultBlockSize) :
			SystemForce(), myCutoff(cutoff), myOneAtomPair(oneAtomPair), myBlockSize(blockSize), myCached(false)
		{
		}

		virtual ~NonbondedFullSystemForce()
		{
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NonbondedFullSystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		void doEvaluate(const GenericTopology*,
		                const Vector3DBlock*,
		                Vector3DBlock*,
		                ScalarStructure*,
		                int i0, int i1, int j0, int j1);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class SystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void evaluate(const GenericTopology* to,
		                      const Vector3DBlock* pos,
		                      Vector3DBlock* f,
		                      ScalarStructure* e)
		{
			doEvaluate(to, pos, f, e, 0, (int)to->atoms.size(), 0, (int)to->atoms.size());
		}

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* pos,
		                              Vector3DBlock* f,
		                              ScalarStructure* e)
		{
			if (!myCached)
				splitRangeArea(static_cast<unsigned int>(Parallel::getAvailableNum()), 0, topo->atoms.size(), myFromRange, myToRange);

			for (int i = 0; i < Parallel::getAvailableNum(); i++)
				if (Parallel::next())
					doEvaluate(topo, pos, f, e, myFromRange[i].first, myFromRange[i].second, myToRange[i].first, myToRange[i].second);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Force
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getKeyword() const
		{
			return keyword;
		}

		virtual unsigned int numberOfBlocks(const GenericTopology*, const Vector3DBlock*)
		{
			return Parallel::getAvailableNum();
		}

	private:
		virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void getParameters(std::vector<Parameter>& parameters) const;

		virtual unsigned int getParameterSize() const
		{
			return 2 + TOneAtomPair::getParameterSize();
		}

		virtual std::string getIdNoAlias() const
		{
			return (TOneAtomPair::getId() + " -algorithm " + keyword);
		}

		virtual void uncache()
		{
			myCached = false;
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Real myCutoff;
		TOneAtomPair myOneAtomPair;
		unsigned int myBlockSize;
		bool myCached;

		std::vector<PairUInt> myFromRange;
		std::vector<PairUInt> myToRange;
		std::vector<Vector3D> myLattice;
	};

	//______________________________________________________________________ INLINES


	template <class TOneAtomPair>
	void NonbondedFullSystemForce<TOneAtomPair>::doEvaluate(const GenericTopology* topo,
	                                                        const Vector3DBlock* positions,
	                                                        Vector3DBlock* forces,
	                                                        ScalarStructure* energies,
	                                                        int i0, int i1, int j0, int j1)
	{
		const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
		if (!myCached)
		{
			myLattice = realTopo->boundaryConditions.buildLatticeVectors(myCutoff);
			myCached = true;
		}
		myOneAtomPair.initialize(realTopo, positions, forces, energies, &myLattice);

		for (int blocki = i0; blocki < i1; blocki += myBlockSize)
		{
			int blocki_max = blocki;
			if (blocki_max < j0) blocki_max = j0;
			for (int blockj = blocki_max; blockj < j1; blockj += myBlockSize)
			{
				int istart = blocki;
				int iend = blocki + myBlockSize;
				if (iend > i1) iend = i1;
				for (int i = istart; i < iend; i++)
				{
					int jstart = blockj;
					if (jstart <= i) jstart = i;
					int jend = blockj + myBlockSize;
					if (jend > j1) jend = j1;
					for (int j = jstart; j < jend; j++)
					{
						myOneAtomPair.doOneAtomPair(i, j);
					}
				}
			}
		}
	}

	//   template<class TOneAtomPair
	//   void NonbondedFullSystemForce<TOneAtomPair>::doEvaluate(GenericTopology* topo, 
	// 								const Vector3DBlock* positions, 
	// 								Vector3DBlock* forces, 
	// 								ScalarStructure* energies, 
	// 								int i0, int i1, int j0, int j1) {

	//     doEvaluate(const_cast<GenericTopology*>(static_cast<const GenericTopology*>(topo)),
	// 	       positions, forces, energies);
	//     if (born)
	//       evaluateBorn(myOneAtomPair, topo, positions, forces, energies);
	//   }


	template <class TOneAtomPair>
	Force* NonbondedFullSystemForce<TOneAtomPair>::doMake(std::string& errMsg, std::vector<Value> values) const
	{
		int n = values.size() - 2;
		return new NonbondedFullSystemForce(values[n], TOneAtomPair::make(errMsg, std::vector<Value>(values.begin(), values.end() - 2)), values[n + 1]);
	}


	template <class TOneAtomPair>
	void NonbondedFullSystemForce<TOneAtomPair>::getParameters(std::vector<Parameter>& parameters) const
	{
		myOneAtomPair.getParameters(parameters);
		parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("algorithm cutoff")));
		parameters.push_back(Parameter("-blocksize", Value(myBlockSize, ConstraintValueType::Positive()), defaultBlockSize));
	}
}
#endif /* NONBONDEDFULLSYSTEMFORCE_H */
