/* -*- c++ -*- */
#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "SemiGenericTopology.h"
#include "CellListEnumerator.h"
#include "buildCellLists.h"
#include "stringutilities.h"
#include "pmconstants.h"

namespace ProtoMol
{
	//_________________________________________________________________ Topology
	/**
	 * Implementation of the topology of a systems with a given boundary conditions 
	 * and cell manager. 
	 */

	template <class TBoundaryConditions, class TCellManager>
	class Topology: public SemiGenericTopology<TBoundaryConditions>
	{
	public:

		typedef CellListEnumerator<TBoundaryConditions, TCellManager> Enumerator;
		typedef TCellManager CellManager;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Topology(): SemiGenericTopology<TBoundaryConditions>()
		{
		}

		Topology(Real csf,
		         const ExclusionType& e,
		         const TBoundaryConditions& b,
		         const TCellManager& c): SemiGenericTopology<TBoundaryConditions>(csf, e, b),
		                                 cellManager(c)
		{
		}

		virtual ~Topology()
		{
		};

		/// marks the the cell list as not valid any more.
		virtual void uncacheCellList()
		{
			cellLists.uncache();
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Topology
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		/// invokces an update of the cell list, if necessary
		void updateCellLists(const Vector3DBlock* positions) const
		{
			if (!cellLists.valid)
			{
				if (this->boundaryConditions.PERIODIC)
				{
					this->min = this->boundaryConditions.getMin();
					this->max = this->boundaryConditions.getMax();
				}
				else
				{
					positions->boundingbox(this->min, this->max);
				}
				buildCellLists(this, positions);
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void getParameters(std::vector<Parameter>& parameters) const
		{
			parameters.push_back(Parameter("coulombScalingFactor", Value(this->coulombScalingFactor), 1.0));
			parameters.push_back(Parameter("exclude",
			                               Value(this->exclude.getString(), ConstraintValueType::NotEmpty()),
			                               Text(std::string("exclusion scheme (") + ExclusionType::getPossibleValues() + std::string(")"))));
			this->boundaryConditions.getParameters(parameters);
			this->cellManager.getParameters(parameters);
		}

		virtual unsigned int getParameterSize() const
		{
			return (this->boundaryConditions.getParameterSize() + this->cellManager.getParameterSize() + 2);
		}

		virtual std::string getIdNoAlias() const
		{
			return std::string(TBoundaryConditions::keyword + TCellManager::keyword);
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From Topology
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string print(const Vector3DBlock* positions = NULL) const
		{
			unsigned int count = 0;
			unsigned int countMult = 0;
			unsigned int countH20 = 0;
			for (unsigned int i = 0; i < this->molecules.size(); ++i)
			{
				count += this->molecules[i].pairs.size();
				countH20 += (this->molecules[i].water ? 1 : 0);
			}
			for (unsigned int i = 0; i < this->dihedrals.size(); ++i)
			{
				countMult += this->dihedrals[i].multiplicity;
			}
			std::string res;
			res +=
				"Atoms                : " + toString(this->atoms.size()) + "\n" +
				"AtomTypes            : " + toString(this->atomTypes.size()) + "\n" +
				"Bonds                : " + toString(this->bonds.size()) + "\n" +
				"Angles               : " + toString(this->angles.size()) + "\n" +
				"Dihedrals            : " + toString(this->dihedrals.size()) + " (" + toString(countMult) + ")\n" +
				"Impropers            : " + toString(this->impropers.size()) + "\n" +
				"Molecules            : " + toString(this->molecules.size()) + "\n" +
				"Water                : " + toString(countH20) + "\n" +
				"Pairs                : " + toString(count) + "\n" +
				"Degree of freedom    : " + toString(this->degreesOfFreedom) + "\n";
			if (positions != NULL)
				res += "Molecule pair dist's : " + toString(this->checkMoleculePairDistances(*positions)) + "\n";
			else
				res += "Molecule pair dist's : " + toString(this->minimalMolecularDistances) + "\n";

			res +=
				"Exclusion pairs      : " + toString(this->exclusions.getTable().size()) + "\n" +
				"Max exclusion dist   : " + toString(this->exclusions.getMaxDelta()) + "\n" +
				"Time                 : " + toString(this->time) + "\n" +
				"CellManager          : " + TCellManager::keyword + "\n" +
				"BoundaryConditions   : " + TBoundaryConditions::keyword + (this->boundaryConditions.isOrthogonal() ? " orthogonal" : "") + "\n";
			std::vector<Parameter> parameters;
			this->getParameters(parameters);
			for (unsigned int i = 0; i < parameters.size(); i++)
				res += getRightFill(parameters[i].keyword, 21) + ": " + parameters[i].value.getString() + "\n";

			if (positions != NULL)
			{
				Real v = this->getVolume(*positions);
				if (v > Constant::EPSILON)
				{
					res += "Atom density         : " + toString(this->atoms.size() / v) + "\n";
					res += "Atom/cell            : " + toString(this->atoms.size() * this->cellManager.getCellVolume() / v) + "\n";
					res += "Real cell size       : (" + toString(this->cellManager.getRealCellSize().x) + "," + toString(this->cellManager.getRealCellSize().y) + "," + toString(this->cellManager.getRealCellSize().z) + ")\n";
					res += "Cell dimension       : (" + toString(this->cellLists.getDimX()) + "," + toString(this->cellLists.getDimY()) + "," + toString(this->cellLists.getDimZ()) + ")\n";
				}
			}

			Vector3D a, b;
			this->getBoundaryConditionsBox(a, b);
			res +=
				"Simulation box       : (" + toString(a.x) + "," + toString(a.y) + "," + toString(a.z) + ")-" +
				"(" + toString(b.x) + "," + toString(b.y) + "," + toString(b.z) + ") " +
				"(" + toString(b.x - a.x) + "," + toString(b.y - a.y) + "," + toString(b.z - a.z) + ")";
			if (positions != NULL)
			{
				this->getBoundingbox(*positions, a, b);
				res +=
					"\nParticle             : (" + toString(a.x) + "," + toString(a.y) + "," + toString(a.z) + ")-" +
					"(" + toString(b.x) + "," + toString(b.y) + "," + toString(b.z) + ") " +
					"(" + toString(b.x - a.x) + "," + toString(b.y - a.y) + "," + toString(b.z - a.z) + ")";
				positions->boundingbox(a, b);
				res +=
					"\nParticle extended    : (" + toString(a.x) + "," + toString(a.y) + "," + toString(a.z) + ")-" +
					"(" + toString(b.x) + "," + toString(b.y) + "," + toString(b.z) + ") " +
					"(" + toString(b.x - a.x) + "," + toString(b.y - a.y) + "," + toString(b.z - a.z) + ")";
			}

			return res;
		}

	private:
		virtual GenericTopology* doMake(std::string& errMsg, std::vector<Value> values) const
		{
			Real csf;
			if (!values[0].get(csf))
				errMsg += " coulombScalingFactor \'" + values[0].getString() + "\' not valid.";
			ExclusionType e(values[1].getString());
			if (!e.valid())
				errMsg += " Exclusion \'" + values[1].getString() + "\' not recognized, possible values are: " + ExclusionType::getPossibleValues(",") + ".";

			if (!errMsg.empty())
				return NULL;

			unsigned int n = TBoundaryConditions::getParameterSize();
			return new Topology<TBoundaryConditions, TCellManager>(csf, e,
			                                                       TBoundaryConditions::make(errMsg, std::vector<Value>(values.begin() + 2, values.begin() + n + 2)),
			                                                       TCellManager::make(errMsg, std::vector<Value>(values.begin() + n + 2, values.end())));
		}


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		CellManager cellManager;
		mutable typename CellManager::CellListStructure cellLists;
	};
}
#endif /* TOPOLOGY_H */
