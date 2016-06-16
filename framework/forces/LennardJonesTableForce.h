/* -*- c++ -*- */
#ifndef LENNARDJONESTABLEFORCE_H
#define LENNARDJONESTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "LookUpTableBase.h"
#include "LennardJonesTableForceBase.h"
#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ LennardJonesTableForce

	template <class TSwitchingFunction, unsigned int PRE, typename TReal=Real>
	class LennardJonesTableForce : public LookUpTableBase<LennardJonesTableForceBase::LookUpValues, PRE, TReal>,
	                               private LennardJonesTableForceBase
	{
	public:
		enum
		{
			DIST_R2=0
		};

		enum
		{
			CUTOFF=1
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		LennardJonesTableForce(): LookUpTableBase<LennardJonesTableForceBase::LookUpValues,
		                                          PRE,
		                                          TReal>(),
		                          myCutoff(0.0),
		                          myCutoff2(0.0),
		                          switchingFunction()
		{
		}

		LennardJonesTableForce(TSwitchingFunction swf): LookUpTableBase<LennardJonesTableForceBase::LookUpValues,
		                                                                PRE,
		                                                                TReal>(0.1,
		                                                                       swf.cutoffSquared(),
		                                                                       2,
		                                                                       LennardJonesTableForceBase::LookUpValues(),
		                                                                       swf,
		                                                                       128),
		                                                myCutoff(swf.cutoff()),
		                                                myCutoff2(swf.cutoffSquared()),
		                                                switchingFunction(swf)
		{
		}

		LennardJonesTableForce(TSwitchingFunction swf,
		                       Real rc): LookUpTableBase<LennardJonesTableForceBase::LookUpValues,
		                                                 PRE,
		                                                 TReal>(0.1,
		                                                        square(rc),
		                                                        2,
		                                                        LennardJonesTableForceBase::LookUpValues(),
		                                                        swf,
		                                                        128),
		                                 myCutoff(rc),
		                                 myCutoff2(rc * rc),
		                                 switchingFunction(swf)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class LennardJonesTableForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void operator()(Real& energy,
		                Real& force,
		                Real distSquared,
		                Real /*rDistSquared*/,
		                const Vector3D& /*diff*/,
		                const GenericTopology* topo,
		                int atom1, int atom2,
		                ExclusionClass excl) const
		{
			const LennardJonesParameters&
				params = topo->lennardJonesParameters(topo->atoms[atom1].type, topo->atoms[atom2].type);

			Real A = (excl != EXCLUSION_MODIFIED ? params.A : params.A14);
			Real B = (excl != EXCLUSION_MODIFIED ? params.B : params.B14);

			Real dt;
			int i;
			this->index(distSquared, i, dt);

			Real a = this->myTable[i + 0] * A + this->myTable[i + 4] * B;
			Real b = this->myTable[i + 1] * A + this->myTable[i + 5] * B;
			Real c = this->myTable[i + 2] * A + this->myTable[i + 6] * B;
			Real d = this->myTable[i + 3] * A + this->myTable[i + 7] * B;

			this->interpolate(a, b, c, d, dt, energy, force);
		}

		static void accumulateEnergy(ScalarStructure* energies, Real energy)
		{
			(*energies)[ScalarStructure::LENNARDJONES] += energy;
		}

		static Real getEnergy(const ScalarStructure* energies)
		{
			return (*energies)[ScalarStructure::LENNARDJONES];
		}

		// Parsing
		static std::string getId()
		{
			return keyword + std::string((!TSwitchingFunction::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunction::getId()));
		}

		static unsigned int getParameterSize()
		{
			return (TSwitchingFunction::CUTOFF ? 0 : 1) + TSwitchingFunction::getParameterSize();
		}

		void getParameters(std::vector<Parameter>& parameters) const
		{
			switchingFunction.getParameters(parameters);
			if (!TSwitchingFunction::CUTOFF)
				parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("cutoff for table look up")));
		}

		static LennardJonesTableForce make(std::string& errMsg, const std::vector<Value>& values)
		{
			if (!TSwitchingFunction::CUTOFF)
				return LennardJonesTableForce(TSwitchingFunction::make(errMsg, std::vector<Value>(values.begin(), values.end() - 1)),
				                              values[values.size() - 1]);
			else
				return LennardJonesTableForce(TSwitchingFunction::make(errMsg, values));
		}

		Real cutoffSquared() const
		{
			return myCutoff2;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Real myCutoff;
		Real myCutoff2;
		TSwitchingFunction switchingFunction;
	};

	//______________________________________________________________________ INLINES
}
#endif /* LENNARDJONESTABLEFORCE_H */
