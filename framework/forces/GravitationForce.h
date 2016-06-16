/* -*- c++ -*- */
#ifndef GRAVITATIONFORCE_H
#define GRAVITATIONFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "ExclusionTable.h"
#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ GravitationForce


	class GravitationForce
	{
	public:
		enum
		{
			DIST_R2=1
		};

		enum
		{
			CUTOFF=0
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		GravitationForce(): myG(0.0)
		{
		}

		GravitationForce(Real g): myG(g)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class GravitationForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void operator()(Real& energy,
		                Real& force,
		                Real /*distSquared*/,
		                Real rDistSquared,
		                const Vector3D&,
		                const GenericTopology* topo,
		                int atom1, int atom2,
		                ExclusionClass) const
		{
			Real m1 = topo->atoms[atom1].scaledMass;
			Real m2 = topo->atoms[atom2].scaledMass;

			//Gravitational force:
			energy = -m1 * m2 * myG * sqrt(rDistSquared);
			force = energy * rDistSquared;
		}

		static void accumulateEnergy(ScalarStructure* energies, Real energy)
		{
			(*energies)[ScalarStructure::OTHER] += energy;
		}

		static Real getEnergy(const ScalarStructure* energies)
		{
			return (*energies)[ScalarStructure::OTHER];
		}

		// Parsing
		static std::string getId()
		{
			return keyword;
		}

		static unsigned int getParameterSize()
		{
			return 1;
		}

		void getParameters(std::vector<Parameter>& parameters) const;
		static GravitationForce make(std::string&, const std::vector<Value>& values);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		Real myG;
	};

	//______________________________________________________________________ INLINES
}
#endif /* GRAVITATIONFORCE_H */
