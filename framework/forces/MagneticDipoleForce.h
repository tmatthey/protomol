/* -*- c++ -*- */
#ifndef MAGNETICDIPOLEFORCE_H
#define MAGNETICDIPOLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "ExclusionTable.h"
#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ MagneticDipoleForce
	// This force calculates the force between two dipoles created by a magnetic
	// field. It returns the absolute value of the force and changes the
	// distance-vector (diff) since the force is usually not radial. It includes the
	// mirror effect (1. aprox, i.e. only one mirror-image on each side of the
	// boundary. It also assumes that the dipoles lie in the middle between the two
	// boundaries.
	// http://www.ife.no/media/504_Thesis.pdf
	// http://www.ife.no/tillegg/index.jsp?projectId=3047&tilleggId=1035
	//

	class MagneticDipoleForce
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
		MagneticDipoleForce();
		MagneticDipoleForce(Real chi, Real radius, Real omega, Real phi, Real Hx, Real Hy, Real Hz, Real D);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class MagneticDipoleForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void operator()(Real& energy,
		                Real& force,
		                Real distSquared,
		                Real rDistSquared,
		                Vector3D& diff,
		                const GenericTopology* topo,
		                int atom1, int atom2,
		                ExclusionClass excl) const;

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
			return 8;
		}

		void getParameters(std::vector<Parameter>& parameters) const;
		static MagneticDipoleForce make(std::string& errMsg, const std::vector<Value>& values);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		Real myChi, myR, myOmega, myPhi, myHx, myHy, myHz, myD;
		//Short cuts
		Real volum, expfactor, realChi, kappa;
	};

	//______________________________________________________________________ INLINES
}
#endif /* MAGNETICDIPOLEFORCE_H */
