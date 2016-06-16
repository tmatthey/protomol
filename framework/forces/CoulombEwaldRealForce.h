/* -*- c++ -*- */
#ifndef COULOMBEWALDREALFORCE_H
#define COULOMBEWALDREALFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include <string>

//#define USE_COULOMBEWALDREAL_EXACT_ERF
namespace ProtoMol
{
	//_________________________________________________________________ CoulombEwaldRealForce
	class CoulombEwaldRealForce
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Replacement potential of the real part for Ewald or PME.
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
	public:
#ifdef USE_COULOMBEWALDREAL_EXACT_ERF
    enum {DIST_R2=1};
#else
		enum
		{
			DIST_R2=0
		};
#endif
		enum
		{
			CUTOFF=0
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		CoulombEwaldRealForce();
		CoulombEwaldRealForce(Real a);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class CoulombEwaldRealForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void operator()(Real& energy,
		                Real& force,
		                Real distSquared,
		                Real
#ifdef USE_COULOMBEWALDREAL_EXACT_ERF
		    rDistSquared
#endif
		                ,
		                const Vector3D&,
		                const GenericTopology* topo,
		                int atom1, int atom2,
		                ExclusionClass excl) const
		{
			Real qq = topo->atoms[atom1].scaledCharge * topo->atoms[atom2].scaledCharge;
			if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED)
				qq *= topo->coulombScalingFactor;
			Real r = sqrt(distSquared);
#ifndef USE_COULOMBEWALDREAL_EXACT_ERF
			Real rr = 1.0 / r;
			Real ar = myAlpha * r;
			Real e = qq * exp(-ar * ar);
			energy = poly5(ar) * e * rr;
			force = ((energy + my2AlphaPI * e) * rr * rr);
#else	  
      Real a = erfc(myAlpha*r)/r;
      energy = qq*a;
      force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*distSquared))*rDistSquared;
#endif
		}

		static void accumulateEnergy(ScalarStructure* energies, Real energy)
		{
			(*energies)[ScalarStructure::COULOMB] += energy;
		}

		static Real getEnergy(const ScalarStructure* energies)
		{
			return (*energies)[ScalarStructure::COULOMB];
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

		void getParameters(std::vector<Parameter>&) const;
		static CoulombEwaldRealForce make(std::string&, const std::vector<Value>&);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		Real myAlpha;
		Real myAlphaSquared;
		Real my2AlphaPI;
	};

	//______________________________________________________________________ INLINES
}
#endif /* COULOMBEWALDREALFORCE_H */
