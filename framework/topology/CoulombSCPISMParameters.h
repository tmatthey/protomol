/* -*- c++ -*- */
#ifndef COULOMBSCPISM_H
#define COULOMBSCPISM_H

#include "Real.h"
#include "AtomType.h"

namespace ProtoMol
{
	//_________________________________________________________________ CoulombSCPISMParameters

	/// The Coulomb SCPISM for an atom type
	struct CoulombSCPISMParameters
	{
	public:

		CoulombSCPISMParameters(Real _al = 0.0, Real _hb = 0.5, Real _r = 0.5, Real _sq = 0.0, Real _rc = 0.37, Real _ga = 0.52, Hbonded _is = NO)
			: alpha_i(_al), hbond_factor(_hb), R_iw(_r), sqrt_alpha_i(_sq), r_cov(_rc), gamma_i(_ga),
			  isHbonded(_is)
		{
		}

		void set(Real _al = 0.0, Real _hb = 0.5, Real _r = 0.5, Real _sq = 0.0, Real _rc = 0.37, Real _ga = 0.0052, Hbonded _is = NO,
		         Real _A = 0.0, Real _B = 0.0, Real _C = 0.0, Real _Rvdw = 0.0)
		{
			alpha_i = _al;
			hbond_factor = _hb;
			R_iw = _r;
			sqrt_alpha_i = _sq;
			r_cov = _rc;
			gamma_i = _ga;
			isHbonded = _is;
			A_i = _A;
			B_i = _B;
			C_i = _C;
			R_vdw = _Rvdw;
		}

		Real alpha_i; // Alpha_i controls slope of D(r) around atom type i
		Real hbond_factor; // hbond_factor (Polar H) * hbond_factor (PA) controls 
		// increment of Born radius to correct h bonding strength
		Real R_iw; // extension of Born radius R_iw to obtain R_ip
		Real sqrt_alpha_i; // Square root of alpha_i
		Real r_cov; // Covalent radius (only values for C,N,O,S,H)
		Real gamma_i; // hydrophobic energy term
		Hbonded isHbonded; // whether atom is involved in H-bonding
		Real A_i; // A_i in (2) on SCPISM.DOC
		Real B_i; // B_i in (2) on SCPISM.DOC
		Real C_i; // C_i in (2) on SCPISM.DOC
		Real R_vdw; // R_i,vdw in (2) on SCPISM.DOC
	};
}
#endif /* not COULOMBSCPISM_H */
