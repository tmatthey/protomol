#include "ReducedHessCoulombSCPISM.h"
#include "GenericTopology.h"


namespace ProtoMol
{
	//_________________________________________________________________ ReducedHessCoulombSCPISM
	Matrix3by3 ReducedHessCoulombSCPISM::operator()(const Real rawEnergy,
	                                                const Real rawForce,
	                                                Real a,
	                                                Real /*rDistSquared*/,
	                                                const Vector3D& rij,
	                                                const GenericTopology* topo,
	                                                int atom1, int atom2,
	                                                const Real switchingValue,
	                                                const Real switchingDeriv,
	                                                const Matrix3by3& switchingHess,
	                                                ExclusionClass excl,
	                                                Real alpha_ij,
	                                                Real es) const
	{
		es = 80;

		Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
		Real scaledCharge_j = topo->atoms[atom2].scaledCharge;
		Real coulombScalingFactor = 1;
		if (excl == EXCLUSION_MODIFIED)
		{
			coulombScalingFactor = topo->coulombScalingFactor;
		}

		Real na = sqrt(a);

		Real tm1 = (coulombScalingFactor * scaledCharge_i * scaledCharge_j) * (-1);

		Real k = (es - 1) * 0.5;
		Real Dsr = (1 + es) / (1 + k * exp(-alpha_ij * na)) - 1;
		Real Ddsr_dr = (alpha_ij / (1 + es)) * (1 + Dsr) * (es - Dsr); //first derivative of Dsr
		Real D2dsr_dr2 = (alpha_ij / (1 + es)) * (es - 1 - 2 * Dsr) * Ddsr_dr; //second derivative of Dsr

		Real P1r = (1 / (na * na * na * Dsr)) + (1 / (na * na * Dsr * Dsr)) * Ddsr_dr;

		//Real P2r = (1/(na*na*Dsr*Dsr))*D2dsr_dr2 - (1/(na*na*na*na*Dsr*Dsr*Dsr*Dsr))*Ddsr_dr*(2*na*Dsr*Dsr+2*na*na*Dsr*Ddsr_dr) - (1/(na*na*na*na*na*na*Dsr*Dsr))*(3*na*na*Dsr+na*na*na*Ddsr_dr);

		Real P2r = (1 / (na * na * na * Dsr * Dsr)) * D2dsr_dr2 - (2 / (na * na * na * na * Dsr * Dsr)) * Ddsr_dr - (2 / (na * na * na * Dsr * Dsr * Dsr)) * (Ddsr_dr * Ddsr_dr) - (3 / (na * na * na * na * na * Dsr)) - (1 / (na * na * na * na * Dsr * Dsr)) * Ddsr_dr;


		Matrix3by3 vec_rij_ij(rij, rij);// outer products of vectors r_ij 


		Matrix3by3 I(1, 0, 0, 0, 1, 0, 0, 0, 1); // now I is identity matrix

		//Matrix3by3 H((I + vec_rij_ij * (3/a))* (tm1 * switchingValue)); 
		Matrix3by3 H((I * P1r + vec_rij_ij * P2r) * (tm1 * switchingValue)); //
		// now done with the electrostatic (Coulomb) part.


		H += (switchingHess * rawEnergy);
		H -= vec_rij_ij * (2 * rawForce * switchingDeriv);

		return H;
	}
}
