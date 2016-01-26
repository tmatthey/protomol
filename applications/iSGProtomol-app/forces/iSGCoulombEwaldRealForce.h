/* -*- c++ -*- */
#ifndef ISGCOULOMBEWALDREALFORCE_H
#define ISGCOULOMBEWALDREALFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include <string>

//#define USE_COULOMBEWALDREAL_EXACT_ERF
namespace ProtoMol {

  //_________________________________________________________________ iSGCoulombEwaldRealForce
  class iSGCoulombEwaldRealForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Replacement potential of the real part for iSGEwald or iSGPME.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGCoulombEwaldRealForce();
    iSGCoulombEwaldRealForce(Real a);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombEwaldRealForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force,
                    Real &deltaMu,
		    Real distSquared,
		    Real 
#ifdef USE_COULOMBEWALDREAL_EXACT_ERF
		    rDistSquared
#endif
		    ,
		    const Vector3D &, 
		    const GenericTopology* topo, 
		    int atom1, int atom2,
		    ExclusionClass excl) const {


      // needed variables
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage /*, atom_stage*/;
      myStage = OldStage /*= atom_stage*/ = 0;
      
      // Use the molecule types to determine if the atoms
      // are being transformed or not
      int M1 = topo->atoms[atom1].molecule;
      int M2 = topo->atoms[atom2].molecule;
      if (topo->molecules[M1].lambda != 0.0) {atom1_scaled = true;}
      if (topo->molecules[M2].lambda != 0.0) {atom2_scaled = true;}

      // integer used to select which interaction type to compute   
      int Choice;

      // determine the value of Choice
      // neither atom is being transformed
      if ( !(atom1_scaled) && !(atom2_scaled) ) {Choice = 0;}
      // intramolecular interaction on a transforming molecule
      else if ( (atom1_scaled) && (atom2_scaled) ) Choice = 2;

      // one of the atoms is being transformed
      else  Choice = 1;
      
      // the scaled charges on each atom
      Real q1 = topo->atoms[atom1].scaledCharge;
      Real q2 = topo->atoms[atom2].scaledCharge;

      // Compute the appropriate electrostatic interaction
      Real qq, rr, ar, e, q_Dq, dmu, myLambda;
      Real qq_correction_add, qq_correction_minus, q_Dq_correction_add, q_Dq_correction_minus;
      Real e_correction, myCorrection, myDMUcorrection;
      qq = rr = ar = e = q_Dq = dmu = myLambda = 0.0;
      qq_correction_add = qq_correction_minus = q_Dq_correction_add = q_Dq_correction_minus = 0.0;
      e_correction = myCorrection = myDMUcorrection = 0.0;
#ifdef USE_COULOMBEWALDREAL_EXACT_ERF
      Real a = 0.0;
#endif         
      // compute the separation distance
      Real r = sqrt(distSquared);
            
      switch (Choice) {
        //----------------------------------------------------------
      case 0:
	// interaction between two non-transforming atoms
	// compute the raw interaction energy and force
              
	// skip if either of the charges is zero
	if (q1 == 0.0 || q2 == 0.0) break;

	// product of charges
	qq = q1*q2;
              
	// scale any 1-4 interactions if desired
	if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED)
          qq *= topo->coulombScalingFactor;
  
	// multiply qq by the complementary error function
	// Approximation Abramowitz & Stegun p299.
	// Energy and force
#ifndef USE_COULOMBEWALDREAL_EXACT_ERF
	rr = 1.0/r;
	ar = myAlpha*r;
	e = qq*exp(-ar*ar);
	energy = poly5(ar)*e*rr;
	force = ((energy+my2AlphaPI*e)*rr*rr);
#else     
	a = erfc(myAlpha*r)/r;
	energy = qq*a;
	force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	break;
	//----------------------------------------------------------
      case 1:
	// interaction between non-transforming atom and a transforming atom
	// compute the raw interaction energy, force, and chemical potential difference
                 
	// skip if either of the charges is zero
	if (q1 == 0.0 || q2 == 0.0) break;

	// product of charges
	qq = q1*q2;

	// derivative of qq with respect to lambda
	q_Dq = (q1 * topo->atoms[atom2].deltaQ
		+ q2 * topo->atoms[atom1].deltaQ);
              
	// multiply qq and q_Dq by the complementary error function
	// Approximation Abramowitz & Stegun p299.
	// Energy
#ifndef USE_COULOMBEWALDREAL_EXACT_ERF
	rr = 1.0/r;
	ar = myAlpha*r;
	e = qq*exp(-ar*ar);
	dmu = q_Dq*exp(-ar*ar);
	energy = poly5(ar)*e*rr;
	deltaMu = poly5(ar)*dmu*rr;
	force = ((energy+my2AlphaPI*e)*rr*rr);
#else     
	a = erfc(myAlpha*r)/r;
	energy = qq*a;
	deltaMu = q_Dq*a;
	force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
#endif
	break;
	//----------------------------------------------------------
      case 2:
	// intramolecular interaction on a transforming molecule
	// get the value of lambda for this molecule
	myLambda = topo->molecules[M1].lambda;

        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
                            
	// Energy
	// The interaction energy is the sum of the interaction in the
	// old state scaled by (myStage - lambda) plus the interaction in the new state
	// scaled by lambda
              
	// product of charges for each identity
	qq = (myStage - myLambda) * topo->atoms[atom1].Qold * topo->atoms[atom2].Qold
	  + (myLambda - OldStage) * topo->atoms[atom1].Qnew * topo->atoms[atom2].Qnew;
           
	// The chemical potential difference is the derivative of the
	// energy above (qq) with respect to lambda
	q_Dq = topo->atoms[atom1].Qnew * topo->atoms[atom2].Qnew
	  - topo->atoms[atom1].Qold * topo->atoms[atom2].Qold;
              
	// correction terms for intramolecular self energy and chemical potential difference
	// these correction terms are needed because the reciprocal space term contains an
	// intramolecular self term that is proportional to atoms[i].scaledCharge * atoms[j].scaledCharge,
	// which is NOT the correct intramolecular self interaction (see qq and q_Dq terms above).
	// Correcting the reciprocal intramolecular interaction amounts to the addition of a factor
	// of qq to the energy and q_Dq to the chemical potential difference, and subtraction of a
	// factor of atoms[i].scaledCharge * atoms[j].scaledCharge from the energy and
	// subtraction of (atoms[i].scaledCharge * atoms[j].deltaQ + atoms[j].scaledCharge * atoms[i].deltaQ)
	// from the chemical potential difference
	qq_correction_add = qq;
	q_Dq_correction_add = q_Dq;
	qq_correction_minus = 0.0;
	q_Dq_correction_minus = 0.0;
              
	// scale any 1-4 interactions if desired
	if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED) {
	  qq *= topo->coulombScalingFactor;
	  q_Dq *= topo->coulombScalingFactor;
	  qq_correction_add *= topo->coulombScalingFactor;
	  q_Dq_correction_add *= topo->coulombScalingFactor;
	}
	else {
	  // this is an unmodified intramolecular interaction, so subtract the total
	  // reciprocal space energy and chemical potential difference
	  qq_correction_minus = q1 * topo->atoms[atom2].scaledCharge;
	  q_Dq_correction_minus = q1 * topo->atoms[atom2].deltaQ + q2 * topo->atoms[atom1].deltaQ;
	}
              
	// multiply qq and q_Dq by the complementary error function
	// Approximation Abramowitz & Stegun p299.
	// Energy
#ifndef USE_COULOMBEWALDREAL_EXACT_ERF
	rr = 1.0/r;
	ar = myAlpha*r;
	e = qq*exp(-ar*ar);
	dmu = q_Dq*exp(-ar*ar);
	energy = poly5(ar)*e*rr;
	deltaMu = poly5(ar)*dmu*rr;
	force = ((energy+my2AlphaPI*e)*rr*rr);
	// Reciprocal space correction term gets
	// multiplied by the error function
	e_correction = erf(myAlpha*r)*rr;
#else     
	a = erfc(myAlpha*r)/r;
	energy = qq*a;
	deltaMu = q_Dq*a;
	force = qq*(a+my2AlphaPI*exp(-myAlphaSquared*rSquared))/rSquared;
	// Reciprocal space correction term gets
	// multiplied by the error function
	e_correction = erf(myAlpha*r)/r;
#endif
	// add in the reciprocal space correction terms
	myCorrection = e_correction * (qq_correction_add - qq_correction_minus);
	myDMUcorrection = e_correction * (q_Dq_correction_add - q_Dq_correction_minus);
	energy += myCorrection;
	deltaMu += myDMUcorrection;
	break;
	//----------------------------------------------------------
      } // end switch structure
    }

    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::COULOMB] += energy;
      (*energies)[ScalarStructure::COULOMB_DELTAMU] += deltaMu;
    }

    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 1;}
    void getParameters(std::vector<Parameter>&) const;
    static iSGCoulombEwaldRealForce make(std::string&, const std::vector<Value>&);
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
#endif /* ISGCOULOMBEWALDREALFORCE_H */
