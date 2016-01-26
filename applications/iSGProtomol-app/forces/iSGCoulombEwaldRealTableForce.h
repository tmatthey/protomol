/* -*- c++ -*- */
#ifndef ISGCOULOMBEWALDREALTABLEFORCE_H
#define ISGCOULOMBEWALDREALTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ CoulombEwaldRealTableForce
  class iSGCoulombEwaldRealTableForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Replacement potential of the real part for Ewald or PME.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    enum {DEFAULT_TABLE_SIZE=20000};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGCoulombEwaldRealTableForce();
    iSGCoulombEwaldRealTableForce(Real a, Real rc, unsigned int n);
    iSGCoulombEwaldRealTableForce(iSGCoulombEwaldRealTableForce const& other);
    iSGCoulombEwaldRealTableForce& operator=(iSGCoulombEwaldRealTableForce const& other);
    ~iSGCoulombEwaldRealTableForce();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombEwaldRealTableForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force,
                    Real &deltaMu, 
		    Real distSquared,
		    Real /*rDistSquared*/,
		    const Vector3D &, 
		    const GenericTopology* topo, 
		    int atom1, int atom2,
		    ExclusionClass excl) const {


      // needed variables
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage, atom_stage;
      myStage = OldStage = atom_stage = 0;
            
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
      else if ( (atom1_scaled) && (atom2_scaled) ) {             

        // determine the current transformation stage of the molecule
        myStage = static_cast<int>(floor(topo->molecules[M1].lambda)) + 1;

        // get the transformation stage #'s for these atoms
        int atom1_stage = topo->atoms[atom1].stageNumber;
        int atom2_stage = topo->atoms[atom2].stageNumber;

        // if the current stage is not the same as either atom1_stage or atom2_stage
        // the neither atom is currently being transformed
        if (myStage != atom1_stage && myStage != atom2_stage) Choice = 0;
        else Choice = 2;
      }
      // one of the atoms is being transformed
      else {
        // determine the current transformation stage of the molecule
        // and get the transformation stage # for this atoms
        if (atom1_scaled) {
          myStage = static_cast<int>(floor(topo->molecules[M1].lambda)) + 1;
          atom_stage = topo->atoms[atom1].stageNumber;}
        else {
          myStage = static_cast<int>(floor(topo->molecules[M2].lambda)) + 1;
          atom_stage = topo->atoms[atom2].stageNumber;}
        
        // if the current stage is not the same as atom_stage
        // then neither atom is currently being transformed
        if (myStage != atom_stage) Choice = 0;
        else Choice = 1;
      }
      
      // the scaled charges on each atom
      Real q1 = topo->atoms[atom1].scaledCharge;
      Real q2 = topo->atoms[atom2].scaledCharge;

      // Compute the appropriate electrostatic interaction
      Real qq, q_Dq, myLambda;
      Real qq_correction_add, qq_correction_minus, q_Dq_correction_add, q_Dq_correction_minus;
      Real e_correction, myCorrection, myDMUcorrection;
      qq = q_Dq = myLambda = 0.0;
      qq_correction_add = qq_correction_minus = q_Dq_correction_add = q_Dq_correction_minus = 0.0;
      e_correction = myCorrection = myDMUcorrection = 0.0;
 
      // determine the ERFC table bin # for this interatomic distance
      Real n = sqrt(sqrt(distSquared))*myFac;
      int index = (int)n;
      Real b = n-index;
      //Report::report << index <<","<<sqrt(distSquared);
      index <<=1;
      //Report::report <<","<< index <<","<<myN<<","<<(unsigned long long)myTable<< Report::endr;
      Real a = 1.0-b;

      // look up the precomputed values of ERFC to get the energy, force, and chemical potential difference
      Real myErfc = myTable[index  ]*a+myTable[index+2]*b;
      Real myErfcForce = myTable[index+1]*a+myTable[index+3]*b;
      
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

	// compute the energy and force
	energy = myErfc * qq;
	force  = myErfcForce * qq;
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

	// compute the energy, force, and deltaMu
	energy = myErfc * qq;
	force  = myErfcForce * qq;
	deltaMu = myErfc * q_Dq;
        break;
        //----------------------------------------------------------
      case 2:
	// intramolecular interaction on a transforming molecule
	// get the value of lambda for this molecule
	myLambda = topo->molecules[M1].lambda;

        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
        //if (myLambda > myStage) myLambda = myStage;
        //else if (myLambda < OldStage) myLambda = OldStage;
                     
	// Energy
	// The interaction energy is the sum of the interaction in the
	// old state scaled by (1 - lambda) plus the interaction in the new state
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
          
	// determine the ERF table bin # for this interatomic distance
	int ERFindex = index/2;

	// compute the energy, force, and deltaMu
	energy = myErfc * qq;
	force  = myErfcForce * qq;
	deltaMu = myErfc * q_Dq;
          
	// look up the precomputed values of ERF to get the energy correction term
	e_correction = myERFTable[ERFindex]*a + myERFTable[ERFindex+1]*b;
          
	// add in the reciprocal space correction terms
	myCorrection = e_correction * (qq_correction_add - qq_correction_minus);
	myDMUcorrection = e_correction * (q_Dq_correction_add - q_Dq_correction_minus);
	energy += myCorrection;
	deltaMu += myDMUcorrection;
        break;
        //----------------------------------------------------------
      } // end switch structure
    } // end operator()

    //_________________________________________________________________ accumulateEnergy
    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::COULOMB] += energy;
      (*energies)[ScalarStructure::COULOMB_DELTAMU] += deltaMu;
    }

    //_________________________________________________________________ getEnergy
    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 3;}
    void getParameters(std::vector<Parameter>&) const;
    static iSGCoulombEwaldRealTableForce make(std::string&, const std::vector<Value>&);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;      
  private:
    Real* myTable;
    Real* myERFTable;
    Real myAlpha;
    Real myRc;
    Real myAlphaSquared;
    Real my2AlphaPI;
    Real myFac;
    unsigned int myN;
  };

  //______________________________________________________________________ INLINES
}
#endif /* ISGCOULOMBEWALDREALTABLEFORCE_H */
