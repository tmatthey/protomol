/* -*- c++ -*- */
#ifndef ISGIDEALGASLENNARDJONESFORCE_H
#define ISGIDEALGASLENNARDJONESFORCE_H

#include "GenericTopology.h"
#include "LennardJonesParameters.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ iSGIdealGasLennardJonesForce
  class iSGIdealGasLennardJonesForce {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGIdealGasLennardJonesForce() {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ISGLennardJonesForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force, 
		    Real &DeltaMu,
		    Real DistSquared,
		    Real rDistSquared,
		    const Vector3D &, 
		    const GenericTopology* topo, 
		    int atom1, int atom2,
		    ExclusionClass excl) const {

      // Use the molecule types to determine if the atoms
      // are being transformed or not
      int M1 = topo->atoms[atom1].molecule;
      int M2 = topo->atoms[atom2].molecule;

      // transforming atom indicators
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage, atom_stage;
      myStage = OldStage = atom_stage = 0;
      
      // if lambda for either molecule is nonzero then the molecule is being transformed
      if (topo->molecules[M1].lambda != 0.0) {atom1_scaled = true;}
      if (topo->molecules[M2].lambda != 0.0) {atom2_scaled = true;}

      // skip if the molecules are different (intermolecular case)
      if (M1 != M2) return;
            
      // Fast LJ computations
      Real r2   = rDistSquared;
      Real r6   = r2*r2*r2;
      Real r12, rij_6;

      // soft-core interaction variables
      Real Vr6_old, Vr6_new, Vr12_old, Vr12_new;
      Real Lambda, alpha;

      // energy and size parameter variables
      LennardJonesParameters params, params_old, params_new;
      Real A_old, A_new, B_old, B_new;

      // LJ parameter bank indices
      unsigned int i, j;

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
        // then neither atom is currently being transformed
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
        else return;
      }
      
      //************************************************************************
      // Compute the appropriate LJ interaction
      DeltaMu = 0.0;
      switch (Choice) {
        //----------------------------------------------------------
      case 0:
	// interaction between two non-transforming atoms

	// get the identities of the two atoms
	// i, j and are used to retrieve the correct set of
	// parameters from the LJ parameter bank
	// to access the proper LJ bank, atomI must always be the
	// atom with the LOWER atom_type number -- TIM
	if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
	  i = topo->molecules[M1].type;
	  j = topo->molecules[M2].type;}
	else {
	  i = topo->molecules[M2].type;
	  j = topo->molecules[M1].type;}

	
	// get the LJ parameters for this pair
	params = topo->isgLJParms(i,j,topo->atoms[atom1].type,topo->atoms[atom2].type);
	A_old = (excl != EXCLUSION_MODIFIED ? params.A : params.A14);
	B_old = (excl != EXCLUSION_MODIFIED ? params.B : params.B14);
  
	// if any of these atoms are dummy atoms (epsilon = 0), then skip
	if (A_old == 0.0) break;

	// attractive LJ term    
	Vr6_old = B_old * r6;
        
	// repulsive LJ term
	r12 = r6 * r6;
	Vr12_old = A_old * r12;
	
	// energy and force
	energy = Vr12_old - Vr6_old;
	force = 12.0*Vr12_old*r2 - 6.0*r2*Vr6_old;
	break;
        //----------------------------------------------------------
      case 2:  
	// intramolecular interaction on a transforming molecule

        // Get the lambda factor for each atom
        Lambda = topo->molecules[M1].lambda;

        // Get the alphaLJ parameter for this interaction
        Real tempAlpha1 = topo->atoms[atom1].alphaLJ;
        Real tempAlpha2 = topo->atoms[atom2].alphaLJ;
        if (tempAlpha1 > tempAlpha2) alpha = tempAlpha1;
        else alpha = tempAlpha2;

        // retrieve the A-A and B-B energy and size parameters
        // j is used to retrieve the correct set of
        // parameters from the LJ parameter bank        
        j = topo->molecules[M1].type;
 
        // get the LJ parameters for this pair when the atoms are in the old state
        params_old = topo->isgLJParms(j,j,topo->atoms[atom1].type,topo->atoms[atom2].type);

        // determine the proper index # into to the bank when the atoms are in the new state
        j = topo->molecules[M1].newtype;

        // get the LJ parameters for this pair when the atoms are in the new state
        params_new = topo->isgLJParms(j,j,topo->atoms[atom1].type,topo->atoms[atom2].type);
        A_old = (excl != EXCLUSION_MODIFIED ? params_old.A : params_old.A14);
        B_old = (excl != EXCLUSION_MODIFIED ? params_old.B : params_old.B14);
        A_new = (excl != EXCLUSION_MODIFIED ? params_new.A : params_new.A14);
        B_new = (excl != EXCLUSION_MODIFIED ? params_new.B : params_new.B14);

        // if both of these atoms are dummy atoms (epsilon = 0), then skip
        if (A_old == 0.0 && A_new == 0.0) break;
        
        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
           
        // LJ soft-core energy scaling factors
        Real LambdaMu = 2.0 * alpha * (Lambda - OldStage) * (myStage - Lambda);
        Real Scale_old = alpha * (Lambda - OldStage) * (Lambda - OldStage);
        Real Scale_new = alpha * (myStage - Lambda) * (myStage - Lambda);

        // soft-core interaction distances
        // if the transforming atom is a dummy atom in one of its identities then
        // we must make sure not to divide by zero
        Real SoftDist_old, SoftDist_new;        
        rij_6 = DistSquared * DistSquared * DistSquared;
        if (A_old != 0.0) SoftDist_old = rij_6 + (A_old / B_old) * Scale_old;
        else SoftDist_old = 1.0;
        if (A_new != 0.0) SoftDist_new = rij_6 + (A_new / B_new) * Scale_new;
        else SoftDist_new = 1.0;

        // attractive LJ terms
        Vr6_old = B_old / SoftDist_old;
        Vr6_new = B_new / SoftDist_new;
 
        // repulsive LJ terms
        Vr12_old = A_old / (SoftDist_old * SoftDist_old);
        Vr12_new = A_new / (SoftDist_new * SoftDist_new);
                        
        // get the transformation stage #'s for these atoms
        int atom1_stage = topo->atoms[atom1].stageNumber;
        int atom2_stage = topo->atoms[atom2].stageNumber;
        int max_stage = max(atom1_stage,atom2_stage);

        // if the transforming atom is a dummy atom in one of its identities then
        // we must make sure not to divide by zer        
        if (A_old == 0.0) {
          A_old = 1.0;
          B_old = 1.0;}
        if (A_new == 0.0) {
          A_new = 1.0;
          B_new = 1.0;}

        // if the current stage is not the same as the largest atom stage number then 
        // we do not compute any deltaMu, and if the current stage is less than the largest
        // atom stage number then we only compute interactions in the old identity, if the
        // current stage is larger than both atom stage numbers then we use only the new identity       
        if (myStage == max_stage) {
          // energy, force, and chemical potential difference
          energy = (myStage - Lambda) * (Vr12_old - Vr6_old) + (Lambda - OldStage) * (Vr12_new - Vr6_new);
          force = 6.0 * r2 * rij_6 * ((myStage - Lambda) * Vr12_old * (2.0 * Vr6_old/B_old - B_old/A_old)
            + (Lambda - OldStage) * Vr12_new * (2.0 * Vr6_new/B_new - B_new/A_new));
          DeltaMu = (Vr12_new - Vr6_new) - (Vr12_old - Vr6_old)
            + LambdaMu * (Vr12_new * (2.0 * Vr6_new*(A_new/(B_new*B_new)) - 1.0)
            - Vr12_old * (2.0 * Vr6_old*(A_old/(B_old*B_old)) - 1.0));
        }
        else if (myStage > max_stage) {
          // energy, force, and chemical potential difference (new identity)
          energy = (Lambda - OldStage) * (Vr12_new - Vr6_new);
          force = 6.0 * r2 * rij_6 * (Lambda - OldStage) * Vr12_new * (2.0 * Vr6_new/B_new - B_new/A_new);
        }
        else {
          // energy, force, and chemical potential difference (old identity)
           energy = (myStage - Lambda) * (Vr12_old - Vr6_old);
          force = 6.0 * r2 * rij_6 * (myStage - Lambda) * Vr12_old * (2.0 * Vr6_old/B_old - B_old/A_old);
        }    
        break;
	//************************************************************************
      } // end switch structure

    }
 

    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
      (*energies)[ScalarStructure::LENNARDJONES_DELTAMU] += deltaMu;
    }
    static Real getEnergy(const ScalarStructure* energies) {
      return  (*energies)[ScalarStructure::LENNARDJONES];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter>&) const{}
    static iSGIdealGasLennardJonesForce make(std::string& , std::vector<Value>) {
      return iSGIdealGasLennardJonesForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;    

  private:
  };
  //______________________________________________________________________ INLINES
}

#endif /* ISGIDEALGASLENNARDJONESFORCE_H */
