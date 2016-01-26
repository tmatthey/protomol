/* -*- c++ -*- */
#ifndef ISGLENNARDJONESFORCE_H
#define ISGLENNARDJONESFORCE_H

#include "GenericTopology.h"
#include "LennardJonesParameters.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include <string>

// uncomment for debugging purposes
//#define DEBUG_NORMAL_LJ
//#define DEBUG_SOFTCORE
//#define DEBUG_LJ_DISTANCE_CHECK

namespace ProtoMol {

  //_________________________________________________________________ ISGLennardJonesForce
  class iSGLennardJonesForce {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGLennardJonesForce() {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGLennardJonesForce
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
#ifdef DEBUG_LJ_DISTANCE_CHECK
      if(rDistSquared > 1.0/0.25)
	Report::report << Report::warning << "iSGLennardJonesForce::operator(): atom "<<atom1<<" and "<<atom2<<" get closer"
	       << " than 0.5 angstroms (="<<sqrt(1.0/rDistSquared)<<")." << Report::endr;   
#endif
 
      // Fast LJ computations
      Real r2   = rDistSquared;
      Real r6   = r2*r2*r2;
      Real r12, rij_6;

      // soft-core interaction variables
      Real Vr6_old = 0.0;
      Real Vr6_new = 0.0;
      Real Vr12_old = 0.0;
      Real Vr12_new = 0.0;
      Real Scale_new, Scale_old, Lambda, alpha;
      Real SoftDist_old, SoftDist_new, LambdaMu;

      // energy and size parameter variables
      LennardJonesParameters params, params_old, params_new;
      Real A_old = 0.0;
      Real A_new = 0.0;
      Real B_old = 0.0;
      Real B_new = 0.0;

      unsigned int type1_old, type1_new;

      // LJ parameter bank indices
      unsigned int type2, i, j;

      // transforming atom indicators
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage, atom_stage;
      myStage = OldStage = atom_stage = 0;

      // Use the molecule types to determine if the atoms
      // are being transformed or not
      int M1 = topo->atoms[atom1].molecule;
      int M2 = topo->atoms[atom2].molecule;

      // if lambda for either molecule is nonzero then the molecule is being transformed
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
        else Choice = 1;
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
        // atom with the LOWER atom_type number since that is how
        // the bank was constructed in the buildISGTopology function
        if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
          i = topo->molecules[M1].type;
          j = topo->molecules[M2].type;}
        else {
          i = topo->molecules[M2].type;
          j = topo->molecules[M1].type;}

        // check to see if atom1 may have already been successfully transformed in a previous stage  
        if (atom1_scaled) {
          // this atom has already been transformed if myStage > atom_stage
          if (myStage > atom_stage) {
            if (topo->atoms[atom1].type <= topo->atoms[atom2].type) i = topo->molecules[M1].newtype;
            else j = topo->molecules[M1].newtype;
          }
        } // end if statement
        // check to see if atom2 may have already been successfully transformed in a previous stage
        if (atom2_scaled) {
          // this atom has already been transformed if myStage > atom_stage
          if (myStage > atom_stage) {
            if (topo->atoms[atom1].type <= topo->atoms[atom2].type) j = topo->molecules[M2].newtype;
            else i = topo->molecules[M2].newtype;
          }
        } // end if statement


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

#ifdef DEBUG_NORMAL_LJ
    	Report::report.precision(6);
	Report::report << "i = " << atom1+1 << ", j = " << atom2+1 << ", " << i << ", " << j
                       << ", eps = " << B_old*B_old / (4.0*A_old) << ", energy = "
                       << energy << Report::endr;
#endif
	break;
        //----------------------------------------------------------
      case 1:  
	// interaction between non-transforming atom and a transforming atom
	// 'old' indicates OldType, 'new' indicates NewType

	if (atom1_scaled) {
	  // get the lambda factor
	  Lambda = topo->molecules[M1].lambda;

	  // get the identities of the two atoms
	  type1_old = topo->molecules[M1].type;
	  type1_new = topo->molecules[M1].newtype;
	  type2 = topo->molecules[M2].type;

          // get the alphaLJ parameter for this atom
          alpha = topo->atoms[atom1].alphaLJ;
	}
	else {
	  // get the lambda factor
	  Lambda = topo->molecules[M2].lambda;

	  // get the identities of the two atoms
	  type1_old = topo->molecules[M2].type;
	  type1_new = topo->molecules[M2].newtype;
	  type2 = topo->molecules[M1].type;

          // get the alphaLJ parameter for this atom
          alpha = topo->atoms[atom2].alphaLJ;
	}

        // i, j and are used to retrieve the correct set of
        // parameters from the LJ parameter bank
        // if atom 1 is the transforming atom...
        if (atom1_scaled) {

          // to access the proper LJ bank, atomI must always be the
          // atom with the LOWER atom_type number -- TIM
          if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
            i = type1_old;
            j = type2;}
          else {
            i = type2;
            j = type1_old;
          }
        }
        // if atom2 is the transforming atom then the procedure
        // below is the opposite of the procedure above.  AtomI must still
        // be the atom with the LOWER atom_type number, but remember that
        // if atom2 is the transforming atom then the variable type2 is equal
        // to the molecule_type of atom 1!!!
        else {

          // if atom1 has a lower atom_type, then i = type2 since now type2
          // is equal to the molecule_type of atom1
          if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
            i = type2;
            j = type1_old;}
          // if atom2 has the lower atom_type, then i = type1A
          else {
            i = type1_old;
            j = type2;
          }
        }

        
        // compute the A-B and B-B epsilon and sigmas
        // retrieve the A & B parameters for these interaction types
        // get the LJ parameters for this pair when the transforming atom is in the old state
        params_old = topo->isgLJParms(i,j,topo->atoms[atom1].type,topo->atoms[atom2].type);

        // determine the new i and j when the transforming atom is in the new state
        if (atom1_scaled) {
          if (topo->atoms[atom1].type <= topo->atoms[atom2].type) i = type1_new;
          else j = type1_new;}
        else {
          if (topo->atoms[atom1].type <= topo->atoms[atom2].type) j = type1_new;
          else i = type1_new;
        }

        // get the LJ parameters for this pair when the transforming atom is in the new state
        params_new = topo->isgLJParms(i,j,topo->atoms[atom1].type,topo->atoms[atom2].type);
        A_old = (params_old.A);
        B_old = (params_old.B);
        A_new = (params_new.A);
        B_new = (params_new.B);

        // if the type2 atom is a dummy atom then skip (epsilon = 0)
        if (A_old == 0.0 && A_new == 0.0) break;

        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
           
        // LJ soft-core energy scaling factors
        LambdaMu = 2.0 * alpha * (Lambda - OldStage) * (myStage - Lambda);
        Scale_old = alpha * (Lambda - OldStage) * (Lambda - OldStage);
        Scale_new = alpha * (myStage - Lambda) * (myStage - Lambda);

        // soft-core interaction distances
        // if the transforming atom is a dummy atom in one of its identities then
        // we must make sure not to divide by zero        
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

        // energy, force, and chemical potential difference
        // if the transforming atom is a dummy atom in one of its identities then
        // we must make sure not to divide by zero
        if (A_old == 0.0) {
          A_old = 1.0;
          B_old = 1.0;}
        if (A_new == 0.0) {
          A_new = 1.0;
          B_new = 1.0;}
        energy = (myStage - Lambda) * (Vr12_old - Vr6_old) + (Lambda - OldStage) * (Vr12_new - Vr6_new);
        force = 6.0 * r2 / r6 * ((myStage - Lambda) * Vr12_old * (2.0 * Vr6_old/B_old - B_old/A_old)
          + (Lambda - OldStage) * Vr12_new * (2.0 * Vr6_new/B_new - B_new/A_new));
        DeltaMu = (Vr12_new - Vr6_new) - (Vr12_old - Vr6_old)
         + LambdaMu * (Vr12_new * (2.0 * Vr6_new*(A_new/(B_new*B_new)) - 1.0)
         - Vr12_old * (2.0 * Vr6_old*(A_old/(B_old*B_old)) - 1.0));

#ifdef DEBUG_SOFTCORE
        Report::report.precision(6);
        Report::report << Report::warning << "i = " << atom1+1 << ", j = " << atom2+1 << ", " << type1_old << "->" << type1_new
                       << ", " << type2 << ", rij = " << 1.0/sqrt(r2) << ", Energy = " << energy
                       << ", force = " << force << ", DeltaMu = " << DeltaMu << Report::endr;
#endif //DEBUG_SOFTCORE
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
        LambdaMu = 2.0 * alpha * (Lambda - OldStage) * (myStage - Lambda);
        Scale_old = alpha * (Lambda - OldStage) * (Lambda - OldStage);
        Scale_new = alpha * (myStage - Lambda) * (myStage - Lambda);

        // soft-core interaction distances
        // if the transforming atom is a dummy atom in one of its identities then
        // we must make sure not to divide by zero        
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
    } // end operator()
 

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
    static iSGLennardJonesForce make(std::string& , const std::vector<Value>&) {
      return iSGLennardJonesForce();
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

#endif /* ISGLENNARDJONESFORCE_H */
