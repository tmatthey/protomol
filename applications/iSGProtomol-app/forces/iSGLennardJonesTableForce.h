/* -*- c++ -*- */
#ifndef ISGLENNARDJONESTABLEFORCE_H
#define ISGLENNARDJONESTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ iSGLennardJonesTableForce
  class iSGLennardJonesTableForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Replacement potential of the real part for  or PME.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    enum {DEFAULT_TABLE_SIZE=20000};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGLennardJonesTableForce();
    iSGLennardJonesTableForce(Real rc, unsigned int n, Real myAlpha);
    iSGLennardJonesTableForce(iSGLennardJonesTableForce const& other);
    iSGLennardJonesTableForce& operator=(iSGLennardJonesTableForce const& other);
    ~iSGLennardJonesTableForce();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGLennardJonesTableForce
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


      // soft-core interaction variables
      Real Vr6A, Vr6B, Vr12A, Vr12B;
      Real ScaleA, ScaleB, Lambda;

      // energy and size parameter variables
      LennardJonesParameters params, params0, params1;
      Real A0, A1, B0, B1;

      // variables that are not needed for an ideal gas simulation
      // soft-core interaction variables
      Real SoftDistA, SoftDistB, sr6A, sr6B;
      Real LambdaMu;

      // energy and size parameter variables
      Real epsilonA, epsilonB, epsilon0, epsilon1;
      Real sig60, sig61;
      unsigned int type1A, type1B;

      // LJ parameter bank indices
      unsigned int type2, i, j;

      // transforming atom indicators
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage, atom_stage;
      myStage = OldStage = atom_stage = 0;
            
      // LJTable look-up variables
      Real n, a, b;
      int index;
      
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
      Real r2, r6;
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
        A0 = (excl != EXCLUSION_MODIFIED ? params.A : params.A14);
        B0 = (excl != EXCLUSION_MODIFIED ? params.B : params.B14); 
  
        // if any of these atoms are dummy atoms (epsilon = 0), then skip
        if (A0 == 0.0) break;

        // determine the DeltaR bin to use
        n = sqrt(distSquared)*myFac;
        index = (int)n;
        b = n-index;
        //Report::report << index <<","<<sqrt(distSquared);
        index <<=2;
        //Report::report <<","<< index <<","<<myN<<","<<(unsigned long long)myTable<< Report::endr;
        a = 1.0-b;

        // use the LJTable bin # to get the precomputed values of r^6 and r^12
        energy = (myTable[index  ]*a+myTable[index+4]*b)*A0 + (myTable[index+1]*a+myTable[index+5]*b)*B0;
        force  = (myTable[index+2]*a+myTable[index+6]*b)*A0 + (myTable[index+3]*a+myTable[index+7]*b)*B0;
        break;
        //----------------------------------------------------------
      case 1:  
        // interaction between non-transforming atom and a transforming atom
        // 'A' indicates OldType, 'B' indicates NewType

        if (atom1_scaled) {
          // get the lambda factor
          Lambda = topo->molecules[M1].lambda;

          // get the identities of the two atoms
          type1A = topo->molecules[M1].type;
          type1B = topo->molecules[M1].newtype;
          type2 = topo->molecules[M2].type;
        }
        else {
          // get the lambda factor
          Lambda = topo->molecules[M2].lambda;

          // get the identities of the two atoms
          type1A = topo->molecules[M2].type;
          type1B = topo->molecules[M2].newtype;
          type2 = topo->molecules[M1].type;
        }

        // i, j and are used to retrieve the correct set of
        // parameters from the LJ parameter bank
        // if atom 1 is the transforming atom...
        if (atom1_scaled) {

          // to access the proper LJ bank, atomI must always be the
          // atom with the LOWER atom_type number -- TIM
          if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
            i = type1A;
            j = type2;}
          else {
            i = type2;
            j = type1A;
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
            j = type1A;}
          // if atom2 has the lower atom_type, then i = type1A
          else {
            i = type1A;
            j = type2;
          }
        }
   
        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
        //if (Lambda > myStage) Lambda = myStage;
        //else if (Lambda < OldStage) Lambda = OldStage;
           
        // LJ soft-core energy scaling factors
        LambdaMu = alpha * (Lambda - OldStage) * (myStage - Lambda);
        ScaleA = alpha * (Lambda - OldStage) * (Lambda - OldStage);
        ScaleB = alpha * (myStage - Lambda) * (myStage - Lambda);
        
        // compute the A-B and B-B epsilon and sigmas
        // retrieve the A & B parameters for these interaction types
        // get the LJ parameters for this pair when the transforming atom is in the '0', or old state
        params0 = topo->isgLJParms(i,j,topo->atoms[atom1].type,topo->atoms[atom2].type);

        // determine the new i and j when the transforming atom is in the '1', or new state
        if (topo->atoms[atom1].type <= topo->atoms[atom2].type) {
          if (atom1_scaled) i = type1B;
          else j = type1B;
        }
        else {
          if (atom1_scaled) j = type1B;
          else i = type1B;
        }  

        // get the LJ parameters for this pair when the transforming atom is in the '1', or new state
        params1 = topo->isgLJParms(i,j,topo->atoms[atom1].type,topo->atoms[atom2].type);
        A0 = (params0.A);
        B0 = (params0.B); 
        A1 = (params1.A);
        B1 = (params1.B); 

        // if the type2 atom is a dummy atom then skip (epsilon = 0)
        if (A0 == 0.0 && A1 == 0.0) break;

        // Use A & B to compute epsilon and sigma
        // for the "0" state
        if (A0 != 0.0) {
          sig60 = A0 / B0;
          epsilon0 = B0 / (4.0 * sig60);}
        else {
          sig60 = 1.0;
          epsilon0 = 0.0;}
        // for the "1" state
        if (A1 != 0.0) {
          sig61 = A1 / B1;
          epsilon1 = B1 / (4.0 * sig61);}
        else {
          sig61 = 1.0;
          epsilon1 = 0.0;}

        // scale epsilon by the appropriate factor of lambda
        epsilonA = (myStage - Lambda) * epsilon0;
        epsilonB = (Lambda - OldStage) * epsilon1;

        // Fast LJ computations
        r2   = 1./distSquared;
        r6   = r2*r2*r2;
      
        // soft-core interaction distances
        sr6A = sig60 * r6;
        sr6B = sig61 * r6;
        SoftDistA = ScaleA + 1.0 / sr6A;
        SoftDistB = ScaleB + 1.0 / sr6B;

        // attractive LJ terms    
        Vr6A = 1.0 / SoftDistA;
        Vr6B = 1.0 / SoftDistB;
 
        // repulsive LJ terms
        Vr12A = Vr6A * Vr6A;
        Vr12B = Vr6B * Vr6B;
  
        // energy, force, and chemical potential difference
        energy = 4.0 * (epsilonA * (Vr12A - Vr6A) + epsilonB * (Vr12B - Vr6B));
        force = 24.0 * r2 * (epsilonA * Vr12A / sr6A * (2.0 * Vr6A - 1.0) +
                             epsilonB * Vr12B / sr6B * (2.0 * Vr6B - 1.0));
        deltaMu = 4.0 * (epsilon1 * (Vr12B - Vr6B) - epsilon0 * (Vr12A - Vr6A))
          + 8.0 * LambdaMu * (epsilon1 * Vr12B * (2.0 * Vr6B - 1.0) -
                              epsilon0 * Vr12A * (2.0 * Vr6A - 1.0));
#ifdef DEBUG_SOFTCORE
        Report::report.precision(6);
        Report::report << "i = " << atom1+1 << ", j = " << atom2+1 << ", " << type1A << "->" << type1B
                       << ", " << type2 << ", rij = " << 1.0/sqrt(r2) 
                       << ", eps0 = " << epsilon0 << ", eps1 = " << epsilon1 << ", Energy = " << energy
                       << ", force = " << force << ", DeltaMu = " << deltaMu << Report::endr;
#endif //DEBUG_SOFTCORE
        break;
        //----------------------------------------------------------
      case 2:  
        // intramolecular interaction on a transforming molecule

        // Get the lambda factor for each atom
        Lambda = topo->molecules[M1].lambda;

        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
        //if (Lambda > myStage) Lambda = myStage;
        //else if (Lambda < OldStage) Lambda = OldStage;
               
        // LJ soft-core energy scaling factors
        Real LambdaMu = alpha * (Lambda - OldStage) * (myStage - Lambda);
        Real ScaleA = alpha * (Lambda - OldStage) * (Lambda - OldStage);
        Real ScaleB = alpha * (myStage - Lambda) * (myStage - Lambda);

        // retrieve the A-A and B-B energy and size parameters
        // j is used to retrieve the correct set of
        // parameters from the LJ parameter bank        
        j = topo->molecules[M1].type;
 
        // get the LJ parameters for this pair when the atoms are in the '0' state
        params0 = topo->isgLJParms(j,j,topo->atoms[atom1].type,topo->atoms[atom2].type);

        // determine the proper index # into to the bank when the atoms are in the '1' state
        j = topo->molecules[M1].newtype;

        // get the LJ parameters for this pair when the atoms are in the '1' state
        params1 = topo->isgLJParms(j,j,topo->atoms[atom1].type,topo->atoms[atom2].type);
        A0 = (excl != EXCLUSION_MODIFIED ? params0.A : params0.A14);
        B0 = (excl != EXCLUSION_MODIFIED ? params0.B : params0.B14);
        A1 = (excl != EXCLUSION_MODIFIED ? params1.A : params1.A14);
        B1 = (excl != EXCLUSION_MODIFIED ? params1.B : params1.B14);

        // if both of these atoms are dummy atoms (epsilon = 0), then skip
        if (A0 == 0.0 && A1 == 0.0) break;

        // Use A & B to compute epsilon and sigma
        // for the "0" state
        if (A0 != 0.0) {
          sig60 = A0 / B0;
          epsilon0 = B0 / (4.0 * sig60);}
        else {
          sig60 = 1.0;
          epsilon0 = 0.0;}
        // for the "1" state
        if (A1 != 0.0) {
          sig61 = A1 / B1;
          epsilon1 = B1 / (4.0 * sig61);}
        else {
          sig61 = 1.0;
          epsilon1 = 0.0;}

        // scale epsilon by the appropriate factor of lambda
        epsilonA = (myStage - Lambda) * epsilon0;
        epsilonB = (Lambda - OldStage) * epsilon1;

        // soft-core interaction distances
        sr6A = sig60 * r6;
        sr6B = sig61 * r6;
        SoftDistA = ScaleA + 1.0 / sr6A;
        SoftDistB = ScaleB + 1.0 / sr6B;

        // attractive LJ terms    
        Vr6A = 1.0 / SoftDistA;
        Vr6B = 1.0 / SoftDistB;
 
        // repulsive LJ terms
        Vr12A = Vr6A * Vr6A;
        Vr12B = Vr6B * Vr6B;
             
        // get the transformation stage #'s for these atoms
        int atom1_stage = topo->atoms[atom1].stageNumber;
        int atom2_stage = topo->atoms[atom2].stageNumber;
        int max_stage = max(atom1_stage,atom2_stage);

        // if the current stage is not the same as the largest atom stage number then 
        // we do not compute any deltaMu, and if the current stage is less than the largest
        // atom stage number then we only compute interactions in the old identity, if the
        // current stage is larger than both atom stage numbers then we use only the new identity
        if (myStage == max_stage) {
          // energy, force, and chemical potential difference
          energy = 4.0 * (epsilonA * (Vr12A - Vr6A) + epsilonB * (Vr12B - Vr6B));
          force = 24.0 * r2 * (epsilonA * Vr12A / sr6A * (2.0 * Vr6A - 1.0) +
                              epsilonB * Vr12B / sr6B * (2.0 * Vr6B - 1.0));
          deltaMu = 4.0 * (epsilon1 * (Vr12B - Vr6B) - epsilon0 * (Vr12A - Vr6A))
            + 8.0 * LambdaMu * (epsilon1 * Vr12B * (2.0 * Vr6B - 1.0) -
                                epsilon0 * Vr12A * (2.0 * Vr6A - 1.0));
        }
        else if (myStage > max_stage) {
          // energy, force, and chemical potential difference (new identity)
          energy = 4.0 * epsilonB * (Vr12B - Vr6B);
          force = 24.0 * r2 * epsilonB * Vr12B / sr6B * (2.0 * Vr6B - 1.0);
        }
        else {
          // energy, force, and chemical potential difference (old identity)
          energy = 4.0 * epsilonA * (Vr12A - Vr6A);
          force = 24.0 * r2 * epsilonA * Vr12A / sr6A * (2.0 * Vr6A - 1.0);
        }
        break;
        //************************************************************************
      } // end switch structure
      
    }  // end operator()

    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
      (*energies)[ScalarStructure::LENNARDJONES_DELTAMU] += deltaMu;
    }

    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::LENNARDJONES];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 3;}
    void getParameters(std::vector<Parameter>&) const;
    static iSGLennardJonesTableForce make(std::string&, const std::vector<Value>&);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;      
  private:
    Real* myTable;
    Real myRc;
    Real myFac;
    unsigned int myN;

    // constant parameter that affects the soft-core interaction
    // if alpha is set to zero then you get the original Lennard-Jones potential
    Real alpha;
  };

  //______________________________________________________________________ INLINES
}
#endif /* ISGLENNARDJONESTABLEFORCE_H */
