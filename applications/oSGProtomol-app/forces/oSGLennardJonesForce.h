/* -*- c++ -*- */
#ifndef OSGLENNARDJONESFORCE_H
#define OSGLENNARDJONESFORCE_H

#include "GenericTopology.h"
#include "LennardJonesParameters.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include <string>

// uncomment for debugging purposes
//#define DEBUG_NORMAL_LJ
//#define DEBUG_SOFTCORE
//#define DEBUG_LJ_DISTANCE_CHECK

namespace ProtoMol {

  //_________________________________________________________________ oSGLennardJonesForce
  class oSGLennardJonesForce {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    oSGLennardJonesForce() {}
    oSGLennardJonesForce(bool mySolute) {
      solute = mySolute;
      if (solute)
        Report::report << Report::plain << "Solute interaction energy will be reported as 'Other'." << Report::endr;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class oSGLennardJonesForce
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
	Report::report << Report::warning << "oSGLennardJonesForce::operator(): atom "<<atom1<<" and "<<atom2<<" get closer"
	       << " than 0.5 angstroms (="<<sqrt(1.0/rDistSquared)<<")." << Report::endr;   
#endif
      
      // Fast LJ computations
      Real r2   = rDistSquared;
      Real r6   = power(r2,3);
      Real r12, rij_6;

      // attractive and repulsive LJ terms
      Real Vr6 = 0.0;
      Real Vr12 = 0.0;

      // energy and size parameter variables
      LennardJonesParameters params;
      Real A = 0.0;
      Real B = 0.0;

      // transforming atom indicators
      bool atom1_scaled = false;
      bool atom2_scaled = false;
      bool deletion = false;
      int myStage = 0;
      int atom_stage = 0;

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
      // we do not scale any of these (intramolecular) terms for 
      // molecule insertions or deletions, so we must treat this case
      // as if neither atom is being transformed 
      else if ( (atom1_scaled) && (atom2_scaled) ) {Choice = 0;}
      
      // one of the atoms is being transformed
      else {
        // determine the current transformation stage of the molecule
        // and get the transformation stage # for this atoms       
        if (atom1_scaled) {
          myStage = static_cast<int>(floor(topo->molecules[M1].lambda)) + 1;
          atom_stage = topo->atoms[atom1].stageNumber;
          if (topo->molecules[M1].type != topo->molecules[M1].newtype) deletion = true;
        }
        else {
          myStage = static_cast<int>(floor(topo->molecules[M2].lambda)) + 1;
          atom_stage = topo->atoms[atom2].stageNumber;
          if (topo->molecules[M2].type != topo->molecules[M2].newtype) deletion = true;
        }
        
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
        
        if  (!(atom1_scaled && atom2_scaled)) {
          /// check to see if atom1 may have already been successfully transformed in a previous stage  
          if (atom1_scaled) {
            // this atom has already been transformed if myStage > atom_stage
            if (myStage > atom_stage) {
              // check to see if we are doing a deletion.  If so, then this atom has already been deleted
              // and we should skip the calculation 
              if (topo->molecules[M1].type != topo->molecules[M1].newtype) deletion = true;
            }
          } // end if statement
          /// check to see if atom2 may have already been successfully transformed in a previous stage
          if (atom2_scaled) {
            // this atom has already been transformed if myStage > atom_stage
            if (myStage > atom_stage) {
              // check to see if we are doing a deletion.  If so, then this atom has already been deleted
              // and we should skip the calculation 
              if (topo->molecules[M2].type != topo->molecules[M2].newtype) deletion = true;
            }
          } // end if statement
        } // end outer if statement
        
	// skip if one of the atoms was deleted in a previous step
	if (deletion) break;
	
	// get the LJ parameters for this pair
        params = topo->lennardJonesParameters(topo->atoms[atom1].type,topo->atoms[atom2].type);
	A = (excl != EXCLUSION_MODIFIED ? params.A : params.A14);
	B = (excl != EXCLUSION_MODIFIED ? params.B : params.B14);

	// if any of these atoms have zero energy parameter then skip
	if (A == 0.0) break;

	// attractive LJ term    
	Vr6 = B * r6;
        
	// repulsive LJ term
	r12 = r6 * r6;
	Vr12 = A * r12;
	
	// energy and force
	energy = Vr12 - Vr6;
	force = 12.0 * Vr12 * r2 - 6.0 * r2 * Vr6;

#ifdef DEBUG_NORMAL_LJ
if (1.0/sqrt(r2) < 3.0) {
    	Report::report.precision(6);
	Report::report << Report::plain << "i = " << atom1+1 << ", j = " << atom2+1 
                       << ", rij = " << 1.0/sqrt(r2) << ", energy = " << energy << Report::endr;
}
#endif
	break;
        //----------------------------------------------------------
      case 1:  
	// interaction between non-transforming atom and a transforming atom
        // soft-core interaction variables
        Real Lambda;
        Real alpha;
        Real E_preFactor;
        Real DMU_outerPreFactor;
        Real DMU_innerPreFactor;
        Real distFactor;

	if (atom1_scaled) {
	  // get the lambda factor
	  Lambda = topo->molecules[M1].lambda;

          // get the alphaLJ parameter for this atom
          alpha = topo->atoms[atom1].alphaLJ;
	}
	else {
	  // get the lambda factor
	  Lambda = topo->molecules[M2].lambda;

          // get the alphaLJ parameter for this atom
          alpha = topo->atoms[atom2].alphaLJ;
	}


        // compute the LJ epsilon and sigmas
        // get the LJ A & B parameters for this pair
        params = topo->lennardJonesParameters(topo->atoms[atom1].type,topo->atoms[atom2].type);
	A = params.A;
	B = params.B;

        // if the either of these atoms has zero energy parameter then skip
        if (A == 0.0) break;

        // determine the current transformation stage of the molecule
        int OldStage = myStage - 1;
           
        // LJ soft-core energy scaling factors
        Real LambdaMu = 2.0 * alpha * (Lambda - OldStage) * (myStage - Lambda);
        if (deletion) { /// for molecule deletion
          E_preFactor = power((myStage - Lambda),4);
          DMU_outerPreFactor = -1.0 * power((myStage - Lambda),3);
          DMU_innerPreFactor = 4.0;
          distFactor = alpha * power((Lambda - OldStage),2);
        }
        else { /// for molecule insertion
          E_preFactor = power((Lambda - OldStage),4);
          DMU_outerPreFactor = power((Lambda - OldStage),3);
          DMU_innerPreFactor = 4.0;
          distFactor = alpha * power((myStage - Lambda),2);
        }

        // soft-core interaction distances
        rij_6 = DistSquared * DistSquared * DistSquared;
        Real SoftDist = rij_6 + (A / B) * distFactor;

        // attractive LJ term 
        Vr6 = B / SoftDist;
 
        // repulsive LJ term
        Vr12 = A / (SoftDist * SoftDist);

        // energy, force, and chemical potential difference
        Real rawEnergy = Vr12 - Vr6;
        energy = E_preFactor * rawEnergy;
        force = 6.0 * r2 / r6 * (E_preFactor * Vr12 * (2.0 * Vr6 / B - B / A));
        DeltaMu = DMU_innerPreFactor * rawEnergy + LambdaMu * Vr12 * (2.0 * Vr6 * A / (B * B) - 1.0);
        DeltaMu *= DMU_outerPreFactor;

#ifdef DEBUG_SOFTCORE
if (1.0/sqrt(r2) < 4.0) {
        Report::report.precision(6);
        Report::report << Report::plain << "i = " << atom1+1 << ", j = " << atom2+1 
                       << ", rij = " << 1.0/sqrt(r2) << ", Energy = " << energy 
                       << ", force = " << force << ", DeltaMu = " << DeltaMu << Report::endr;
}
#endif //DEBUG_SOFTCORE
	break;
	//************************************************************************
      } // end switch structure
    } // end operator()
 

    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
      (*energies)[ScalarStructure::LENNARDJONES_DELTAMU] += deltaMu;
    }
    
    /// separate accumulateEnergy function for computing coarse-grained potentials
    /// here all solute-solvent interactions are reported in ScalarStructure::OTHER, while
    /// all solvent-solvent interactions are reported in the usual ScalarStructure::LENNARDJONES
    static void accumulateSoluteEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::OTHER] += energy;
      (*energies)[ScalarStructure::LENNARDJONES_DELTAMU] += deltaMu;
    }
    
    static Real getEnergy(const ScalarStructure* energies) {
      return  (*energies)[ScalarStructure::LENNARDJONES];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 1;}
    void getParameters(std::vector<Parameter>& parameters) const{
      parameters.push_back(Parameter("-soluteEnergy", Value(solute), false));
      //if (!(solute == "yes") || !(solute == "no")) 
        //Report::report << Report::error << "parameter soluteEnergy: " 
        //               << solute << " for " << getId() << " not valid." << Report::endr;
    }
    static oSGLennardJonesForce make(std::string& , const std::vector<Value>& values) {
      return oSGLennardJonesForce(values[0]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;    
    bool solute;
    
  private:

  };
  //______________________________________________________________________ INLINES
}

#endif /* OSGLENNARDJONESFORCE_H */
