/* -*- c++ -*- */
#ifndef OSGONEATOMPAIRFULL_H
#define OSGONEATOMPAIRFULL_H

#include "Topology.h"
#include "Parameter.h"
#include "oneAtomContraints.h"

namespace ProtoMol {
  //_________________________________________________________________ oSGOneAtomPairFull

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce,
	   class TConstraint=NoConstraint>
  class oSGOneAtomPairFull {
    // Computes the interaction for a given force between two atoms.

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef TBoundaryConditions BoundaryConditions;
    // Make the boundary conditions visible

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    oSGOneAtomPairFull(): switchingFunction(), nonbondedForceFunction(){};
    oSGOneAtomPairFull(TNonbondedForce nF,
		       TSwitchingFunction sF): 
      switchingFunction(sF), 
      nonbondedForceFunction(nF){};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class oSGOneAtomPairFull
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<TBoundaryConditions>* topo, 
		    const Vector3DBlock* pos, 
		    Vector3DBlock* f, 
		    ScalarStructure* e,
		    const std::vector<Vector3D>* l){
      realTopo  = topo;
      positions = pos;
      forces    = f;
      energies  = e;
      lattice   = l;
    }

    void doOneAtomPair(const int i, const int j);
    // Computes the force and energy for atom i and j.

    void getParameters(std::vector<Parameter>& parameters) const{
      nonbondedForceFunction.getParameters(parameters);
      switchingFunction.getParameters(parameters);
    }
    static unsigned int getParameterSize() {
      return TNonbondedForce::getParameterSize()+TSwitchingFunction::getParameterSize();
    }


    static oSGOneAtomPairFull make(std::string& errMsg, std::vector<Value> values) {
      unsigned int n = TNonbondedForce::getParameterSize();
      return oSGOneAtomPairFull(TNonbondedForce::make(errMsg,std::vector<Value>(values.begin(),values.begin()+n)),
				TSwitchingFunction::make(errMsg,std::vector<Value>(values.begin()+n,values.end())));
    }

    static std::string getId() {
      return TConstraint::getPrefixId() + TNonbondedForce::getId() + TConstraint::getPostfixId() + std::string((!TSwitchingFunction::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunction::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<TBoundaryConditions>* realTopo;
    const Vector3DBlock* positions;
    Vector3DBlock* forces;
    ScalarStructure* energies;
    TSwitchingFunction switchingFunction;
    TNonbondedForce nonbondedForceFunction;
    const std::vector<Vector3D>* lattice;
  };
  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce,
	   class TConstraint>
  inline void oSGOneAtomPairFull<TBoundaryConditions,
				 TSwitchingFunction,
				 TNonbondedForce,
				 TConstraint>::doOneAtomPair(const int i, const int j) {

    if(TConstraint::PRE_CHECK)
      if(!TConstraint::check(realTopo,i,j))
        return;

    // Do we have something to do?
    bool same = (i==j);
    if(same && lattice->empty())
      return;

    Vector3D diffMinimal(realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j]));
    if(!same){
      // Get atom distance.
      Real distSquared = diffMinimal.normSquared();
      // Do switching function rough test, if necessary.
      if (TSwitchingFunction::USE && !switchingFunction.roughTest(distSquared))
	return;

      // Check for an exclusion.
      ExclusionClass excl = realTopo->exclusions.check(i,j);
      if (excl != EXCLUSION_FULL){ 
	// Calculate the force, energy, and chemical potential difference.
	Real rawEnergy, rawForce, rawDeltaMu;
        rawEnergy = rawForce = rawDeltaMu = 0.0;
	Real rDistSquared = 1.0/distSquared;
	nonbondedForceFunction(rawEnergy, rawForce, rawDeltaMu, distSquared, rDistSquared, diffMinimal, realTopo, i, j, excl);
	// Calculate the switched force, energy, and chemical potential.
	Real energy, force, deltaMu;
	if (TSwitchingFunction::USE) {
	  Real switchingValue, switchingDeriv;
	  switchingFunction(switchingValue, switchingDeriv, distSquared);
	  energy = rawEnergy * switchingValue;
	  deltaMu = rawDeltaMu * switchingValue;
	  // This has a - sign because the force is the negative of the 
	  // derivative of the energy (divided by the distance between the atoms).
	  force = rawForce * switchingValue - rawEnergy * switchingDeriv;
	} else {
	  energy = rawEnergy;
	  deltaMu = rawDeltaMu;
	  force = rawForce;
	}
	
	// Add this energy into the total system energy.
        if (nonbondedForceFunction.solute && 
             (realTopo->molecules[ realTopo->atoms[i].molecule ].type == 1 ||
              realTopo->molecules[ realTopo->atoms[j].molecule ].type == 1) ) 
          nonbondedForceFunction.accumulateSoluteEnergy(energies, energy, deltaMu);
        else nonbondedForceFunction.accumulateEnergy(energies, energy, deltaMu);
    
	// Add this force into the atom forces.
	Vector3D fij = -diffMinimal*force;
	(*forces)[i] += fij;
	(*forces)[j] -= fij;

	// compute the vector between molecular centers of mass
	int mi = realTopo->atoms[i].molecule;
	int mj = realTopo->atoms[j].molecule;
	if(mi != mj){
	  Vector3D molDiff = realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
									    realTopo->molecules[mj].position);
	  
	  // Add to the atomic and molecular virials
	  energies->addVirial(fij,-diffMinimal,-molDiff);
	}
	else {
	  energies->addVirial(fij,-diffMinimal);
	}
        if(TConstraint::POST_CHECK)
          TConstraint::check(realTopo,i,j,diffMinimal,energy,fij);
          
      }
    }

    for(unsigned int k=0;k<lattice->size();k++){
      Vector3D diff(diffMinimal+(*lattice)[k]);
      // Get atom distance.
      Real distSquared = diff.normSquared();
      // Do switching function rough test, if necessary.
      if (TSwitchingFunction::USE && !switchingFunction.roughTest(distSquared))
	continue;

      // Calculate the force, energy, and chemical potential difference.
      Real rawEnergy, rawForce, rawDeltaMu;
      rawEnergy = rawForce = rawDeltaMu = 0.0;
      Real rDistSquared = 1.0/distSquared;
      nonbondedForceFunction(rawEnergy, rawForce, rawDeltaMu, distSquared, rDistSquared, diff, realTopo, i, j, EXCLUSION_NONE);
      // Calculate the switched force, energy, and chemical potential difference.
      Real energy, force, deltaMu;
      if (TSwitchingFunction::USE) {
	Real switchingValue, switchingDeriv;
	switchingFunction(switchingValue, switchingDeriv, distSquared);
	energy = rawEnergy * switchingValue;
	deltaMu = rawDeltaMu * switchingValue;
	// This has a - sign because the force is the negative of the 
	// derivative of the energy (divided by the distance between the atoms).
	force = rawForce * switchingValue - rawEnergy * switchingDeriv;
      } else {
	energy = rawEnergy;
	deltaMu = rawDeltaMu;
	force = rawForce;
      }
      // Correct the energy by factor 1/2 when same atom since
      // there is only one pair (i,j) with i==j, where 
      // there are two pairs with same contribution with i !=j
      if(same){
	energy /= 2;
	deltaMu /= 2;
      }
      else {
	// Add this force into the atom forces.
	Vector3D fij = -diff*force;
	(*forces)[i] += fij;
	(*forces)[j] -= fij;	  
	
	// compute the vector between molecular centers of mass
	int mi = realTopo->atoms[i].molecule;
	int mj = realTopo->atoms[j].molecule;
	if(mi != mj){
	  Vector3D molDiff = realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
									    realTopo->molecules[mj].position);
	  
	  // Add to the atomic and molecular virials
	  energies->addVirial(fij,-diff,-molDiff);
	}
	else {
	  energies->addVirial(fij,-diff);
	}
      }
      // Add this energy into the total system energy.
      if (nonbondedForceFunction.solute && 
           (realTopo->molecules[ realTopo->atoms[i].molecule ].type == 1 ||
            realTopo->molecules[ realTopo->atoms[j].molecule ].type == 1) ) 
        nonbondedForceFunction.accumulateSoluteEnergy(energies, energy, deltaMu);
      else nonbondedForceFunction.accumulateEnergy(energies, energy, deltaMu);
      
      if(TConstraint::POST_CHECK)
        TConstraint::check(realTopo,i,j,diff,energy,-diff*force);
    }
  }
}
#endif /* ONEATOMPAIR_H */
