/* -*- c++ -*- */
#ifndef ISGONEATOMPAIR_H
#define ISGONEATOMPAIR_H

#include "Topology.h"
#include "Parameter.h"
#include "oneAtomContraints.h"

namespace ProtoMol {
  //_________________________________________________________________ ISGOneAtomPair

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce,
	   class TConstraint=NoConstraint>
  class iSGOneAtomPair {
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
    iSGOneAtomPair(): switchingFunction(), nonbondedForceFunction(){};
    iSGOneAtomPair(TNonbondedForce nF,
		   TSwitchingFunction sF): 
      switchingFunction(sF), 
      nonbondedForceFunction(nF),mySquaredCutoff(sF.cutoffSquared()){};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ISGOneAtomPair
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<TBoundaryConditions>* topo, 
		    const Vector3DBlock* pos, Vector3DBlock* f, ScalarStructure* e){
      realTopo  = topo;
      positions = pos;
      forces    = f;
      energies  = e;
    }

    void doOneAtomPair(const int i, const int j);
    // Computes the force, energy, and chemical potential difference for atom i and j.

    void getParameters(std::vector<Parameter>& parameters) const{
      nonbondedForceFunction.getParameters(parameters);
      switchingFunction.getParameters(parameters);
    }

    static unsigned int getParameterSize() {
      return TNonbondedForce::getParameterSize()+TSwitchingFunction::getParameterSize();
    }


    static iSGOneAtomPair make(std::string& errMsg, std::vector<Value> values) {
      unsigned int n = TNonbondedForce::getParameterSize();
      return iSGOneAtomPair(TNonbondedForce::make(errMsg,std::vector<Value>(values.begin(),values.begin()+n)),
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
    Real mySquaredCutoff;
  };
  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce,
	   class TConstraint>
  inline void iSGOneAtomPair<TBoundaryConditions,
			     TSwitchingFunction,
			     TNonbondedForce,
			     TConstraint>::doOneAtomPair(const int i, const int j) {

    if(TConstraint::PRE_CHECK)
      if(!TConstraint::check(realTopo,i,j))
        return;

    // Get atom distance.
    Vector3D diff = realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j]);
    Real distSquared = diff.normSquared();
    // Do switching function rough test, if necessary.
    if (TSwitchingFunction::USE && distSquared > mySquaredCutoff)
      return;

    // Check for an exclusion.
    int mi = realTopo->atoms[i].molecule;
    int mj = realTopo->atoms[j].molecule;
    bool same = (mi==mj);
    ExclusionClass excl = (same?realTopo->exclusions.check(i,j):EXCLUSION_NONE);
    if (excl == EXCLUSION_FULL) 
      return;

    // Calculate the force, energy, and chemical potential difference.
    Real energy, force, deltaMu;
    energy = force = deltaMu = 0.0;
    Real rDistSquared = 1.0/distSquared;
    nonbondedForceFunction(energy, force, deltaMu, distSquared, rDistSquared, diff, realTopo, i, j, excl);

    // Calculate the switched force, energy, and chemical potential difference.
    if (TSwitchingFunction::USE) {
      Real switchingValue, switchingDeriv;
      switchingFunction(switchingValue, switchingDeriv, distSquared);
      // This has a - sign because the force is the negative of the 
      // derivative of the energy (divided by the distance between the atoms).
      force = force * switchingValue - energy * switchingDeriv;
      energy = energy * switchingValue;
      deltaMu = deltaMu * switchingValue;
    }
    
    // Add this energy into the total system energy.
    nonbondedForceFunction.accumulateEnergy(energies, energy, deltaMu);
    // Add this force into the atom forces.
    Vector3D fij(diff*force);
    (*forces)[i] -= fij;
    (*forces)[j] += fij;

    if(!same && energies->molecularVirial()){
      
      // Add to the atomic and molecular virials
      energies->addVirial(fij,diff,
                          realTopo->boundaryConditions.minimalDifference(realTopo->molecules[mi].position,
                                                                         realTopo->molecules[mj].position));
    }
    else if(energies->virial()) {
      energies->addVirial(fij,diff);
    }
    // End of force computation.
    if(TConstraint::POST_CHECK)
      TConstraint::check(realTopo,i,j,diff,energy,fij);

  }
}
#endif /* ISGONEATOMPAIR_H */
