/* -*- c++ -*- */
#ifndef ONEMOLLYPAIR_H
#define ONEMOLLYPAIR_H

#include "Topology.h"
#include "Parameter.h"
#include "ReducedHessAngle.h"
#include "ReducedHessTraits.h"

namespace ProtoMol {
  //_________________________________________________________________ OneMollyPair

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce>
  class OneMollyPair {
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
    OneMollyPair(): switchingFunction(), nonbondedForceFunction(){};
    OneMollyPair(TNonbondedForce nF,
		 TSwitchingFunction sF): 
      switchingFunction(sF), 
      nonbondedForceFunction(nF){};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneMollyPair
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const SemiGenericTopology<TBoundaryConditions>* topo, 
		    const Vector3DBlock* pos, 
		    std::vector< ReducedHessAngle > *angleFilter){
      realTopo      = topo;
      positions     = pos;
      myAngleFilter = angleFilter;
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


    static OneMollyPair make(std::string& errMsg, std::vector<Value> values) {
      unsigned int n = TNonbondedForce::getParameterSize();
      return OneMollyPair(TNonbondedForce::make(errMsg,std::vector<Value>(values.begin(),values.begin()+n)),
			  TSwitchingFunction::make(errMsg,std::vector<Value>(values.begin()+n,values.end())));
    }

    static std::string getId() {
      return "Mollify"+TNonbondedForce::getId() + std::string((!TSwitchingFunction::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunction::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<TBoundaryConditions>* realTopo;
    const Vector3DBlock* positions;
    std::vector< ReducedHessAngle > *myAngleFilter;
    TSwitchingFunction switchingFunction;
    TNonbondedForce nonbondedForceFunction; 
    typename ReducedHessTraits<TNonbondedForce>::Hessian nonbondedHessian;
  };
  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions,
	   class TSwitchingFunction,
	   class TNonbondedForce>
  inline void OneMollyPair<TBoundaryConditions,
			   TSwitchingFunction,
			   TNonbondedForce>::doOneAtomPair(const int i, const int j) {

    Real M1 = realTopo->atoms[i].scaledMass;
    Real M2 = realTopo->atoms[j].scaledMass;
    // get the atom mass

    if(equal(M1,M2)) // mollify method excludes the H-H, O-O, and alike.
      return;
 
    // Get atom distance.
    Vector3D diff = realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j]);
    Real distSquared = diff.normSquared();
    // Do switching function rough test, if necessary.
    if (TSwitchingFunction::USE && !switchingFunction.roughTest(distSquared))
      return;

    // Check for an exclusion.
    ExclusionClass excl = realTopo->exclusions.check(i,j);
    if (excl == EXCLUSION_FULL) 
      return;

    // Calculate the force and energy.
    Real rawEnergy, rawForce;
    Real rDistSquared = 1.0/distSquared;
    nonbondedForceFunction(rawEnergy, rawForce, distSquared, rDistSquared, diff, realTopo, i, j, excl);

    Real switchingValue, switchingDeriv;
    switchingFunction(switchingValue, switchingDeriv, distSquared);

    Matrix3by3 hessian = nonbondedHessian(rawEnergy,
					  rawForce,
					  distSquared, 
					  rDistSquared, 
					  diff, 
					  realTopo, 
					  i, j, 
					  switchingValue, 
					  switchingDeriv,
					  switchingFunction.hessian(diff,distSquared),
					  excl);

    if(M1 < 2){ // now atom1 is Hydrogen
      int angleIndex = static_cast<int>(floor((i+0.1)/3)); // for water ...
      if(realTopo->angles[angleIndex].atom1 == i)
	(*myAngleFilter)[angleIndex].accumulateTo(0,0,hessian);
      else
	(*myAngleFilter)[angleIndex].accumulateTo(2,2,hessian);
    }else{
      int angleIndex = static_cast<int>(floor((j+0.1)/3)); // for water ...
      if(realTopo->angles[angleIndex].atom1 == j)
	(*myAngleFilter)[angleIndex].accumulateTo(0,0,hessian);
      else
	(*myAngleFilter)[angleIndex].accumulateTo(2,2,hessian);
    }

    // End of force computation.

  }
}
#endif /* ONEATOMPAIR_H */
