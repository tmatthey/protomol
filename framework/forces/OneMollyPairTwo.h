/* -*- c++ -*- */
#ifndef ONEMOLLYPAIRTWO_H
#define ONEMOLLYPAIRTWO_H

#include "Topology.h"
#include "Parameter.h"
#include "ReducedHessAngle.h"
#include "ReducedHessTraits.h"

namespace ProtoMol {

  //_____________________________________________________________ OneMollyPairTwo

  template<class TBoundaryConditions,
	   class TSwitchingFunctionFirst,
	   class TNonbondedForceFirst, 
	   class TSwitchingFunctionSecond,
	   class TNonbondedForceSecond>
  class OneMollyPairTwo {
    // Computes the interaction for a given force between to atoms.
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
    OneMollyPairTwo(){}
    OneMollyPairTwo(TNonbondedForceFirst f1, 
		    TNonbondedForceSecond f2,		   
		    TSwitchingFunctionFirst sF1, 
		    TSwitchingFunctionSecond sF2):
      switchingFunctionFirst(sF1),
      nonbondedForceFunctionFirst(f1),
      switchingFunctionSecond(sF2),
      nonbondedForceFunctionSecond(f2){};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OneMollyPairTwo
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
      nonbondedForceFunctionFirst.getParameters(parameters);
      nonbondedForceFunctionSecond.getParameters(parameters);
      switchingFunctionFirst.getParameters(parameters);
      switchingFunctionSecond.getParameters(parameters);
    }

    static unsigned int getParameterSize(){
      return 
	TNonbondedForceFirst::getParameterSize()+
	TNonbondedForceSecond::getParameterSize()+
	TSwitchingFunctionFirst::getParameterSize()+					    
	TSwitchingFunctionSecond::getParameterSize();
    }

    static OneMollyPairTwo make(std::string& errMsg, std::vector<Value> values) {
      unsigned int l1 = TNonbondedForceFirst::getParameterSize();
      unsigned int l2 = TNonbondedForceSecond::getParameterSize()+l1;
      unsigned int l3 = TSwitchingFunctionFirst::getParameterSize()+l2;
      return OneMollyPairTwo(TNonbondedForceFirst::make(errMsg,std::vector<Value>(values.begin(),values.begin()+l1)),			    
			     TNonbondedForceSecond::make(errMsg,std::vector<Value>(values.begin()+l1,values.begin()+l2)),
			     TSwitchingFunctionFirst::make(errMsg,std::vector<Value>(values.begin()+l2,values.begin()+l3)),
			     TSwitchingFunctionSecond::make(errMsg,std::vector<Value>(values.begin()+l3,values.end())));

    }

    static std::string getId() {
      return "Mollify"+TNonbondedForceFirst::getId() + " Mollify"+TNonbondedForceSecond::getId() + std::string((!(TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE)) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunctionFirst::getId() + " -switchingFunction " +TSwitchingFunctionSecond::getId()));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    const SemiGenericTopology<TBoundaryConditions>* realTopo;
    const Vector3DBlock* positions;
    std::vector< ReducedHessAngle > *myAngleFilter;
    TSwitchingFunctionFirst switchingFunctionFirst;
    TNonbondedForceFirst nonbondedForceFunctionFirst;  
    typename ReducedHessTraits<TNonbondedForceFirst>::Hessian nonbondedHessianFirst;
    TSwitchingFunctionSecond switchingFunctionSecond;
    TNonbondedForceSecond nonbondedForceFunctionSecond;  
    typename ReducedHessTraits<TNonbondedForceSecond>::Hessian nonbondedHessianSecond;
  };
  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions,
	   class TSwitchingFunctionFirst,
	   class TNonbondedForceFirst, 
	   class TSwitchingFunctionSecond,
	   class TNonbondedForceSecond>
  inline void OneMollyPairTwo<TBoundaryConditions,
			      TSwitchingFunctionFirst,
			      TNonbondedForceFirst, 
			      TSwitchingFunctionSecond,
			      TNonbondedForceSecond>::doOneAtomPair(const int i, const int j) {

    Real M1 = realTopo->atoms[i].scaledMass;
    Real M2 = realTopo->atoms[j].scaledMass;
    // get the atom mass

    if(equal(M1,M2)) // mollify method excludes the H-H, O-O, and alike.
      return;
 
    // Get atom distance.
    Vector3D diff = realTopo->boundaryConditions.minimalDifference((*positions)[i], (*positions)[j]);
    Real distSquared = diff.normSquared();
    // Do switching function rough test, if necessary.
    if ((TSwitchingFunctionFirst::USE || TSwitchingFunctionSecond::USE) && !switchingFunctionFirst.roughTest(distSquared) 
	&& !switchingFunctionSecond.roughTest(distSquared))
      return;

    // Check for an exclusion.
    ExclusionClass excl = realTopo->exclusions.check(i,j);
    if (excl == EXCLUSION_FULL) 
      return;

    // Calculate the force and energy.
    Real rawEnergy1, rawForce1, rawEnergy2, rawForce2;
    Real rDistSquared = 1.0/distSquared;
    nonbondedForceFunctionFirst(rawEnergy1, rawForce1, distSquared, rDistSquared, diff, realTopo, i, j, excl);
    nonbondedForceFunctionSecond(rawEnergy2, rawForce2, distSquared, rDistSquared, diff, realTopo, i, j, excl);

    // Calculate the switched force and energy.
    Real switchingValue1, switchingDeriv1;
    Real switchingValue2, switchingDeriv2;
    switchingFunctionFirst(switchingValue1, switchingDeriv1, distSquared);
    switchingFunctionSecond(switchingValue2, switchingDeriv2, distSquared);

    Matrix3by3 hessian = 
      nonbondedHessianFirst(rawEnergy1,
			    rawForce1,
			    distSquared, 
			    rDistSquared, 
			    diff, 
			    realTopo, 
			    i, j, 
			    switchingValue1, 
			    switchingDeriv1,
			    switchingFunctionFirst.hessian(diff,distSquared),
			    excl) 
      +
      nonbondedHessianSecond(rawEnergy2,
			     rawForce2,
			     distSquared, 
			     rDistSquared, 
			     diff, 
			     realTopo, 
			     i, j, 
			     switchingValue2, 
			     switchingDeriv2,
			     switchingFunctionSecond.hessian(diff,distSquared),
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

#endif /* ONEATOMPAIRTWO_H */
