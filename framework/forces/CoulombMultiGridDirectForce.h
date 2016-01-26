/* -*- c++ -*- */
#ifndef COULOMBMULTIGRIDDIRECTFORCE_H
#define COULOMBMULTIGRIDDIRECTFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "CoulombMultiGridDirectForceBase.h"
#include "Parameter.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ CoulombMultiGridDirectForce
  template<class TKernel>  
  class CoulombMultiGridDirectForce : private CoulombMultiGridDirectForceBase {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Potential to compute the direct part of MultiGrid
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  public:
    enum {DIST_R2=0};
    enum {CUTOFF=0};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CoulombMultiGridDirectForce():myS(0.0){}
    CoulombMultiGridDirectForce(Real s):myS(s),myRS(1.0/s){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombMultiGridDirectForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force, 
		    Real distSquared,
		    Real /*rDistSquared*/,
		    const Vector3D &, 
		    const GenericTopology* topo, 
		    int atom1, int atom2,
		    ExclusionClass excl) const {

      Real qq = topo->atoms[atom1].scaledCharge * topo->atoms[atom2].scaledCharge;
      if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED) 
	qq *= topo->coulombScalingFactor;

      Real r = sqrt(distSquared);
      Real rr = 1.0/r;      
      force  = (-TKernel::dKernelR(rr) + TKernel::dSmooth(r,myS,myRS))*rr*qq;
      energy = ( TKernel::kernelR(rr)  - TKernel::smooth(r,myS,myRS))*qq; 
    }

    static void accumulateEnergy(ScalarStructure* energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {return keyword+" -kernel "+TKernel::keyword;}
    static unsigned int getParameterSize() {return 1;}
    void getParameters(std::vector<Parameter>&parameters) const{parameters.push_back(Parameter("-s",Value(myS,ConstraintValueType::Positive()),Text("MG smoothing distance")));}
    static CoulombMultiGridDirectForce make(std::string&, const std::vector<Value>& values){return CoulombMultiGridDirectForce(static_cast<Real>(values[0]));}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myS;
    Real myRS;
  };

  //______________________________________________________________________ INLINES
}
#endif /* COULOMBMULTIGRIDDIRECTFORCE_H */
