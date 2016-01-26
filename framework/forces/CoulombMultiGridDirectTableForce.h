/* -*- c++ -*- */
#ifndef COULOMBMULTIGRIDDIRECTTABLEFORCE_H
#define COULOMBMULTIGRIDDIRECTTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "LookUpTableBase.h"
#include "CoulombMultiGridDirectTableForceBase.h"
#include "CutoffSwitchingFunction.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ CoulombMultiGridDirectTableForce

  template<class TKernel, unsigned int PRE, typename TReal=Real>
  class CoulombMultiGridDirectTableForce : public LookUpTableBase<CoulombMultiGridDirectTableForceBase::LookUpValues<TKernel>,PRE,TReal>, 
					   private CoulombMultiGridDirectTableForceBase {
  public:
    enum {DIST_R2=0};
    enum {CUTOFF=1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CoulombMultiGridDirectTableForce():LookUpTableBase<CoulombMultiGridDirectTableForceBase::LookUpValues<TKernel>,
						       PRE,
						       TReal>(),
				       myCutoff(0.0),
				       myCutoff2(0.0),
				       myS(-1.0){}

    CoulombMultiGridDirectTableForce(Real rc, 
				     Real s):LookUpTableBase<CoulombMultiGridDirectTableForceBase::LookUpValues<TKernel>,
							     PRE,
							     TReal>(0.1,
								    square(rc),
								    2,
								    CoulombMultiGridDirectTableForceBase::LookUpValues<TKernel>(s),
								    CutoffSwitchingFunction(rc),
								    128),
					     myCutoff(rc),
					     myCutoff2(rc*rc),
					     myS(s){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombMultiGridDirectTableForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force, 
		    Real distSquared,
		    Real /*rDistSquared*/,
		    const Vector3D& /*diff*/, 
		    const GenericTopology* topo, 
		    int atom1, int atom2, 
		    ExclusionClass excl) const{

      Real q = topo->atoms[atom1].scaledCharge * topo->atoms[atom2].scaledCharge *
	((topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED) ?topo->coulombScalingFactor:1.0);

      Real dt;
      int i;
      this->index(distSquared,i,dt);

      Real a = this->myTable[i+0]*q;
      Real b = this->myTable[i+1]*q;
      Real c = this->myTable[i+2]*q;
      Real d = this->myTable[i+3]*q;

      this->interpolate(a,b,c,d,dt,energy,force);
    }

    static void accumulateEnergy(ScalarStructure* energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {
      return keyword + " -kernel "+TKernel::keyword;
    }

    static unsigned int getParameterSize() {
      return 2;
    }
    void getParameters(std::vector<Parameter>& parameters) const{
      parameters.push_back(Parameter("-cutoff",Value(myCutoff,ConstraintValueType::Positive()),Text("cutoff for table look up")));
      parameters.push_back(Parameter("-s",Value(myS,ConstraintValueType::Positive()),Text("MG smoothing distance")));
    }
    static CoulombMultiGridDirectTableForce make(std::string&, const std::vector<Value>& values){
      return CoulombMultiGridDirectTableForce(values[0],values[1]);

    }
    Real cutoffSquared() const{return myCutoff2;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff;
    Real myCutoff2;
    Real myS;
  };

  //______________________________________________________________________ INLINES
}
#endif /* COULOMBMULTIGRIDDIRECTTABLEFORCE_H */
