/* -*- c++ -*- */
#ifndef COULOMBTABLEFORCE_H
#define COULOMBTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "LookUpTableBase.h"
#include "CoulombTableForceBase.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ CoulombTableForce

  template<class TSwitchingFunction, unsigned int PRE, typename TReal=Real>
  class CoulombTableForce : public LookUpTableBase<CoulombTableForceBase::LookUpValues,PRE,TReal>, 
			    private CoulombTableForceBase {
  public:
    /// no need for reciprocal squared distance
    enum {DIST_R2=0};
    /// has a valid cutoff
    enum {CUTOFF=1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CoulombTableForce():LookUpTableBase<CoulombTableForceBase::LookUpValues,
					PRE,
					TReal>(),
			myCutoff(0.0),
			myCutoff2(0.0),
			switchingFunction(){}

    /// Constructor with a cutoff switching function
    CoulombTableForce(TSwitchingFunction swf):LookUpTableBase<CoulombTableForceBase::LookUpValues,
					       PRE,
					       TReal>(0.1,
						      swf.cutoffSquared(),
						      2,
						      CoulombTableForceBase::LookUpValues(),
						      swf,
						      128),
					      myCutoff(swf.cutoff()),
					      myCutoff2(swf.cutoffSquared()),
					      switchingFunction(swf){}

    /// Constructor without a cutoff switching function
    CoulombTableForce(TSwitchingFunction swf, 
		      Real rc):LookUpTableBase<CoulombTableForceBase::LookUpValues,
					       PRE,
					       TReal>(0.1,
						      square(rc),
						      2,
						      CoulombTableForceBase::LookUpValues(),
						      swf,
						      128),
			       myCutoff(rc),
			       myCutoff2(rc*rc),
			       switchingFunction(swf){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombTableForce
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
      // get index and interpolation delta dt
      this->index(distSquared,i,dt);

      // compute the coefficients
      Real a = this->myTable[i+0]*q;
      Real b = this->myTable[i+1]*q;
      Real c = this->myTable[i+2]*q;
      Real d = this->myTable[i+3]*q;

      // interpolate energy and force
      this->interpolate(a,b,c,d,dt,energy,force);
    }

    static void accumulateEnergy(ScalarStructure* energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure* energies) {return  (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {
      return keyword + std::string((!TSwitchingFunction::USE) ? std::string("") : std::string(" -switchingFunction " + TSwitchingFunction::getId()));
    }
    static unsigned int getParameterSize() {
      return (TSwitchingFunction::CUTOFF?0:1)+TSwitchingFunction::getParameterSize();
    }
    void getParameters(std::vector<Parameter>& parameters) const{
      switchingFunction.getParameters(parameters);
      if(!TSwitchingFunction::CUTOFF)
	parameters.push_back(Parameter("-cutoff",Value(myCutoff,ConstraintValueType::Positive()),Text("cutoff for table look up")));
    }
    static CoulombTableForce make(std::string& errMsg, const std::vector<Value>& values){
      if(!TSwitchingFunction::CUTOFF)
	return CoulombTableForce(TSwitchingFunction::make(errMsg,std::vector<Value>(values.begin(),values.end()-1)),
				 values[values.size()-1]);
      else
	return CoulombTableForce(TSwitchingFunction::make(errMsg,values));

    }
    Real cutoffSquared() const{return myCutoff2;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff;
    Real myCutoff2;
    TSwitchingFunction switchingFunction;
  };

  //______________________________________________________________________ INLINES
}
#endif /* COULOMBTABLEFORCE_H */
