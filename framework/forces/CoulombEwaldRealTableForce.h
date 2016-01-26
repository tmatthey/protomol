/* -*- c++ -*- */
#ifndef COULOMBEWALDREALTABLEFORCE_H
#define COULOMBEWALDREALTABLEFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include "mathutilities.h"
#include "LookUpTableBase.h"
#include "CoulombEwaldRealTableForceBase.h"
#include "CutoffSwitchingFunction.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ CoulombEwaldRealTableForce

  template<class TSwitchingFunction, unsigned int PRE, typename TReal=Real>
  class CoulombEwaldRealTableForce : public LookUpTableBase<CoulombEwaldRealTableForceBase::LookUpValues,PRE,TReal>, 
				     private CoulombEwaldRealTableForceBase {
  public:
    enum {DIST_R2=0};
    enum {CUTOFF=1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CoulombEwaldRealTableForce():LookUpTableBase<CoulombEwaldRealTableForceBase::LookUpValues,
						 PRE,
						 TReal>(),
				 myCutoff(0.0),
				 myCutoff2(0.0),
				 myAlpha(-1.0),
				 switchingFunction(){}

    CoulombEwaldRealTableForce(TSwitchingFunction swf, 
			       Real a):LookUpTableBase<CoulombEwaldRealTableForceBase::LookUpValues,
						       PRE,
						       TReal>(0.1,
							      swf.cutoffSquared(),
							      2,
							      CoulombEwaldRealTableForceBase::LookUpValues(a),
							      swf,
							      128),
				       myCutoff(swf.cutoff()),
				       myCutoff2(swf.cutoffSquared()),
				       myAlpha(a),
				       switchingFunction(swf){}

    CoulombEwaldRealTableForce(TSwitchingFunction swf, 
			       Real rc, 
			       Real a):LookUpTableBase<CoulombEwaldRealTableForceBase::LookUpValues,
						       PRE,
						       TReal>(0.1,
							      square(rc),
							      2,
							      CoulombEwaldRealTableForceBase::LookUpValues(a),
							      swf,
							      128),
				       myCutoff(rc),
				       myCutoff2(rc*rc),
				       myAlpha(a),
				       switchingFunction(swf){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombEwaldRealTableForce
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
      return keyword +  std::string((TSwitchingFunction::getId() != CutoffSwitchingFunction::getId()) ? 
				    std::string(std::string(" -switchingFunction " + TSwitchingFunction::getId())) : std::string(""));
    }
    static unsigned int getParameterSize() {
      return 1+(TSwitchingFunction::CUTOFF?0:1)+TSwitchingFunction::getParameterSize();
    }
    void getParameters(std::vector<Parameter>& parameters) const{
      switchingFunction.getParameters(parameters);
      if(!TSwitchingFunction::CUTOFF)
	parameters.push_back(Parameter("-cutoff",Value(myCutoff,ConstraintValueType::Positive()),Text("cutoff for table look up")));
      parameters.push_back(Parameter("-alpha",Value(myAlpha,ConstraintValueType::Positive()),Text("splitting")));
    }
    static CoulombEwaldRealTableForce make(std::string& errMsg, const std::vector<Value>& values){
      if(!TSwitchingFunction::CUTOFF)
	return CoulombEwaldRealTableForce(TSwitchingFunction::make(errMsg,std::vector<Value>(values.begin(),values.end()-2)),
					  values[values.size()-2],values[values.size()-1]);
      else
	return CoulombEwaldRealTableForce(TSwitchingFunction::make(errMsg,std::vector<Value>(values.begin(),values.end()-1)),
					  values[values.size()-1]);

    }
    Real cutoffSquared() const{return myCutoff2;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff;
    Real myCutoff2;
    Real myAlpha;
    TSwitchingFunction switchingFunction;
  };

  //______________________________________________________________________ INLINES
}
#endif /* COULOMBEWALDREALTABLEFORCE_H */
