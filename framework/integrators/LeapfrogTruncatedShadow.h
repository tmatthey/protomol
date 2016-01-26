/*  -*- c++ -*-  */
#ifndef LEAPFROGTRUNCATEDSHADOW_H
#define LEAPFROGTRUNCATEDSHADOW_H

#include "STSIntegrator.h"
#include "Force.h"
#include <fstream>
#include "Vector3DBlock.h"
using namespace std;

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class ReducedHessAngle;


  //__________________________________________________ LeapfrogTruncatedShadow
  class LeapfrogTruncatedShadow : public STSIntegrator {

    struct BondIndex {
      BondIndex(int a, int b):b1(std::min(a,b)),b2(std::max(a,b)){}
      int b1,b2;
      bool operator<(const BondIndex& a) const {
	if(b1 < a.b1)
	  return true;
	else if(b1 > a.b1)
	  return false;
	else
	  return (b2 < a.b2);
      }
    };

    struct AngleIndex {
      AngleIndex(int a, int b, int c):angle(a),bond1(b),bond2(c){}
      int angle,bond1,bond2;
    };

    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LeapfrogTruncatedShadow();
    LeapfrogTruncatedShadow(Real timestep, ForceGroup *overloadedForces);
    ~LeapfrogTruncatedShadow(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LeapfrogTruncatedShadow
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doKickdoDrift();
    void doHalfKickdoDrift();
    Real calcNearHamiltonian();
    Real outputHessian();
    Real calcPairInteractionHess(Real doSwitch,Real switchon,Real cutoff,Real order,Real switchoff,int lj);
	void TorsionHess(const Torsion& currTorsion, double * hessD);//void Dihedral(int dh, double * hessD);
	double *rotateV3D(double *aRot, double *mf);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 1;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
			    Vector3DBlock   *positions,
			    Vector3DBlock   *velocities, 
			    ScalarStructure *energies);
    virtual void run(int numTimesteps);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg)const;
  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    bool myBond;
    bool myAngle;
    bool myCoulomb;
    bool myCoulombDielec;
    bool myCoulombSCPISM;
	bool myCoulombBornRadii;
    bool myLennardJones;
	bool myDihedral;
	bool myImproper;
    Real cCutoff, cSwitchon, cSwitch, cOrder, cSwitchoff;
    Real lCutoff, lSwitchon, lSwitch, lOrder, lSwitchoff;
    Real D, S, epsi;
	int swt;

};
}
#endif


