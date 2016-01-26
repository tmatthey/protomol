/*  -*- c++ -*-  */
#ifndef HESSIANINT_H
#define HESSIANINT_H

#include "STSIntegrator.h"
#include "Force.h"
#include <fstream>
#include "Vector3DBlock.h"
#include "Hessian.h"
#include "typeSelection.h"

using namespace std;

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class ReducedHessAngle;


  //__________________________________________________ HessianInt
  class HessianInt : public STSIntegrator {

    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    HessianInt();
    HessianInt(Real timestep, std::string evec_s, std::string eval_s, 
	       std::string hess_s, bool sorta, int fm, bool tef, bool masswt, ForceGroup *overloadedForces);
    ~HessianInt(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HessianInt
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
	int diagHessian(double *eigVecO, double *eigValO);
	void outputDiagHess();
	void absSort();
  protected:
    void doKickdoDrift();
    void doHalfKickdoDrift();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 8;}
	virtual void getParameters(std::vector<Parameter>& parameters) const;

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
    typedef TypeSelection::Int<4>::type int32;
	double *eigVec, *eigVal;
	int *eigIndx;
	int totStep;
	unsigned int sz;
	std::string evecfile, evalfile, hessfile;
	bool sortOnAbs;
    unsigned int fixedModes;
	bool textEig, massWeight;
	Hessian hsn;

};
}
#endif



