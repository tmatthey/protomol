/*  -*- c++ -*-  */
#ifndef NORMMODEINT_H
#define NORMMODEINT_H

#include "MTSIntegrator.h"
#include "Vector3DBlock.h"

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //__________________________________________________ NormModeInt
  class NormModeInt : public MTSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormModeInt();
    NormModeInt(int cycles, int fixmodes, Real gamma, int seed, Real temperature, int nve, Real bertau, int fdof, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    ~NormModeInt(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormModeInt
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
	void drift();
  public:
	Vector3DBlock *subspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce);
	Vector3DBlock *subspaceVelocity(Vector3DBlock * force, Vector3DBlock * iPforce);
	Vector3DBlock *nonSubspaceForce(Vector3DBlock * force, Vector3DBlock * iPforce);
	Vector3DBlock *nonSubspacePosition(Vector3DBlock * force, Vector3DBlock * iPforce);
	void subSpaceSift();
	Real expansionPe(Real *hatpe, Real *coupl, int typ);

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
    virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;
  public:
	double *vector3DBlockTOvect(Vector3DBlock* blkDat, double* vecDat);
	Vector3DBlock *vectTOvector3DBlock(double* vecDat, Vector3DBlock* blkDat);

    void restoreState(CheckpointInputStream&);
    void saveState(CheckpointOutputStream&);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    int _N, _m, _rfM, _3N, avItrs, itrs, avMinForceCalc, minForceCalc;
	double *tmpFX, *tmpC;
	Vector3DBlock gaussRandCoordM;

  private:
    Vector3DBlock* ex0;
	int fixMod;
	Real myGamma;
	int mySeed;
	Real myTemp;
	int myNVE;
	Real myBerendsen;
	int fDof;
	double *invSqrtMass, *sqrtMass;
    int numSteps;
    Vector3DBlock gaussRandCoord1, gaussRandCoord2, posUpdt;
	Real orgEnergy;

  };

}

#endif


