/*  -*- c++ -*-  */
#ifndef BSPLINEMOLLYINTEGRATOR_H
#define BSPLINEMOLLYINTEGRATOR_H


#include "MOLLYIntegrator.h"
#include "BSplineType.h"
#include <vector>
#include <map>

namespace ProtoMol {

  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;
  class StandardIntegrator;
  class ReducedHessAngle;
  
  //__________________________________________________ BSplineMOLLYIntegrator
  /**
     This integrator is developped using an approximation to the dynamics:
     let the heavy atoms always be fixed in their position such that we could
     use the reduced forms of Hessians. This is a physically and mathematically
     correct approximation.
  */

  class BSplineMOLLYIntegrator: public MOLLYIntegrator {


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
    BSplineMOLLYIntegrator();
    BSplineMOLLYIntegrator(int cycles,
			   const BSplineType& typeOfBSpline,
			   Real mollyStepsize,
			   ForceGroup *overloadedForces,
			   StandardIntegrator *nextIntegrator);
    
    virtual ~BSplineMOLLYIntegrator();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 3;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MOLLYIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Vector3DBlock *doAveragingPositions();
    virtual void doMollification(Vector3DBlock *preprocessedPositions); 
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
			    Vector3DBlock *positions,
			    Vector3DBlock *velocities,
			    ScalarStructure *energies);

  private:
    virtual void doUncache();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class MTSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods for class BSplineMOLLYIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void calcMOLLYForcesHalfKickOneDrift();
    void calcMOLLYForcesOneKickOneDrift();
    void calculateMOLLYForcesBonded();
    void calculateMOLLYForcesHBonded();
    void setMyNumIterMOLLYStepsize();
    void updateXx() ;
    void updateB_Bx_Px_for1Kick();
    void updateB_Bx();
    void updateB_Bx_Px_for1stHalfKick();
    void calcHessiansBondsAnglesHBonds();
    void init();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    bool myCached;

    ForceGroup *myMOLLYForcesBonded;
    // could be bond or bond+angle
    
    ForceGroup *myMOLLYForcesHBonded;
    // could be LJ or LJ+Coulomb. We want to separate the ForceGroup into bonded
    // and Hydrogen bonded
    
    ForceGroup *myHBondForces;
    //this is the range force (2.5 ~ 4 Angstrom) to be mollified.
    Vector3DBlock   *myMOLLYPositions;
    Vector3DBlock   *myMOLLYVelocities;
    Vector3DBlock   *myMOLLYForces;
    ScalarStructure *myMOLLYEnergies;
    
    Vector3DBlock *myB;
    std::vector<ReducedHessAngle> *myAngleFilter;  
    std::vector<ReducedHessAngle> *myPxAngle;  
    std::vector<ReducedHessAngle> *myXxAngle;
    std::vector<ReducedHessAngle> *myBxAngle;
    
    BSplineType myTypeOfBSpline;
    Real myMOLLYStepsize;
    Real myMOLLYStepsizeP;
    Real myMOLLYStepsizeSquare;
    unsigned int myNumIter;

    bool myBond;
    bool myAngle;
    bool myCoulomb;
    bool myLennardJones;
    std::vector<AngleIndex> myAngleIndexes;
  }; 
}
#endif
