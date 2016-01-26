/* -*- c++ -*- */
#ifndef EQUILIBRIUMMOLLYINTEGRATOR_H
#define EQUILIBRIUMMOLLYINTEGRATOR_H

#include "MOLLYIntegrator.h"
#include <vector>
namespace ProtoMol {

  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;

  //__________________________________________________ EquilibriumMOLLYIntegrator
  
  class EquilibriumMOLLYIntegrator: public MOLLYIntegrator {
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Static adnd types
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    struct HydrogenBond {
      // Defining a structur to keep all needed information for one constraint
      int a1;     // Global index of atom 1, GenericTopology)
      int a2;     // Global index of atom 1, GenericTopology)
      short i1;   // Local index of atom 1 myHydrogenAtomGroups
      short i2;   // Local index of atom 1 myHydrogenAtomGroups
      Real lambda;// Lambda of the constraint 
      Real l12;   // squared rest lenght of the constrain
    };

    enum { maxGDim=4};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    EquilibriumMOLLYIntegrator();
    EquilibriumMOLLYIntegrator(int cycles, ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    virtual ~EquilibriumMOLLYIntegrator();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class MOLLYIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Vector3DBlock *doAveragingPositions();
    virtual void doMollification(Vector3DBlock *preprocessedPositions);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class MTSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
			    Vector3DBlock *positions,
			    Vector3DBlock *velocities,
			    ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // new methods for class EquilibriumMOLLYIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void luDcmp(Real (&m)[maxGDim][maxGDim], int dim, int (&index)[maxGDim], Real& d) const;
    void luBksb(Real (&m)[maxGDim][maxGDim], int dim, const int (&index)[maxGDim], Real (&b)[maxGDim]) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3DBlock* myAveragedPositions;
    Real myMOLLYStepsize;
        
    std::vector<std::vector<HydrogenBond> > myHydrogenConstraintGroups;
    // Collection of all constraints for each hydrogen group, where
    // a hydrogen group is all bonds containing at least one hydrogen
    // connected together with by bonds.
    std::vector<std::vector<int> > myHydrogenAtomGroups;
  }; 
  
  //______________________________________________________________________ INLINES
}
#endif
