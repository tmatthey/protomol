/*  -*- c++ -*-  */
#ifndef DIHEDRALHMCINTEGRATOR_H
#define DIHEDRALHMCINTEGRATOR_H

#include "MTSIntegrator.h"

namespace ProtoMol {
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;

  //____________________________________________________________ DihedralHMCIntegrator

  class DihedralHMCIntegrator: public MTSIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // totally overload the function.  Use no defaults (unless the lists passed
    //   are empty lists of forces.  Then the default forces will be assigned);
    //   however, the lists still must be passed nonetheless.
    DihedralHMCIntegrator();  
    DihedralHMCIntegrator(int cycles,
			Real initialTemperature,
			bool randomCycLen,
			bool dihset,
			std::string dsetfile,
                          bool angset,
                          std::string asetfile,
			ForceGroup *overloadedForces,
			StandardIntegrator *nextIntegrator);



    virtual ~DihedralHMCIntegrator();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 6;}


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo,
			    Vector3DBlock *positions,
			    Vector3DBlock *velocities,
			    ScalarStructure *energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void run (int numTimesteps);

  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class MTSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods for class DihedralHMCIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    bool metropolisTest(Real newEnergy, 
			Real oldEnergy,
			Real theTemperature,
			Real& acceptProb);
    void saveValues();
    void restoreValues();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    Vector3DBlock*   myOldPositions;
    Vector3DBlock*   myOldVelocities;
    ScalarStructure* myOldEnergies;
    Real myInitialTemperature;
    bool myDihedralsSet;
    std::string myDhmcDiSetFile;
    bool myAnglesSet;    
    std::string myDhmcAnSetFile;
    int myDihedralIndex;
    std::vector< int >* myDihedrals;
    std::vector< Real >* myAngles;
  };
}
#endif
