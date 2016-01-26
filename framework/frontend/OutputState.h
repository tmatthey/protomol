/*  -*- c++ -*-  */
#ifndef OUTPUTSTATE_H
#define OUTPUTSTATE_H

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "OutputFile.h"
#include "Vector3DBlock.h"


namespace ProtoMol {

  class XYZTrajectoryWriter;

  //________________________________________________________ OutputState
  /** 
      Writes the dihedral values to a file at given freqeuncy.
  */
  class OutputState : public OutputFile {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputState();
    OutputState(const std::string& filename, int freq, int cacheFreq, int cacheSize, 
		    Real closeTime, std::string statesInFile);
    ~OutputState();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputState
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class OutputFile
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual void doRunCached(int step);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const;
    virtual void doInitialize();
    virtual void doFinalize(int step);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{ return keyword;}
    // Returns the identification string

    virtual unsigned int getParameterSize() const {return 6;}
    virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;
    virtual void getParameters(std::vector<Parameter> &parameter) const;
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    static const std::string keyword;
    
    struct orderParam{
      int dihID;
      Real lowBound;
      Real highBound;
    };

  private:
    bool AObserved;
    bool BObserved;
    std::string myStatesFile;
    std::vector< orderParam > stateA;
    std::vector< orderParam > stateB;
    XYZTrajectoryWriter* myVELA;
    XYZTrajectoryWriter* myVELB;
    XYZTrajectoryWriter* myPOSA;
    XYZTrajectoryWriter* myPOSB;
  };
}
#endif
