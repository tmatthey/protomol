/* -*- c++ -*- */
#ifndef CUTOFFSWITCHINGFUNCTION_H
#define CUTOFFSWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "CutoffSwitchingFunctionBase.h"
namespace ProtoMol {
  //_________________________________________________________________ CutoffSwitchingFunction

  /**
   * Cutoff switching function, implements a simple truncation of the potential.
   */
  class CutoffSwitchingFunction : private CutoffSwitchingFunctionBase {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE=1};
    enum {MODIFY=0};
    enum {CUTOFF=1};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CutoffSwitchingFunction();
    CutoffSwitchingFunction(Real cutoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CutoffSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// simple and fast test if we should apply the switching function
    bool roughTest(Real distSquared) const{return (distSquared <= myCutoff2);}
    Real cutoffSquared() const{return myCutoff2;}
    Real cutoff() const{return myCutoff;}
    void operator()(Real &value, Real &derivOverD, Real distSquared) const{
      derivOverD=0; 
      value = (distSquared > myCutoff2 ? 0.0 : 1.0);
    }

    static const std::string& getId() {return keyword;}
    void getParameters(std::vector<Parameter>& parameters) const;
    static unsigned int getParameterSize() {return 1;}
    static CutoffSwitchingFunction make(std::string& errMsg,std::vector<Value> values);
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff, myCutoff2;
  };

  //______________________________________________________________________ INLINES

}
#endif /* CUTOFFSWITCHINGFUNCTION_H */
