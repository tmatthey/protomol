/* -*- c++ -*- */
#ifndef SHIFTSWITCHINGFUNCTION_H
#define SHIFTSWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "ShiftSwitchingFunctionBase.h"
namespace ProtoMol {
  //_________________________________________________________________ ShiftSwitchingFunction

  /**
   * Shift switching function, moving the potential up or down such that 
   * potential is 0 for the value cutoff.
   */
  class ShiftSwitchingFunction : private ShiftSwitchingFunctionBase {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE=1};
    enum {MODIFY=1};
    enum {CUTOFF=1};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ShiftSwitchingFunction();
    ShiftSwitchingFunction(Real cutoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ShiftSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// simple and fast test if we should apply the switching function
    bool roughTest(Real distSquared) const{return (distSquared <= myCutoff2);}
    Real cutoffSquared() const{return myCutoff2;}
    Real cutoff() const{return myCutoff;}
    void operator()(Real &value, Real &derivOverD, Real distSquared) const{
      if (distSquared > myCutoff2) {
	value=0.0;
	derivOverD=0.0;
	return;
      }
      Real x = (1.0 - distSquared*myCutoff_2);
      value=x*x;
      derivOverD=my4Cutoff_2*x;
    }

    static const std::string& getId() {return keyword;}
    void getParameters(std::vector<Parameter>& parameters) const;
    static unsigned int getParameterSize() {return 1;}
    static ShiftSwitchingFunction make(std::string& errMsg,std::vector<Value> values);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myCutoff;
    Real myCutoff2;
    Real myCutoff_2;
    Real my4Cutoff_2;
  };

  //______________________________________________________________________ INLINES

}
#endif /* SHIFTSWITCHINGFUNCTION_H */
