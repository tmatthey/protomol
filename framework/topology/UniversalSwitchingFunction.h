/* -*- c++ -*- */
#ifndef UNIVERSALSWITCHINGFUNCTION_H
#define UNIVERSALSWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "UniversalSwitchingFunctionBase.h"
#include "pmconstants.h"

namespace ProtoMol {
  //_________________________________________________________________ UniversalSwitchingFunction

  /**
   * Switching function that always returns 1 and 0 for the derivative, acts 
   * as a placeholder. Is intend to be use when non modification and / or truncation
   * of the potential is required. Defines USE = 0.
   */
  class UniversalSwitchingFunction : private UniversalSwitchingFunctionBase {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// defines a dummy switching function
    enum {USE=0};
    /// no modification
    enum {MODIFY=0};
    /// no cutoff
    enum {CUTOFF=0};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class UniversalSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// simple and fast test if we should apply the switching function
    bool roughTest(Real /*distSquared*/) const {return true;};
    Real cutoffSquared() const{return Constant::MAXREAL;}
    Real cutoff() const{return Constant::MAXREAL;}
    void operator()(Real &value, Real &deriv, Real /*distSquared*/) const{ value=1.0; deriv=0.0;}

    static const std::string& getId() {return keyword;}
    void getParameters(std::vector<Parameter>&) const{};
    static unsigned int getParameterSize() {return 0;}
    static UniversalSwitchingFunction make(std::string&, std::vector<Value>) {
      return UniversalSwitchingFunction();
    }

  };
}
//______________________________________________________________________ INLINES

#endif /* UNIVERSALSWITCHINGFUNCTION_H */
