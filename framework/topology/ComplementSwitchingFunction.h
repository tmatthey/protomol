/* -*- c++ -*- */
#ifndef COMPLEMENTSWITCHINGFUNCTION_H
#define COMPLEMENTSWITCHINGFUNCTION_H

#include "Parameter.h"
#include "ComplementSwitchingFunctionBase.h"
namespace ProtoMol {
  //_________________________________________________________________ ComplementSwitchingFunction

  /**
   * Defines the complement of a given switching function (TOriginalSwitchingFunction).
   */
  template<class TOriginalSwitchingFunction>
  class ComplementSwitchingFunction : private ComplementSwitchingFunctionBase {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {USE=1};
    enum {MODIFY=1};
    enum {CUTOFF=0};

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ComplementSwitchingFunction(){}
    ComplementSwitchingFunction(const TOriginalSwitchingFunction sw) : myOrigFunc(sw){};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ComplementSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// simple and fast test if we should apply the switching function
    bool roughTest(Real /*distSquared*/) const{return true;}
    Real cutoffSquared() const{return Constant::MAXREAL;}
    Real cutoff() const{return Constant::MAXREAL;}
    void operator()(Real &value, Real &derivOverD, Real distSquared) const{
      myOrigFunc(value,derivOverD,distSquared);
      value      = 1.0-value;
      derivOverD = -derivOverD;
    }

    static std::string getId()  {return (keyword + TOriginalSwitchingFunction::getId());}

    void getParameters(std::vector<Parameter>& parameters) const{
      myOrigFunc.getParameters(parameters);
    }
    static unsigned int getParameterSize() {return TOriginalSwitchingFunction::getParameterSize();}

    static ComplementSwitchingFunction make(std::string& errMsg,std::vector<Value> values) {
      return ComplementSwitchingFunction(TOriginalSwitchingFunction::make(errMsg,values));
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TOriginalSwitchingFunction myOrigFunc;

  };

  //______________________________________________________________________ INLINES
}
#endif /* COMPLEMENTSWITCHINGFUNCTION_H */
