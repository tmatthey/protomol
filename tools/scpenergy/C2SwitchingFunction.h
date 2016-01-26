/* -*- c++ -*- */
#ifndef C2SWITCHINGFUNCTION_H
#define C2SWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "C2SwitchingFunctionBase.h"
#include "Matrix3by3.h"

namespace ProtoMol {
  //_________________________________________________________________ C2SwitchingFunction

  /**
   * The switching function provide C2 continues
   */
  class C2SwitchingFunction : private C2SwitchingFunctionBase {
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
    C2SwitchingFunction();
    C2SwitchingFunction(Real switchon,Real cutoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class C2SwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const{return (distSquared <= myCutoff2);}
    Real cutoffSquared() const{return myCutoff2;}
    Real cutoff() const{return myCutoff;}
    void operator()(Real &value, Real &deriv, Real distSquared) const{
      deriv=0.0; 
      if (distSquared > myCutoff2) {
	value=0.0;
      }
      else if (distSquared >= mySwitchon2) {
	Real c2 = myCutoff2-distSquared;
	Real c4 = c2*(mySwitch2 + 2.0*distSquared);
	value = mySwitch1*(c2*c4);
	deriv = mySwitch3*(c2*c2-c4);
      }
      else {
	value = 1.0;
      }
    }
    Matrix3by3 hessian(const Vector3D& rij, Real distSquared) const;

    static const std::string& getId() {return keyword;}
    void getParameters(std::vector<Parameter>& parameters) const;
    static unsigned int getParameterSize() { return 2;}
    static C2SwitchingFunction make(std::string& errMsg,std::vector<Value> values) ;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real mySwitchon, mySwitchon2,  myCutoff, myCutoff2, mySwitch1, mySwitch2, mySwitch3;
  };
  //______________________________________________________________________ INLINES

}
#endif /* C2SWITCHINGFUNCTION_H */
