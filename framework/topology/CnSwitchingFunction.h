/* -*- c++ -*- */
#ifndef CnSWITCHINGFUNCTION_H
#define CnSWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "CnSwitchingFunctionBase.h"
#include "Matrix3by3.h"


namespace ProtoMol {
  //_________________________________________________________________ CnSwitchingFunction

  /**
   * The switching function provide C2, 3, 4 0r 6 continuous
   */
  class CnSwitchingFunction : private CnSwitchingFunctionBase {

#define	MAXEQNN	7
#define NUMSW	5

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
    CnSwitchingFunction();
    CnSwitchingFunction(Real switchon,Real cutoff,Real order,Real switchoff);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CnSwitchingFunction
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool roughTest(Real distSquared) const{return (distSquared <= myCutoff2);}
    Real cutoffSquared() const{return myCutoff2;}
    Real cutoff() const{return myCutoff;}
    void operator()(Real &value, Real &deriv, Real distSquared) const;
    Matrix3by3 hessian(const Vector3D& rij, Real distSquared) const;

    static const std::string& getId() {return keyword;}
    void getParameters(std::vector<Parameter>& parameters) const;
    static unsigned int getParameterSize() { return 4;}
    static CnSwitchingFunction make(std::string& errMsg,std::vector<Value> values) ;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
	  Real mySwitchon, mySwitchon2, myCutoff, myCutoff2, myOrder, mySwitchoff;
	  Real myIRange[MAXEQNN];
	  int ordInt, ordIdx;

    static Real swcoef[][MAXEQNN],
                dswcoef[][MAXEQNN],
                d2swcoef[][MAXEQNN];

  };
  //______________________________________________________________________ INLINES

}
#endif /* CnSWITCHINGFUNCTION_H */
