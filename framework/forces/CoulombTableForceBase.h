/*  -*- c++ -*-  */
#ifndef COULOMBTABLEFORCEBASE_H
#define COULOMBTABLEFORCEBASE_H

#include<string>
#include "mathutilities.h"
namespace ProtoMol {
  //_________________________________________________________________ CoulombTableForceBase

  class CoulombTableForceBase {
  public:
    //_________________________________________________________________ LookUpValues
    /**
     * Defines the function(s) to be tabulated.
     */
    class LookUpValues {
    public:
      enum {ENTRIES=1};
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class LookUpValues
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      template<typename TReal>
      void assign(Real r, Real r1, int ex, Real v, Real d, TReal* val) const{
	Real e = (ex == -1 ? r : power<-1>(r1)); // value
	Real f = (ex == -3 ? r : power<-3>(r1)); // gradient
	val[0] =        e*v;
	val[1] = -0.5 * (f*v-e*d);
      }
    };

  public:
    static const std::string keyword;
  };
}
#endif /* COULOMBTABLEFORCEBASE_H */

