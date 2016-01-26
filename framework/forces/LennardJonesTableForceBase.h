/*  -*- c++ -*-  */
#ifndef LENNARDJONESTABLEFORCEBASE_H
#define LENNARDJONESTABLEFORCEBASE_H

#include<string>
#include "mathutilities.h"
namespace ProtoMol {
  //_________________________________________________________________ LennardJonesTableForceBase

  class LennardJonesTableForceBase {
  public:
    //_________________________________________________________________ LookUpValues
    /**
     * Defines the function(s) to be tabulated.
     */
    class LookUpValues {
    public:
      enum {ENTRIES=2};
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class LookUpValues
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      template<typename TReal>
      void assign(Real r, Real r1, int ex, Real v, Real d, TReal* val) const{
	Real e =      (ex == -12 ? r : power<-12>(r1)); // value
	Real f = 12.0*(ex == -14 ? r : power<-14>(r1)); // gradient
	val[0] =        e*v;
	val[1] = -0.5 * (f*v-e*d);
	e =     -(ex == -6 ? r : power<-6>(r1)); // value
	f = -6.0*(ex == -8 ? r : power<-8>(r1)); // gradient
	val[4] =        e*v;
	val[5] = -0.5 * (f*v-e*d);
      }
    };



  public:
    static const std::string keyword;
  };
}
#endif /* LENNARDJONESTABLEFORCEBASE_H */

