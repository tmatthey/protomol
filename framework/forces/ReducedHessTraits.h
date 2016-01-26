/*  -*- c++ -*-  */
#ifndef REDUCEDHESSTRAITS_H
#define REDUCEDHESSTRAITS_H

#include "ReducedHessLennardJones.h"
#include "ReducedHessCoulomb.h"
#include "LennardJonesForce.h"
#include "CoulombForce.h"

namespace ProtoMol {

  // Trait struct 
  template<typename T> 
  struct ReducedHessTraits;

  template<>                     
  struct ReducedHessTraits<LennardJonesForce>{
    typedef ReducedHessLennardJones Hessian;

  }; 

  template<>                     
  struct ReducedHessTraits<CoulombForce>{
    typedef ReducedHessCoulomb Hessian;
  }; 


}
#endif /* REDUCEDHESSTRAITS_H */
