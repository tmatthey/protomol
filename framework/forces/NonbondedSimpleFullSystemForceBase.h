/*  -*- c++ -*-  */
#ifndef NONBONDEDSIMPLEFULLSYSTEMFORCEBASE_H
#define NONBONDEDSIMPLEFULLSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ NonbondedSimpleFullSystemForceBase

  class NonbondedSimpleFullSystemForceBase {
  public:
    virtual ~NonbondedSimpleFullSystemForceBase(){}
    static const unsigned int defaultBlockSize;
    static const std::string keyword;
  };
}
#endif /* NONBONDEDSIMPLEFULLSYSTEMFORCEBASE_H */

