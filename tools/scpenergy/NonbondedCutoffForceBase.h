/*  -*- c++ -*-  */
#ifndef NONBONDEDCUTOFFFORCEBASE_H
#define NONBONDEDCUTOFFFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ NonbondedCutoffForceBase

  class NonbondedCutoffForceBase {
  public:
    virtual ~NonbondedCutoffForceBase(){}
    static const std::string keyword;
  };
}
#endif /* NONBONDEDCUTOFFFORCEBASE_H */

