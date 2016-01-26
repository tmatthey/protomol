/*  -*- c++ -*-  */
#ifndef ELECTRICFIELDFORCEBASE_H
#define ELECTRICFIELDFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ ElectricFieldSystemForceBase

  class ElectricFieldSystemForceBase {
  public:
    virtual ~ElectricFieldSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ELECTRICFIELDFORCEBASE_H */

