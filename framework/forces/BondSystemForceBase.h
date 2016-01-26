/*  -*- c++ -*-  */
#ifndef BONDSYSTEMFORCEBASE_H
#define BONDSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ BondSystemForceBase

  class BondSystemForceBase {
  public:
    virtual ~BondSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* BONDSYSTEMFORCEBASE_H */
