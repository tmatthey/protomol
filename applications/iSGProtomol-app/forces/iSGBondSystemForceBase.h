/*  -*- c++ -*-  */
#ifndef ISGBONDSYSTEMFORCEBASE_H
#define ISGBONDSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ iSGBondSystemForceBase

  class iSGBondSystemForceBase {
  public:
    virtual ~iSGBondSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ISGBONDSYSTEMFORCEBASE_H */
