/*  -*- c++ -*-  */
#ifndef ANGLESYSTEMFORCEBASE_H
#define ANGLESYSTEMFORCEBASE_H

#include <string>

namespace ProtoMol {
  //_________________________________________________________________ AngleSystemForceBase

  class AngleSystemForceBase {
  public:
    virtual ~AngleSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ANGLESYSTEMFORCEBASE_H */

