/*  -*- c++ -*-  */
#ifndef ISGANGLESYSTEMFORCEBASE_H
#define ISGANGLESYSTEMFORCEBASE_H

#include <string>

namespace ProtoMol {
  //_________________________________________________________________ iSGAngleSystemForceBase

  class iSGAngleSystemForceBase {
  public:
    virtual ~iSGAngleSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ISGANGLESYSTEMFORCEBASE_H */

