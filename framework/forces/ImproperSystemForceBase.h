/*  -*- c++ -*-  */
#ifndef IMPROPERSYSTEMFORCEBASE_H
#define IMPROPERSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ ImproperSystemForceBase

  class ImproperSystemForceBase {
  public:
    virtual ~ImproperSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* IMPROPERSYSTEMFORCEBASE_H */

  
