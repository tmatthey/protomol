/*  -*- c++ -*-  */
#ifndef ISGIMPROPERSYSTEMFORCEBASE_H
#define ISGIMPROPERSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ iSGImproperSystemForceBase

  class iSGImproperSystemForceBase {
  public:
    virtual ~iSGImproperSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ISGIMPROPERSYSTEMFORCEBASE_H */

  
