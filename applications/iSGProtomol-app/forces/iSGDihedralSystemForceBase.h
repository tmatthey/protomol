/*  -*- c++ -*-  */
#ifndef ISGDIHEDRALSYSTEMFORCEBASE_H
#define ISGDIHEDRALSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol {
  //_________________________________________________________________ iSGDihedralSystemForceBase

  class iSGDihedralSystemForceBase {
  public:
    virtual ~iSGDihedralSystemForceBase(){}
    static const std::string keyword;
  };
}
#endif /* ISGDIHEDRALSYSTEMFORCEBASE_H */

