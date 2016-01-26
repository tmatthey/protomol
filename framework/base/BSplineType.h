/*  -*- c++ -*-  */
#ifndef BSPLINETYPE_H
#define BSPLINETYPE_H

#include "AbstractEnumType.h"

namespace ProtoMol {
  //_____________________________________________________ BSplineEnum

  /// BSpline MOLLY types
  class BSplineEnum  {
  public:
    virtual ~BSplineEnum(){}
    enum Enum {
      FIRST = 0,       // Only internal purpose
      UNDEFINED  = 0,  // Value returned when no string matches
      SHORT,
      LONG,
      LAST             // Only internal purpose
    };
    static const std::string str[];
  };

  //_____________________________________________________ BSplineType

  typedef AbstractEnumType<BSplineEnum> BSplineType;
}
#endif //  BSPLINETYPE_H
