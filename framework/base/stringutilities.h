/*  -*- c++ -*-  */
#ifndef STRINGUTILITIES_H
#define STRINGUTILITIES_H

#include <string>
#ifdef HAVE_NO_SSTREAM
#include "sstream_local.h"
#else
#include <sstream>
#endif

#include "Vector3D.h"
#include <vector>

namespace ProtoMol {

  //_____________________________________________________________________ uppercase
  std::string uppercase (const std::string& word);
  //_____________________________________________________________________ lowercase
  std::string lowercase (const std::string& word);
  
  
  //_____________________________________________________________________ equal
  bool equal(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalNocase
  bool equalNocase(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalBegin
  bool equalBegin(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalBeginNocase
  bool equalBeginNocase(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalStart
  bool equalStart(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalStartNocase
  bool equalStartNocase(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalEnd
  bool equalEnd(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalEndNocase
  bool equalEndNocase(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalTerminate
  bool equalTerminate(const std::string& s1, const std::string& s2);
  //_____________________________________________________________________ equalTerminateNocase
  bool equalTerminateNocase(const std::string& s1, const std::string& s2);
  
  
  
  //_____________________________________________________________________ isReal
  bool isReal(const std::string& word);
  //_____________________________________________________________________ toReal
  bool toReal(const std::string& word, Real &r);
  //_____________________________________________________________________ toReal
  Real toReal(const std::string& word);
  //_____________________________________________________________________ isInt
  bool isInt(const std::string& word);
  //_____________________________________________________________________ toInt
  bool toInt(const std::string& word, int &i);
  //_____________________________________________________________________ toInt
  int  toInt(const std::string& word);
  //_____________________________________________________________________ isUInt
  bool isUInt(const std::string& word);
  //_____________________________________________________________________ toUInt
  bool toUInt(const std::string& word, unsigned int &i);
  //_____________________________________________________________________ toUInt
  unsigned int toUInt(const std::string& word);
  //_____________________________________________________________________ isBool
  bool isBool(const std::string& word);
  //_____________________________________________________________________ toBool
  bool toBool(const std::string& word, bool &b);
  //_____________________________________________________________________ toBool
  bool toBool(const std::string& word);
  //_____________________________________________________________________ isVector3D
  bool isVector3D(const std::string& word);
  //_____________________________________________________________________ toVector3D
  bool toVector3D(const std::string& word, Vector3D &c);
  //_____________________________________________________________________ toVector3D
  Vector3D toVector3D(const std::string& word);
  //_____________________________________________________________________ isVector
  bool isVector(const std::string& word);
  //_____________________________________________________________________ toVector
  std::vector<Real> toVector(const std::string& word);
  //_____________________________________________________________________ toVector
  bool toVector(const std::string& word, std::vector<Real> &c);
  //_____________________________________________________________________ toString
  std::string toString(Real x);
  //_____________________________________________________________________ toString
  std::string toString(Real x,unsigned int n, unsigned int m);
  //_____________________________________________________________________ toString
  std::string toString(bool x);
  //_____________________________________________________________________ toString
  std::string toString(const Vector3D& x);
  //_____________________________________________________________________ toString
  std::string toString(const std::vector<Real>& x);
  //_____________________________________________________________________ toString   
  /// NB: Needed for symmetry reason!
  inline const std::string & toString(const std::string& x) {return x;}
  //_____________________________________________________________________ toString !NB
  /// NB: Template to catch other types ...
  template <class T> 
  inline std::string toString(T x){
    // http://www.bespecific.com/dialog/becodetalk/archive/980405/0058.html
    std::stringstream ss;
    ss << x;
    return std::string(ss.str());
  }

  
  //_____________________________________________________________________ isBlank
  bool isBlank(const std::string& word);
  //_____________________________________________________________________ isblank
  bool isblankchar(char c);
  //_____________________________________________________________________ isPrintable
  bool isPrintable(const std::string& word);
  //_____________________________________________________________________ isprintablechar
  bool isprintablechar(char c);

  //_____________________________________________________________________ getBegin
  std::string getBegin(const std::string& s, std::string::size_type n);
  //_____________________________________________________________________ getEnd
  std::string getEnd(const std::string& s, std::string::size_type n);

  //_____________________________________________________________________ getRightFill
  std::string getRightFill(const std::string& s, std::string::size_type n);
  //_____________________________________________________________________ getLeftFill
  std::string getLeftFill(const std::string& s, std::string::size_type n);

  //_____________________________________________________________________ removeBeginEndBlanks
  std::string removeBeginEndBlanks(const std::string& s);

  //_____________________________________________________________________ ltstrNocase
  struct ltstrNocase {bool operator()(const std::string& s1, const std::string& s2) const;};
  //_____________________________________________________________________ ltstrNocaseOp
  bool ltstrNocaseOp (const std::string& s1, const std::string& s2);

  //_____________________________________________________________________ equalWildcard
  /**
   * Wildcard specifications:
   * - * : matches any string of characters (including none),
   * - \% : matches any single character,
   * - # : matches any string of digits (including none),
   * - + : matches any single digit.
   *
   * Return:
   * - 2 : match without wildcards
   * - 1 : match with wildcards
   * - 0 : no match at all
   */
  int equalWildcard(const std::string& wildcard, const std::string& name);

  //_____________________________________________________________________ splitString
  std::vector<std::string> splitString(const std::string& str);
  //_____________________________________________________________________ mergeString
  std::string mergeString(const std::vector<std::string>& str);
  //_____________________________________________________________________ normalizeString
  std::string normalizeString(const std::string& str);
  //_____________________________________________________________________ headString
  std::string headString(const std::string& str);
  //_____________________________________________________________________ tailString
  std::string tailString(const std::string& str);
}
#endif
  
