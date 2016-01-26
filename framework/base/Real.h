/*  -*- c++ -*-  */
#ifndef REAL_H
#define REAL_H

// The standard Real number type, settable by the user.
// 
// Define double as default type
#if !defined(USE_REAL_IS_FLOAT) && !defined(USE_REAL_IS_DOUBLE)
#define USE_REAL_IS_DOUBLE
#endif

namespace ProtoMol {


  /// The standard Real number type, settable by the user: double
#ifdef USE_REAL_IS_DOUBLE
  typedef double Real;
#endif

  /// The standard Real number type, settable by the user: float
#ifdef USE_REAL_IS_FLOAT
  typedef float Real;
#endif

}
#endif /*REAL_H*/
