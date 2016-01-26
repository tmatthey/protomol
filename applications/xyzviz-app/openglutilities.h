/*  -*- c++ -*-  */
#ifndef OPENGLUTILITIES_H
#define OPENGLUTILITIES_H

#include <ostream>
#include "PPM.h"
#include "PGM.h"

namespace ProtoMol {

  unsigned int openglToEPS(std::ostream& output,void (*display)());
  unsigned int openglToPlain(std::ostream& output,void (*display)());
  void openglToPPM(PPM& ppm);
  void openglToPGM(PGM& pgm);

}
#endif
