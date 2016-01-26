/*  -*- c++ -*-  */
#ifndef ANM_H
#define ANM_H

#include "PDBReader.h"

namespace ProtoMol {

  class ANM {
    
  public:
    ANM(std::string pdbfilename, std::string outfilename, Real c, Real gm);
    ~ANM();
    void connectivity_matrix();
  private:
    PDBReader pdbReader;    
    Real cutoff;
    Real gamma;
    std::fstream myFile;
    std::string fname;
    
  };

}

#endif

