/*  -*- c++ -*-  */
#ifndef GNM_H
#define GNM_H

#include "PDBReader.h"

namespace ProtoMol {

  class GNM {
    
  public:
    GNM(std::string pdbfilename, std::string outfilename, Real c, Real gm);
    ~GNM();
    void connectivity_matrix();
  private:
    PDBReader pdbReader;
    Real cutoff;
    Real gamma;
    std::string fname;
    std::fstream myFile;	
	std::fstream acFile; /* file which contains the coordinates of the alpha carbons */
    
  };

}

#endif
