/*  -*- c++ -*-  */
#ifndef TRAJECTORYANALYZER_H
#define TRAJECTORYANALYZER_H


#include "PDBReader.h"
#include "DCDTrajectoryReader.h"

namespace ProtoMol {

  //__________________________________TrajectoryAnalyzer
  //A simple class which contains the code to extract data
  //from MD trajectory for analysis to obtain dynamics 
  //information
  class TrajectoryAnalyzer {
  public:
    TrajectoryAnalyzer(std::string dcdfilenames, std::string pdbfilename, std::string out_filename);
    ~TrajectoryAnalyzer();

    void bond_time_correlation(int natom, int hatom, int records, int frequency);
    void collective_mode_analysis(int records);	
    
    
  private:
	void getDcdFiles();
	void copyXYZ(XYZ& to, XYZ from);
	
    //DCDTrajectoryReader dcdReader;
    //PDBReader pdbReader;
    //PDB pdb;
    std::string fname;
	std::string dcdFileSetName;
	std::fstream dcdFileSet;
    std::fstream myFile;
	std::vector<std::string> dcdFiles;	
    DCDTrajectoryReader dcdReader;
    PDBReader pdbReader;
    PDB pdb;
    
  };

}

#endif /*TRAJECTORYANALYZER_H*/
