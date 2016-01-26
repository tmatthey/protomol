#include "PDBReader.h"
#include "PDBWriter.h"
#include "Report.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {


  // Pdbse
  if(argc < 3)
    report << quit << "usage: " << argv[0] << " <PDB input file> <PDB output file>"<< endr;

  // Open
  PDBReader in(argv[1]);
  if(!in)
    report << error << "Could not open '" << argv[1] << "'." << endr;

  // Read
  PDB pdb;
  if(!(in >> pdb))
    report << error << "Could not read pdb file '" << argv[1] << "'." << endr;

  // Open
  PDBWriter out(argv[2]);
  if(!out)
    report << error << "Could not open '" << argv[2] << "'." << endr;

  if(!(out << pdb))
    report << error << "Could not write pdb file '" << argv[2] << "'." << endr;
  
  return 0;
}
