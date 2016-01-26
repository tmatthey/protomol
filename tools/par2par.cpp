#include "PARReader.h"
#include "PARWriter.h"
#include "Report.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {


  // Parse
  if(argc < 3)
    report << quit << "usage: " << argv[0] << " <PAR input file> <PAR output file>"<< endr;

  // Open
  PARReader in(argv[1]);
  if(!in)
    report << error << "Could not open '" << argv[1] << "'." << endr;

  // Read
  PAR par;
  if(!(in >> par))
    report << error << "Could not read par file '" << argv[1] << "'." << endr;

  // Open
  PARWriter out(argv[2]);
  if(!out)
    report << error << "Could not open '" << argv[2] << "'." << endr;

  if(!(out << par))
    report << error << "Could not write par file '" << argv[2] << "'." << endr;
  
  return 0;
}
