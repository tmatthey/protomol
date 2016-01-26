#include "DCDTrajectoryReader.h"
#include "DCDTrajectoryWriter.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {


  // Parse
  if(argc < 3)
    report << quit << "usage: " << argv[0] << " <dcd input file> <dcd output file>"<< endr;

  // Open
  DCDTrajectoryReader in(argv[1]);
  if(!in)
    report << error << "Could not open '" << argv[1] << "'." << endr;

  // Read
  std::vector<Vector3DBlock> trajectory;
  Vector3DBlock xyz;
  while((in >> xyz))
    trajectory.push_back(xyz);

  //Process
  if(trajectory.size() < 1)
    report << error << "Did not find any frames." << endr;
  Vector3D a0,b0;
  trajectory[0].boundingbox(a0,b0);
  for(unsigned int i=1;i<trajectory.size();++i){
    Vector3D a,b;
    trajectory[i].boundingbox(a,b);
    Vector3D d(a0-a);
    for(unsigned int j=0;j<trajectory[i].size();++j)
      trajectory[i][j] += d;
  }

  // Open
  DCDTrajectoryWriter out(argv[2]);
  if(!out)
    report << error << "Could not open '" << argv[2] << "'." << endr;
  for(unsigned int i=0;i<trajectory.size();++i)
    out << trajectory[i];

  // Done
  report << hint << "Wrote DCD with " << trajectory.size() << " frame(s) and " << trajectory[0].size() << " element(s)." << endr;  
  
  return 0;
}
