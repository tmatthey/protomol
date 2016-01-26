#include "DCDTrajectoryReader.h"
#include "DCDTrajectoryWriter.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "XYZBinReader.h"
#include "XYZBinWriter.h"
#include "XYZReader.h"
#include "XYZTrajectoryReader.h"
#include "XYZTrajectoryWriter.h"
#include "XYZWriter.h"
#include "XFigWriter.h"
#include "stringutilities.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::vector;
using std::string;

static bool cmpSpherical (const Vector3D& v1, const Vector3D& v2);
//_____________________________________________________________________ cmpSpherical
bool cmpSpherical (const Vector3D& v1, const Vector3D& v2){
  return (v1.normSquared() < v2.normSquared());
}

//_____________________________________________________________________ dcd2dcd
struct Format{
  enum Type {UNDEF,
	     XFIG,
	     XYZ,
	     PDB,
	     DCD,
	     XYZBIN,
	     XYZTRA};
};

int main(int argc, char **argv) {


  // Parse
  if(argc < 3)
    report << quit << "usage: " << argv[0] << " [-xyz] [-xfig] [-xyzbin] [-pdb] [-dcd] [-xyztra] [-real4] [-real8] [-bigEndian] [-littleEndian] [-range <from> <to>] [-scale <float>] [-translate <float> <float><float>] <coordinate input file> ... <coordinate input file> <coordinate output file>"<< endr;

  vector<string> in;
  string out(argv[argc-1]);
  Format::Type format = Format::UNDEF;
  bool real8  = false;
  bool real4  = false;
  bool opt    = true;
  bool high   = false;
  bool low    = false;
  bool check  = false;
  Real scale  = 1.0;
  Vector3D delta(0.0,0.0,0.0);
  int from    = 0;
  int to      = Constant::MAX_INT;
  for(int i=1;i<argc-1;++i){
    string arg(argv[i]);
    if(opt){
      if(arg == "-xyzbin")
	format = Format::XYZBIN;
      else if(arg == "-xfig")
	format = Format::XFIG;
      else if(arg == "-xyz")
	format = Format::XYZ;
      else if(arg == "-pdb")
	format = Format::PDB;
      else if(arg == "-dcd")
	format = Format::DCD;
      else if(arg == "-xyztra")
	format = Format::XYZTRA;
      else if(arg == "-real4")
	real4 = true;
      else if(arg == "-real8")
	real8 = true;
      else if(arg == "-bigEndian")
	high = true;
      else if(arg == "-littleEndian")
	low = true;
      else if(arg == "-check"){
	check = true;
	from = Constant::MAX_INT;
	to = Constant::MAX_INT;
      }
      else if(arg == "-range" && i<argc-3){
	i++;
	from = atoi(argv[i]);
	i++;
	to = atoi(argv[i]);
      }
      else if(arg == "-translate" && i<argc-4){
	i++;
	delta.x = atoi(argv[i]);
	i++;
	delta.y = atoi(argv[i]);
	i++;
	delta.z = atoi(argv[i]);
      }
      else if(arg == "-scale"){
	i++;
	scale = atof(argv[i]);
      }
      else
	opt = false;
    }
    if(!opt)
      in.push_back(arg);
  }
  if(check){
      in.push_back(out);
      out = "";    
  }
  

  // Read frames
  vector<XYZ> trajectory;
  XYZ xyz;
  PDB pdb;
  vector<string> names;
  int count=0;

  for(unsigned int i=0;i<in.size();++i){
    if(PDBReader(in[i]).tryFormat()){
      PDBReader reader(in[i]);
      PDB tmp;
      if(reader >> tmp){
	xyz.clear();
	xyz.coords = tmp.coords;
	for(unsigned int j=0;j<tmp.atoms.size();++j)
	  xyz.names.push_back(tmp.atoms[j].elementName);			
	names = xyz.names;
	pdb = tmp;
	if(count >= from && count < to){
	  trajectory.push_back(xyz);
	}
      }
      count++;
    }
    else if(XYZReader(in[i]).tryFormat()){
      XYZReader reader(in[i]);
      if(reader >> xyz){
	names = xyz.names;
	if(count >= from && count < to){
	  trajectory.push_back(xyz);
	}
      }
      count++;
    }
    else if(XYZBinReader(in[i]).tryFormat()){
      if(count >= from && count < to){
	XYZBinReader reader(in[i]);
	if(reader >> xyz){
	  trajectory.push_back(xyz);
	}
      }
      count++;
    }
    else if(DCDTrajectoryReader(in[i]).tryFormat()){
      DCDTrajectoryReader reader(in[i]);
      while((reader >> xyz)){
	if(count >= from && count < to){
	  trajectory.push_back(xyz);
	}
	count++;
      }
    }
    else if(XYZTrajectoryReader(in[i]).tryFormat()){
      XYZTrajectoryReader reader(in[i]);
      while((reader >> xyz)){
	if(count >= from && count < to){
	  trajectory.push_back(xyz);
	  names = xyz.names;
	}
	count++;
      }
    }
    else if(!isAccessible(in[i])){
      report << "Can not open \'"<<in[i]<<"\', skipping."<<endr;
    }
    else {
      report << "Can not figure out format of \'"<<in[i]<<"\', skipping."<<endr;
    }
  }

  //Check & statistic
  if(check){
    report << quit << "Read "<<count<<" frame(s). Nothing to write back.";
    report <<endr;
  }

  report << "Read "<<count<<" frame(s), writing back "<<trajectory.size()<<" frame(s) ("<<from<<","<<to<<").";
  if(trajectory.size()<1)
    report << quit << " Nothing to write back.";
  report <<endr;

  if(format == Format::UNDEF){
    format = (trajectory.size()<2?Format::XYZ:Format::XYZTRA);
  }

  if(scale != 1.0 || delta.norm()> 0){
    for(unsigned int i=0;i<trajectory.size();++i){
      for(unsigned int j=0;j<trajectory[i].size();++j){
	trajectory[i].coords[j] += delta;
	trajectory[i].coords[j].x *= scale;    
	trajectory[i].coords[j].y *= scale;    
	trajectory[i].coords[j].z *= scale;    
      }
    }
  }

  if(trajectory.size() > 0  && names.size() == 0){
    names = vector<string>(trajectory[0].names.size(),"X");
  }

  // Overwrite names with the last one read
  for(unsigned int i=0;i<trajectory.size();++i){
    trajectory[i].names = names;
  }
  
  bool multiple = trajectory.size()> 1;

  // Write
  switch (format){
  case Format::XYZ:
    {
      if(multiple){
	for(unsigned int i=0;i<trajectory.size();++i){
	  XYZWriter writer(out+"."+toString(i));
	  if(writer << trajectory[i])
	    report << "Wrote XYZ file \'"<<writer.getFilename()<<"\'"<<endr;
	  else
	    report << error << "Could not write XYZ file \'"<<writer.getFilename()<<"\'"<<endr;	    
	}	
      }
      else{
	XYZWriter writer(out);
	if(writer << trajectory[0])
	  report << "Wrote XYZ file \'"<<writer.getFilename()<<"\'"<<endr;
	else
	  report << error << "Could not write XYZ file \'"<<writer.getFilename()<<"\'"<<endr;	    
      }

      break;
    }
  case Format::XFIG:
    {
      if(multiple){
	for(unsigned int i=0;i<trajectory.size();++i){
	  XFigWriter writer(out+"."+toString(i));
	  Vector3D a,b;
	  trajectory[i].coords.boundingbox(a,b);
	  if(equal(a.z,b.z)){
	    std::sort(trajectory[0].coords.begin(),trajectory[0].coords.end(),cmpSpherical);
	    Matrix3by3 mat;
	    Real beta = atan2(trajectory[0].coords[0].x,trajectory[0].coords[0].y);
	    report << "XFig z-rotation of "<<rtod(beta)<<"."<<endr;
	    mat.rotate(Vector3D(0.0,0.0,1.0),-beta);
	    writer.setTransformation(mat);
	  }
	  if(writer << trajectory[i])
	    report << "Wrote XFig file \'"<<writer.getFilename()<<"\'"<<endr;
	  else
	    report << error << "Could not write XFig file \'"<<writer.getFilename()<<"\'"<<endr;	    
	}	
      }
      else{
	XFigWriter writer(out);
	Vector3D a,b;
	trajectory[0].coords.boundingbox(a,b);
	if(equal(a.z,b.z)){
	  std::sort(trajectory[0].coords.begin(),trajectory[0].coords.end(),cmpSpherical);
	  Matrix3by3 mat;
	  Real beta = atan2(trajectory[0].coords[0].x,trajectory[0].coords[0].y);
	  report << "XFig z-rotation of "<<rtod(beta)<<"."<<endr;
	  mat.rotate(Vector3D(0.0,0.0,1.0),-beta);
	  writer.setTransformation(mat);
	}
	if(writer << trajectory[0])
	  report << "Wrote XFig file \'"<<writer.getFilename()<<"\'"<<endr;
	else
	  report << error << "Could not write XFig file \'"<<writer.getFilename()<<"\'"<<endr;	    
      }

      break;
    }
  case Format::PDB:
    {

      if(multiple){
	for(unsigned int i=0;i<trajectory.size();++i){
	  PDBWriter writer(out+"."+toString(i));
	  pdb.coords = trajectory[i].coords;
	  if(writer << pdb)
	    report << "Wrote PDB file \'"<<writer.getFilename()<<"\'"<<endr;
	  else
	    report << error << "Could not write PDB file \'"<<writer.getFilename()<<"\'"<<endr;	    
	}	
      }
      else{
	PDBWriter writer(out);
	pdb.coords = trajectory[0].coords;
	if(writer << pdb)
	  report << "Wrote PDB file \'"<<writer.getFilename()<<"\'"<<endr;
	else
	  report << error << "Could not write PDB file \'"<<writer.getFilename()<<"\'"<<endr;
	    
      }

      break;
    }
  case Format::DCD:
    {
      DCDTrajectoryWriter writer(out);
      if(!writer)
	report << error << "Could not open '" << writer.getFilename() << "'." << endr;
      if(low)
	writer.setLittleEndian(true);
      else if(high)
	writer.setLittleEndian(false);
      for(unsigned int i=0;i<trajectory.size();++i)
	writer << trajectory[i];
      report << "Wrote DCD file \'"<<writer.getFilename()<<"\'"<<endr;
      break;
    }
  case Format::XYZBIN:
    {

      if(multiple){
	for(unsigned int i=0;i<trajectory.size();++i){
	  XYZBinWriter writer(out+"."+toString(i));
	  if(low)
	    writer.setLittleEndian(true);
	  else if(high)
	    writer.setLittleEndian(false);
	  if(real4)
	    writer.setSize(4);
	  else if(real8)
	    writer.setSize(8);	      
	  if(writer << trajectory[i])
	    report << "Wrote XYZBin file \'"<<writer.getFilename()<<"\'"<<endr;
	  else
	    report << error << "Could not write XYZBin file \'"<<writer.getFilename()<<"\'"<<endr;	    
	}	
      }
      else{
	XYZBinWriter writer(out);
	if(low)
	  writer.setLittleEndian(true);
	else if(high)
	  writer.setLittleEndian(false);
	if(real4)
	  writer.setSize(4);
	else if(real8)
	  writer.setSize(8);	      
	if(writer << trajectory[0])
	  report << "Wrote XYZBin file \'"<<writer.getFilename()<<"\'"<<endr;
	else
	  report << error << "Could not write XYZBin file \'"<<writer.getFilename()<<"\'"<<endr;	    
      }

      break;
    }
  case Format::XYZTRA:
    {

      XYZTrajectoryWriter writer(out);
      if(!writer)
	report << error << "Could not open '" << writer.getFilename() << "'." << endr;
      for(unsigned int i=0;i<trajectory.size();++i)
	writer << trajectory[i];
      report << "Wrote XYZTrajectory file \'"<<writer.getFilename()<<"\'"<<endr;

      break;
    }


  default:
    report << "Ops undefined output format!"<<endr;
  }

  return 0;
}
