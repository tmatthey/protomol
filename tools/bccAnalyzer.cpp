#include "Array.h"
#include "Matrix3by3.h"
#include "Vector3DBlock.h"
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
#include "stringutilities.h"
#include "mathutilities.h"
#include <algorithm>
#include <set>
using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::vector;
using std::string;
using std::sort;
using std::set;


//_____________________________________________________________________ Vector3DIdx
class Vector3DIdx {
public:
  Vector3D v;  
  int i;
  Vector3DIdx():v(Vector3D(0.0,0.0,0.0)),i(-1){}
  //Vector3DIdx(Vector3DIdx a):v(a.v),i(a.i){}
  Vector3DIdx(const Vector3D& a):v(a),i(-1){}
  Vector3DIdx(const Vector3D& a, int b):v(a),i(b){}
  //operator Vector3D() const{return v;}
};

//_____________________________________________________________________ radialCmpVector3D
class radialCmpVector3D {
private:
  const Vector3D c;
public:
  radialCmpVector3D():c(Vector3D(0.0,0.0,0.0)){}
  radialCmpVector3D(const Vector3D& a):c(a){}
  bool operator()(const Vector3D& v1, const Vector3D& v2) const{
    return ((v1-c).normSquared() < (v2-c).normSquared() );
  }
  bool operator()(const Vector3DIdx& v1, const Vector3DIdx& v2) const{
    return ((v1.v-c).normSquared() < (v2.v-c).normSquared() );
  }

};

//_____________________________________________________________________ PairVector3D
class PairVector3D {

public:
  Vector3DIdx v0,c,v1;
  PairVector3D(){}
  PairVector3D(const Vector3DIdx& a, const Vector3DIdx& b, const Vector3DIdx& c):v0(a),c(b),v1(c){}
  bool operator()(const PairVector3D& a, const PairVector3D& b) const {
    Real d = ((a.v0.v+a.v1.v)*0.5-a.c.v).normSquared();
    Real e = ((b.v0.v+b.v1.v)*0.5-b.c.v).normSquared();
    if(fabs((d-e)/std::max(d,e)) <1e-9)
      return (a.v0.v-a.v1.v).normSquared() < (b.v0.v-b.v1.v).normSquared();
    return (d < e);
  }
};

//_____________________________________________________________________ dcd2dcd
struct Format{
  enum Type {UNDEF,
	     XYZ,
	     PDB,
	     DCD,
	     XYZBIN,
	     XYZTRA};
};

static bool isBCC(XYZ xyz, Real eps, Real rest, Vector3D& normal, Real& d, Real& err, Real& minimal, Real& maximal,
		  int& j0, int& j1, int& j2, int& j3, int& j4, int& j5, Vector3DIdx& origin, Vector3DIdx& center){

  Vector3D perfect(rtod(acos(1.0/sqrt(3.0))),
		   rtod(acos(1.0/sqrt(3.0))),
		   rtod(acos(1.0/     3.0 )));

  Vector3D sum;
  vector<Vector3DIdx> coords(xyz.coords.size());
  for(unsigned int i=0;i<coords.size();i++){
    coords[i].v = xyz.coords[i];
    sum += coords[i].v/coords.size();
    coords[i].i = i;
  }
  sort(coords.begin(),coords.end(),radialCmpVector3D(sum));
  center = coords[0];

  sort(coords.begin(),coords.end(),radialCmpVector3D());
  origin = coords[0];

  for(unsigned int i=0;i<coords.size();i++){
    coords[i].v-=origin.v;
    xyz.coords[i]-=origin.v;
  }
  //xyz.names[origin.i] = "O";


  sort(coords.begin(),coords.end(),radialCmpVector3D());

  minimal = coords[1].v.norm();
  Real dmax = coords[1].v.normSquared()*power<2>(eps);
  maximal = coords[coords.size()-1].v.norm();

  vector<Vector3DIdx> inner;

  for(unsigned int i=1;i<coords.size();i++){
    if(coords[i].v.normSquared()> dmax)
      break;
    inner.push_back(coords[i]);
  }
       
  report << hint <<inner.size()<<" particle(s) within "<<sqrt(dmax)<<" with epsilon "<<eps<<". ";


  //  for(unsigned int i=0;i<inner.size();i++)
  //  xyz.names[inner[i].i] = "Q";

  vector<PairVector3D> all;
  for(unsigned int i=0;i<inner.size();i++){
    for(unsigned int j=i+1;j<inner.size();j++){
      all.push_back(PairVector3D(inner[i],Vector3D(0.0,0.0,0.0),inner[j]));
    }
  }
  sort(all.begin(),all.end(),PairVector3D());
  report <<all.size()<<" pairs(s). ";

  for(unsigned int i=0;i<all.size();i++){
    if((all[i].v0.v+all[i].v1.v).norm()/(2*minimal) >rest){
      all.resize(i);
      break;
    }

  }
    
  //all.resize(std::min((int)all.size(),(argc >3 ? atoi(argv[4]):4)));
  report <<"Keeping "<<all.size()<<" pairs(s), "<<(all[all.size()-1].v0.v+all[all.size()-1].v1.v).norm()/(2*minimal)<<"."<<endr;
  Real z = Constant::MAXREAL;
  j0 =-1;
  j1 =-1;
  j2 =-1;
  j3 =-1;
  j4 =-1;
  j5 =-1;
  Vector3DBlock points(7);
  points[0] = Vector3D(0.0,0.0,0.0);
  for(unsigned int i=0;i<all.size();i++){
    points[1] = all[i].v0.v;
    points[2] = all[i].v1.v;
    for(unsigned int j=i+1;j<all.size();j++){
      points[3] = all[j].v0.v;
      points[4] = all[j].v1.v;
      for(unsigned int k=j+1;k<all.size();k++){
	set<int> index;
	index.insert(all[i].v0.i);
	index.insert(all[i].v1.i);
	index.insert(all[j].v0.i);
	index.insert(all[j].v1.i);
	index.insert(all[k].v0.i);
	index.insert(all[k].v1.i);
	if(index.size() <6)
	  continue;
	points[5] = all[k].v0.v;
	points[6] = all[k].v1.v;
	Vector3D n;
	Real e,dc;


	if(!points.fitplane(n,dc,e))
	  continue;

	Vector3D a(all[i].v0.v - all[i].v1.v);
	Vector3D b(all[j].v0.v - all[j].v1.v);
	Vector3D c(all[k].v0.v - all[k].v1.v);
	a.normalize();
	b.normalize();
	c.normalize();
	vector<Real> ang(3);
	ang[0] = rtod(acos(std::min(fabs(a.dot(b)),1.0)));
	ang[1] = rtod(acos(std::min(fabs(a.dot(c)),1.0)));
	ang[2] = rtod(acos(std::min(fabs(b.dot(c)),1.0)));
	sort(ang.begin(),ang.end());
	Vector3D angles(ang[0],ang[1],ang[2]);

	if((angles-perfect).norm() < z && e/minimal < 0.1){
	  z = (angles-perfect).norm();
	  err = e;
	  normal = n;
	  d = dc;
	  j0 = all[i].v0.i;
	  j1 = all[i].v1.i;
	  j2 = all[j].v0.i;
	  j3 = all[j].v1.i;
	  j4 = all[k].v0.i;
	  j5 = all[k].v1.i;
	  // 	  report << "A=[0 0 0;"
	  // 		 <<all[i].v0.v.x/minimal<<" "<<all[i].v0.v.y/minimal<<" "<<all[i].v0.v.z/minimal<<";"
	  // 		 <<all[i].v1.v.x/minimal<<" "<<all[i].v1.v.y/minimal<<" "<<all[i].v1.v.z/minimal<<";"
	  // 		 <<all[j].v0.v.x/minimal<<" "<<all[j].v0.v.y/minimal<<" "<<all[j].v0.v.z/minimal<<";"
	  // 		 <<all[j].v1.v.x/minimal<<" "<<all[j].v1.v.y/minimal<<" "<<all[j].v1.v.z/minimal<<";"
	  // 		 <<all[k].v0.v.x/minimal<<" "<<all[k].v0.v.y/minimal<<" "<<all[k].v0.v.z/minimal<<";"
	  // 		 <<all[k].v1.v.x/minimal<<" "<<all[k].v1.v.y/minimal<<" "<<all[k].v1.v.z/minimal<<"];C = [A ones(size(A,1),1)];[u d v] = svd(C);n= v(1:3,4);d=v(4,4)/norm(n);n=n/norm(n); r=norm(A*n+d)"<<endr;
	  //	  report << hint << "err="<<err/minimal <<", n="<< n <<", d="<<dc<<endr;
	  report << hint << "err="<<e/minimal <<", n="<< n <<", d="<<dc<<", -rot "<<toString(n.x)<<" "<<toString(n.y)<<" "<<toString(n.z)
		 << "  "<<angles<<endr;
	  
	}
      }
    }
  }

  return (j0 >=0);

}

int main(int argc, char **argv) {


  // parse
  if (argc < 2 || (argc >= 2 && (string(argv[1]) =="-h" ||
                                 string(argv[1]) =="--help" ))){
    report << plain << quit << "usage: "<<argv[0]<<"-n <frame=0> -eps <distance epsilon=1.6> -r <normal distance epsilon=0.5> input output "<<endr;
  }


  int cur = 1;
  int frame = 0;
  Real eps  = 1.6;
  Real rest = 0.5;
  bool doAll = false;
  while (cur<(argc-1) && argv[cur][0]=='-') {

    string str(argv[cur]);

    if (str == "-n") {
      frame = toInt(argv[++cur]);
      cur++;
      continue;
    }
    if (str == "-eps") {
      eps = toReal(argv[++cur]);
      cur++;
      continue;
    }
    if (str == "-all") {
      doAll = true;
      cur++;
      continue;
    }
    if (str == "-r") {
      rest = toReal(argv[++cur]);
      cur++;
      continue;
    }

    break;
  }

  string in;
  string out;
  if(cur < argc)
    in = argv[cur++];
  if(cur < argc)
    out = argv[cur];
  if(out.empty())
    report << hint << "Missing output file."<<endr;
  if(in.empty())
    report << error << "Missing input file."<<endr;
	    

  // read frames
  vector<XYZ> trajectory;
  XYZ xyz;

  if(PDBReader(in).tryFormat()){
    PDBReader reader(in);
    reader >> xyz;
    trajectory.push_back(xyz);
  }
  else if(XYZReader(in).tryFormat()){
    XYZReader reader(in);
    reader >> xyz;
    trajectory.push_back(xyz);
  }
  else if(XYZBinReader(in).tryFormat()){
    XYZBinReader reader(in);
    reader >> xyz;
    trajectory.push_back(xyz);
  }
  else if(DCDTrajectoryReader(in).tryFormat()){
    DCDTrajectoryReader reader(in);
    while(reader >> xyz)
      trajectory.push_back(xyz);
  }
  else if(XYZTrajectoryReader(in).tryFormat()){
    XYZTrajectoryReader reader(in);
    while(reader >> xyz)
      trajectory.push_back(xyz);
  }
  else if(!isAccessible(in)){
    report << error<< "Can not open \'"<<in<<"\'."<<endr;
  }
  else {
    report << error<< "Can not figure out format of \'"<<in<<"\', skipping."<<endr;
  }

  // do all?
  int j0,j1,j2,j3,j4,j5;
  Vector3D    normal;
  Vector3DIdx origin, center;
  Real err,d,minimal,maximal;
  vector<Real> errors(trajectory.size(),Constant::MAXREAL);
  if(doAll){
    report << donthint;
    for(unsigned int i=0;i<trajectory.size();i++){
      //      report << "Frame : \t"<<i<<endr;
      xyz = trajectory[i];
      isBCC(xyz,eps, rest,normal,d,err,minimal,maximal,j0, j1, j2, j3, j4, j5, origin,center);
//	report << recoverable << "Nope, could not find BCC structure!"<<endr;
//       report << plain 
// 	     << "e="<<err/minimal<<"\n"
// 	     << "p="<<toString(normal.x)<<","<<toString(normal.y)<<","<<toString(normal.z)
// 	     <<","<<d<<"\n"
// 	     << "m="<<(xyz.coords[j0]+xyz.coords[j1]-origin.v).norm()/(2*minimal)<<","
// 	     <<(xyz.coords[j2]+xyz.coords[j3]-origin.v).norm()/(2*minimal)<<","
// 	     <<(xyz.coords[j4]+xyz.coords[j5]-origin.v).norm()/(2*minimal)<<"\n"
// 	     << "d="<<(xyz.coords[j0]-origin.v).norm()/minimal<<","
// 	     << (xyz.coords[j1]-origin.v).norm()/minimal<<","
// 	     << (xyz.coords[j2]-origin.v).norm()/minimal<<","
// 	     << (xyz.coords[j3]-origin.v).norm()/minimal<<","
// 	     << (xyz.coords[j4]-origin.v).norm()/minimal<<","
// 	     << (xyz.coords[j5]-origin.v).norm()/minimal<<"\n"
// 	     << "i="<<j0<<","
// 	     << j1<<","
// 	     << j2<<","
// 	     << j3<<","
// 	     << j4<<","
// 	     << j5<<endr;
      
      errors[i] = err/minimal;
    }    
    report << dohint;
    sort(errors.begin(),errors.end());
    Real last = 0.0;
    report << "Error:";
    for(unsigned int i=0;i<errors.size();i++){
      if(errors[i] >= 0.005 && last < 0.005){
	report << " 0.005:"<<i;
      }
      if(errors[i] >= 0.01 && last < 0.01){
	report << " 0.01:"<<i;

      }
      if(errors[i] >= 0.05 && last < 0.05){
	report << " 0.05:"<<i;

      }
      if(errors[i] >= 0.1 && last < 0.1){
	report << " 0.1:"<<i;
      }
      if(errors[i] >= 1 && last < 1){
	report << " 1:"<<i;
      }
      last = errors[i];
    }
    report << endr;
  }
  
  // select frame
  if(frame < 0 || frame >= static_cast<int>(trajectory.size()))
    report << error<< "Frame "<<frame<<"< out of range [0,"<<trajectory.size()<<"."<<endr;
  report << hint << "Read "<<xyz.coords.size()<<" particle(s) from frame "<<frame<<"."<<endr;
  xyz = trajectory[frame];

  //Check & statistic
  if(!isBCC(xyz,eps, rest,normal,d,err,minimal,maximal,j0, j1, j2, j3, j4, j5, origin,center))
    report << error << "Nope, could not find BCC structure!"<<endr;

  report << plain << "m="<<(xyz.coords[j0]+xyz.coords[j1]-origin.v).norm()/(2*minimal)<<","
	 <<(xyz.coords[j2]+xyz.coords[j3]-origin.v).norm()/(2*minimal)<<","
	 <<(xyz.coords[j4]+xyz.coords[j5]-origin.v).norm()/(2*minimal)<<"\n"
	 << "d="<<(xyz.coords[j0]-origin.v).norm()/minimal<<","
	 << (xyz.coords[j1]-origin.v).norm()/minimal<<","
	 << (xyz.coords[j2]-origin.v).norm()/minimal<<","
	 << (xyz.coords[j3]-origin.v).norm()/minimal<<","
	 << (xyz.coords[j4]-origin.v).norm()/minimal<<","
	 << (xyz.coords[j5]-origin.v).norm()/minimal<<"\n"
	 << "i ="<<j0<<","
	 << j1<<","
	 << j2<<","
	 << j3<<","
	 << j4<<","
	 << j5<<endr;
  report << plain <<"-rot "<<toString(normal.x)<<" "<<toString(normal.y)<<" "<<toString(normal.z)<<endr;
  
  Vector3DBlock n111s;
  n111s.push_back(xyz.coords[j0]-origin.v);
  n111s.push_back(xyz.coords[j1]-origin.v);
  n111s.push_back(xyz.coords[j2]-origin.v);
  n111s.push_back(xyz.coords[j3]-origin.v);
  n111s.push_back(xyz.coords[j4]-origin.v);
  n111s.push_back(xyz.coords[j5]-origin.v);
  sort(n111s.begin(),n111s.end(),radialCmpVector3D());
  for(unsigned int i=0;i<n111s.size();i++){
    report << plain << "-rot "<<toString(n111s[i].x)<<" "<<toString(n111s[i].y)<<" "<<toString(n111s[i].z)<<endr;
  }  
#if 0
  XYZ tmp; 
  tmp.coords.push_back(xyz.coords[origin.i]);
  tmp.coords.push_back(xyz.coords[j0]);
  tmp.coords.push_back(xyz.coords[j1]);
  tmp.coords.push_back(xyz.coords[j2]);
  tmp.coords.push_back(xyz.coords[j3]);
  tmp.coords.push_back(xyz.coords[j4]);
  tmp.coords.push_back(xyz.coords[j5]);
  tmp.names.push_back(xyz.names[origin.i]);
  tmp.names.push_back(xyz.names[j0]);
  tmp.names.push_back(xyz.names[j1]);
  tmp.names.push_back(xyz.names[j2]);
  tmp.names.push_back(xyz.names[j3]);
  tmp.names.push_back(xyz.names[j4]);
  tmp.names.push_back(xyz.names[j5]);
  XYZWriter writer("plane.xyz");
  writer << tmp;

  Matrix3by3 rot;
  rot.rotate(normal,Vector3D(0.0,0.0,1.0));
  for(unsigned int i=0;i<tmp.coords.size();i++){
    tmp.coords[i] += normal*d-origin.v;
    tmp.coords[i] = rot*tmp.coords[i];
  }
  writer.open("plane.rot.xyz");
  writer << tmp;

  writer.open("out.xyz");
  writer << xyz;
#endif

  if(!out.empty()){
    // Rotate
    Matrix3by3 rot;
    rot.rotate(normal,Vector3D(0.0,0.0,1.0));
     
    for(unsigned int j=0;j<trajectory.size();j++){
      for(unsigned int i=0;i<trajectory[j].coords.size();i++){
	trajectory[j].coords[i] += normal*d-origin.v;
	trajectory[j].coords[i] = rot*trajectory[j].coords[i];
      }
    }
     
    // Write back
    DCDTrajectoryWriter writer(out);
    if(!writer)
      report << error << "Could not open '" << out << "'." << endr;
    for(unsigned int i=0;i<trajectory.size();++i)
      writer << trajectory[i];
  }

  return 0;
}
