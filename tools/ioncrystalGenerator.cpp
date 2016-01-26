#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
using std::min;
using std::max;


#ifdef WIN32
#define drand48() (double(rand())/double(RAND_MAX))
#define srand48(seed) srand(seed)

#define	M_PI	3.14159265358979323846
#define M_PI_2  1.57079632679489661923	// pi/2
#define	M_2PI	6.283185307179586476926 
#endif

#ifdef SVR4
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/param.h>
#elif WIN32
#include <time.h>
#include <windows.h>  // mmsystem.h cannot exist without this.
#include <mmsystem.h> // timeGetTime()
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

/* 
   CC -LANG:std -O2 -LANG:restrict=ON -LANG:ansi-for-init-scope=ON ioncrystalGenerator.cpp -o ioncrystalGenerator -lm
   g++ -O9 -ffast-math -finline-functions -funroll-loops ioncrystalGenerator.cpp -o ioncrystalGenerator -lm
   xlC -qnolm -q64 -b64 -O3 -qmaxmem=-1 -qalign=natural -qcache=auto -qansialias -qarch=auto -qtune=auto -qrtti=dyna ioncrystalGenerator.cpp -o ioncrystalGenerator
*/


//_____________________________________________________________________ 
const std::string copyright="Ioncrystal Generator, 2002, matthey@ii.uib.no, "+std::string(__DATE__)+".";
const double EPSILON = 1e-6;
static const double SQRTCOULOMBCONSTANT = 18.2226123264;
//_____________________________________________________________________ Atom
struct Atom {
  Atom():x(0.0),y(0.0),z(0.0),i(0),tag(0),type(0){};
  Atom(double a, double b, double c, long d, short e, short f):x(a),y(b),z(c),i(d),tag(e),type(f){};
  double x,y,z;
  long i;
  short tag,type;
};

//_____________________________________________________________________ Lattice
enum Lattice {
  SC  = 0,
  FCC = 1,
  BCC = 2,
  HCP = 3
};

//_____________________________________________________________________ AtomType
struct AtomType {
  long         number;
  std::string  name;
  float        mass;
  float        charge;
  AtomType():number(0),name(""),mass(0.0),charge(0.0){}
  AtomType(long n, std::string s, float m, float c):number(n),name(s),mass(m),charge(c){}
};

//_____________________________________________________________________ 
static void parseCmdLine(int argc, char* argv[],
			 std::string& outputname,
			 long& seed,
			 bool& dim3d,
			 double& temperature,
			 double& finalTemperature,
			 double& thermal,
			 double& distance,
			 double& omegar,
			 double& omegaz,
			 double& perturbation,
			 Lattice& lattice,
			 bool& cubicShape,
			 bool& strictAssign,
			 std::vector<AtomType>& atomTypes,
			 bool& multigrid,
			 long& levels,
                         long& order,
			 double& cutoff,
			 double& gridsize,
			 double& cellsize,
			 long& numsteps,
			 double& timestep,
			 long& outputfreq,
                         double& maxRadius);
static void usage(char* name);
static bool equalNocase(const std::string& s1, const std::string& s2);
static std::string formated(long i, int n);
static std::string formated(std::string s, int n);
static std::string formated(float x, int n);
static double globalTime();
static bool equal(double a, double b);
static double generate_driver (double sdv);
static bool cmpCubic (const Atom& a1, const Atom& a2);
static bool cmpSpherical (const Atom& a1, const Atom& a2);
static bool cmpType (const Atom& a1, const Atom& a2);
static bool cmpTag (const Atom& a1, const Atom& a2);
static bool cmpIndex (const Atom& a1, const Atom& a2);

//_____________________________________________________________________ 
int main (int argc, char* argv[]){

  std::cerr.precision(13);
  std::cerr <<copyright<<std::endl<<std::endl;

  std::string outputname = "";          // Outputname
  double time = globalTime()*10000;     // Seed
  long seed = (long)((time - floor(time))*100000);
  bool dim3d = true;                   // Dimension
  double temperature = 1e-3;
  double finalTemperature = 1e-5;
  double thermal = 1e-17;
  double distance = 0.0;
  double omegar = 0.0;
  double omegaz = 0.0;
  bool cubicShape = false;
  double perturbation = 0.0;
  std::vector<AtomType> atomTypes;
  double cutoff = 500000;
  double gridsize = 250000;
  double cellsize = 500000;
  long order = 6;
  long numsteps = 10000;
  double timestep = 100000000;
  long outputfreq = 1;
  bool multigrid = false;
  long levels = 3;
  Lattice lattice = SC;
  bool strictAssign = false;
  double maxRadius = 0;


  std::string cmdLine = "";
  for(int j=0;j<argc;j++)
    cmdLine += std::string(argv[j])+(j<argc-1?" ":"");

  //
  // Parse input
  //
  parseCmdLine(argc,argv,
	       outputname,
	       seed,
	       dim3d,
	       temperature,
	       finalTemperature,
	       thermal,
	       distance,
	       omegar,
	       omegaz,
	       perturbation,
	       lattice,
	       cubicShape,
	       strictAssign,
	       atomTypes,
	       multigrid,
	       levels,
               order,
	       cutoff,
	       gridsize,
	       cellsize,
	       numsteps,
	       timestep,
	       outputfreq,
               maxRadius);

  // Omega and distance
  if(maxRadius > 0.0)
    distance = 1.0;
  if(omegar == 0.0 ){
    std::cerr << "Omega must be defined."<<std::endl;
    usage(argv[0]);
    exit(0);
  }

  // Seed
  srand48(seed);


  // Atom definition
  if(atomTypes.size() < 1){
      std::cerr << "Atom definition missing."<<std::endl;
      usage(argv[0]);
      exit(0);
  }
  long total = 0;
  int atomTypeSize = atomTypes.size();
  for(int i=0;i<atomTypeSize;i++){
    if(atomTypes[i].number < 0){
      std::cerr << "Atom definition "<<i+1<<": number >= 0."<<std::endl;
      usage(argv[0]);
      exit(0);
    }
    total += atomTypes[i].number;
  }
  if(atomTypeSize==1)
    strictAssign = false;




  double paulQ = 0.0;
  double paulM = 0.0;
  for(unsigned int i=0;i<atomTypeSize;i++){
    if(total>1)
      paulQ += atomTypes[i].charge*atomTypes[i].charge*atomTypes[i].number*(atomTypes[i].number-1.0)/2.0;
    else
      paulQ += atomTypes[i].charge*atomTypes[i].charge*atomTypes[i].number;
    paulM += atomTypes[i].mass*atomTypes[i].number;
    for(unsigned int j=i+1;j<atomTypeSize;j++){
      paulQ += atomTypes[i].charge*atomTypes[j].charge*atomTypes[i].number*atomTypes[j].number;
    }
  }
  if(total>1)
    paulQ = SQRTCOULOMBCONSTANT*sqrt(paulQ/((double)total*((double)total-1.0)/2.0));
  else
    paulQ = SQRTCOULOMBCONSTANT*sqrt(paulQ);

  std::vector<double> omegaR(atomTypeSize,0.0);
  std::vector<double> omegaZ(atomTypeSize,0.0);
  for(unsigned int i=0;i<atomTypeSize;i++){
    double f =atomTypes[i].charge/atomTypes[0].charge*atomTypes[0].mass/atomTypes[i].mass;
    omegaR[i] = sqrt((omegar*omegar + 0.5*omegaz*omegaz)*f*f - 0.5*omegaz*omegaz*f);
    omegaZ[i] = sqrt(omegaz*omegaz*f);
  }


  double omega = 0.0;
  for(unsigned int i=0;i<atomTypeSize;i++){
    omega += static_cast<double>(atomTypes[i].number)*(2.0*omegaR[i]*omegaR[i] + omegaZ[i]*omegaZ[i])/3.0;
  }
  omega = sqrt(omega/static_cast<double>(total));



  double paulK    = paulM/(double)total*omega*omega*1.0e7/4184.0;
  double paulA    = pow(paulQ*paulQ/paulK,1.0/3.0);
  double paulR    = pow((double)total*paulQ*paulQ/paulK,1.0/3.0);
  double paulUHom = 9.0/10.0*pow((double)total,5.0/3.0)*paulQ*paulQ/paulA;
  double paulF    = 1.0/((double)total*paulQ*paulQ/paulA);

  double packFactor = 1;
  switch(lattice){
  case SC:
    break;
  case FCC:
    packFactor = (dim3d?4.0:2.0);
    break;
  case BCC:
    packFactor = (dim3d?2.0:1.0);
    break;
  case HCP:
    packFactor = (dim3d?sqrt(2.0):sqrt(3.0)/4.0*2.0);
    break;
  default:
    break;
  }	 

  double paulD    = pow(packFactor*(dim3d?4.0*M_PI/3.0:M_PI),1.0/3.0)*paulA;

  if(distance <= 0.0){
    distance = paulD;
  }



  double distanceX = distance;
  double distanceY = distance;
  double distanceZ = distance;
  double kx = 0.0;
  double ky = 0.0;
  double kz = 0.0;

  const double hcpDist  = 1.0/sqrt(3.0);

  switch(lattice){
  case SC:
    break;
  case FCC:
    break;
  case BCC:
    kx = 1.0/distanceX;
    ky = 1.0/distanceY;
    break;
  case HCP:
    distanceX = distance;
    distanceY = distance*sqrt(3.0);
    distanceZ = distance*2.0*sqrt(2.0/3.0);
    kz = 2.0/distanceZ;
    break;
  default:
    break;
  }	 

  long nRadius;           
  if(dim3d){
    nRadius = (long)(ceil((1+perturbation/2)*pow((float)(total/3.0/packFactor+1),(float)(1.0/3.0))));
  }
  else {
    nRadius = (long)(ceil((1+perturbation/2)*pow((float)(total/packFactor+1),(float)(1.0/2.0))));
  }
  std::cerr << "Creating the following:\n";
  for(int i=0;i<atomTypeSize;i++){
    if(atomTypes[i].number < 0){
      std::cerr << "Atom definition "<<i+1<<": number >= 0."<<std::endl;
      usage(argv[0]);
      exit(0);
    }
    if(atomTypes[i].name == ""){
      std::cerr << "Atom definition "<<i+1<<": name."<<std::endl;
      usage(argv[0]);
      exit(0);
    }
    if(atomTypes[i].mass <= 0.0){
      std::cerr << "Atom definition "<<i+1<<": mass > 0."<<std::endl;
      usage(argv[0]);
      exit(0);
    }
    std::cerr <<"Atom definition    : "<<i+1<<std::endl
	      <<"   Number of atom  : "<<atomTypes[i].number<<std::endl
	      <<"   Name            : \'"<<atomTypes[i].name<<"\'"<<std::endl
	      <<"   Mass            : "<<atomTypes[i].mass<<std::endl
	      <<"   Charge          : "<<atomTypes[i].charge<<std::endl
	      <<std::endl;
  }

  std::cerr <<"Outputname          : \'"<<outputname<<"\'"<<std::endl;
  std::cerr <<"Seed                : "<<seed<<std::endl;
  std::cerr <<"Configuration       : "<<(dim3d?"3D":"2D")<<std::endl;
  std::cerr <<"Temperature         : "<<temperature<<std::endl;
  std::cerr <<"Final temperature   : "<<finalTemperature<<std::endl;
  std::cerr <<"Thermal factor      : "<<thermal<<std::endl;
  if(maxRadius > 0.0)
    std::cerr <<"Maximal radius      : "<<maxRadius<<std::endl;
  else
    std::cerr <<"Distance            : "<<distance<<std::endl;
  std::cerr <<"k                   : "<<paulK<<std::endl;
  std::cerr <<"a                   : "<<paulA<<std::endl;
  std::cerr <<"R                   : "<<paulR<<std::endl;
  std::cerr <<"d                   : "<<paulD<<std::endl;
  std::cerr <<"m                   : "<<paulM/total<<std::endl;
  std::cerr <<"q                   : "<<paulQ/SQRTCOULOMBCONSTANT<<std::endl;
  std::cerr <<"n                   : "<<packFactor/pow(distance,3.0)<<std::endl;
  std::cerr <<"n_a                 : "<<1.0/((dim3d?4.0*M_PI/3.0:M_PI)*paulA*paulA*paulA)<<std::endl;
  std::cerr <<"Uhomo               : "<<paulUHom<<std::endl;
  std::cerr <<"Packing             : "<<packFactor<<std::endl;
  std::cerr <<"Omega               : "<<omega<<std::endl;
  for(unsigned int i=0;i<atomTypeSize;i++){
    std::cerr << "                    : w_"<<i<<"  = "
	      << sqrt((2.0*omegaR[i]*omegaR[i] + omegaZ[i]*omegaZ[i])/3.0)
	      <<"[fs-1], w_r = "<<omegaR[i]<<"[fs-1], "
	      << "w_z = "<<omegaZ[i]<<"[fs-1]"<<std::endl;
  }
  std::cerr <<"Perturbation        : "<<perturbation<<std::endl;
  std::cerr <<"Slab shape          : "<<(cubicShape?"Cubic":"Spherical")<<std::endl;
  std::cerr <<"Lattice             : ";
  switch(lattice){
  case SC:
    std::cerr <<"SC - Simple Cubic";
    break;
  case FCC:
    std::cerr <<"FCC - Face Centered Cubic";
    break;
  case BCC:
    std::cerr <<"BCC - Body Centered Cubic";
    break;
  case HCP:
    std::cerr <<"HCP - Hexagonal Closest Packing";
    break;
  default:
    std::cerr <<"Unknown!";
    break;
  }	 
  std::cerr<<std::endl <<"Assign type         : "<<(atomTypeSize>1?(strictAssign?"Strict":"Random"):"-")<<std::endl;
  if(multigrid){
    std::cerr <<"Algorithm           : Multigrid"<<std::endl; 
    std::cerr <<"Number of MG levels : "<<levels<<std::endl; 
    std::cerr <<"Interpolation order : "<<order<<std::endl;  
    std::cerr <<"Cutoff              : "<<cutoff<<std::endl;
    std::cerr <<"Gridsize            : "<<gridsize<<std::endl;
    std::cerr <<"Cellsize            : "<<cellsize<<std::endl;
  }
  else {
    std::cerr <<"Algorithm           : Direct method"<<std::endl;     
  }
  std::cerr <<"Number of steps     : "<<numsteps<<std::endl;
  std::cerr <<"Timestep            : "<<timestep<<std::endl;
  std::cerr <<"Output frequency    : "<<outputfreq<<std::endl;
  std::cerr <<"Total of atoms      : "<<total<<std::endl;
  std::cerr <<"Max. spherical dim. : "<<nRadius*2+1<<"x"<<nRadius*2+1<<"x"<<(dim3d?(nRadius*2+1):1)<<" = "<<(nRadius*2+1)*(nRadius*2+1)*(dim3d?(nRadius*2+1):1)<<std::endl;
  std::cerr <<"Integrator          : STS: Nose-NVT Leapfrog"<<std::endl; 
  std::cerr <<std::endl;


  // Create positions
  std::cerr <<"Creating positions ";
  std::vector<Atom> xyzPos;
  xyzPos.clear();
  long n = 0;
  for(long i=-nRadius;i<=nRadius;i++){
    for(long j=-nRadius;j<=nRadius;j++){
      for(long k=(dim3d?-nRadius:0);k<=(dim3d?nRadius:0);k++){
	double x = distanceX*i;
	double y = distanceY*j;
	double z = (dim3d?distanceZ*k:0);	    
	xyzPos.push_back(Atom(x,y,z,n,0,0));
	n++; if(n%100000 == 0) std::cerr <<".";
	switch(lattice){
	case FCC:
	  xyzPos.push_back(Atom(x+distanceX*0.5,y+distanceY*0.5,z,n,1,0));
	  n++; if(n%100000 == 0) std::cerr <<".";
	  if(!dim3d)
	    break;	  
	  xyzPos.push_back(Atom(x+distanceX*0.5,y,z+distanceZ*0.5,n,2,0));
	  n++; if(n%100000 == 0) std::cerr <<".";
	  xyzPos.push_back(Atom(x,y+distanceY*0.5,z+distanceZ*0.5,n,3,0));
	  n++; if(n%100000 == 0) std::cerr <<".";
	  break;
	case BCC:
	  if(!dim3d)
	    break;
	  xyzPos.push_back(Atom(x+distanceX*0.5,y+distanceY*0.5,z+distanceZ*0.5,n,1,0));
	  n++; if(n%100000 == 0) std::cerr <<".";	  
	  break;
	case HCP:
	  if(!dim3d){
	    xyzPos.push_back(Atom(x+distanceX*0.5,y+distanceY*0.5,z,n,1,0));
	    n++; if(n%100000 == 0) std::cerr <<".";	  
	    break;
	  }
	  xyzPos.push_back(Atom(x,y+distance*hcpDist,z+distanceZ*0.5,n,1,0));
	  n++; if(n%100000 == 0) std::cerr <<".";	  
	  xyzPos.push_back(Atom(x+distanceX*0.5,y+distanceY*0.5,z,n,2,0));
	  n++; if(n%100000 == 0) std::cerr <<".";	  
	  xyzPos.push_back(Atom(x+distanceX*0.5,y+distanceY*0.5+distance*hcpDist,z+distanceZ*0.5,n,3,0));
	  n++; if(n%100000 == 0) std::cerr <<".";	  	  
	  break;
	default:
	  break;
	}	 
      }
    }
  }
  std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;

  
  // Add perturbation  
  if(perturbation != 0.0){
    std::cerr <<"Add perturbation ";
    n = 0;
    for(long i=0;i<xyzPos.size();i++){
      xyzPos[i].x += distance*(drand48()*perturbation-perturbation*0.5);
      xyzPos[i].y += distance*(drand48()*perturbation-perturbation*0.5);
      xyzPos[i].z += (dim3d?distance*(drand48()*perturbation-perturbation*0.5):0);
      n++; if(n%100000 == 0) std::cerr <<".";
    }
    std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;
  }
  
  
  if(strictAssign && atomTypeSize > 1 && total > atomTypeSize){
    std::cerr <<"Strict assignment linear."<<std::endl;
    for(long i=0;i<xyzPos.size();i++)
      xyzPos[i].type = i%atomTypeSize;

    // Sort according the wanted structure and ...
    std::cerr <<"Sorting after initial structure ";
    if(cubicShape){
      std::sort(xyzPos.begin(),xyzPos.end(),cmpCubic);
    }
    else {
      std::sort(xyzPos.begin(),xyzPos.end(),cmpSpherical);
    }
    std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;
    // ... remove the additional atoms
    std::cerr <<"Cave out initial structure and assign atom type(s) ";
    std::stable_sort(xyzPos.begin(),xyzPos.end(),cmpType);
    std::cerr <<".";
    std::vector<Atom> tmp(total);
    long m = 0;
    long l = 0;
    n = 0;
    for(long i=0;i<xyzPos.size();i++){
      if (l >= atomTypes[m].number){
	l = 0;
	m++;
      }
      //std::cerr <<xyzPos[i].x<<" "<<xyzPos[i].y<<" "<<xyzPos[i].z<<" "<<xyzPos[i].type<<std::endl;
      if(xyzPos[i].type == m){
	tmp[n] = xyzPos[i];
	n++; if(n%100000 == 0) std::cerr <<".";
	l++;
      }
    }
    if(tmp.size() == total){
      std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;
      xyzPos=tmp;
    }
    else {
      strictAssign = false;
      std::cerr <<std::endl<<" ... could not find enough atoms for a given type, random assignment."<<std::endl;
    }
  }

  if(!strictAssign){
    // Sort according the wanted structure and ...
    std::cerr <<"Random assignment."<<std::endl;
    std::cerr <<"Cave out initial structure ";
    if(cubicShape){
      std::sort(xyzPos.begin(),xyzPos.end(),cmpCubic);
    }
    else {
      std::sort(xyzPos.begin(),xyzPos.end(),cmpSpherical);
    }
    // ... remove atoms
    std::cerr <<".";
    xyzPos.resize(total);
    std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;

    // Assign atom type
    if(atomTypeSize > 1){
      std::cerr <<"Assign atom type(s) by random ";
      std::random_shuffle(xyzPos.begin(),xyzPos.end());
      std::cerr <<".";
      std::random_shuffle(xyzPos.begin(),xyzPos.end());
      std::cerr <<".";
      std::random_shuffle(xyzPos.begin(),xyzPos.end());
      std::cerr <<".";
      std::random_shuffle(xyzPos.begin(),xyzPos.end());
      std::cerr <<".";
      std::random_shuffle(xyzPos.begin(),xyzPos.end());
      std::cerr <<".";
      long m = 0;
      long l = 0;
      n = 0;
      for(long i=0;i<xyzPos.size();i++){
	if (l >= atomTypes[m].number){
	  l = 0;
	  m++;
	}
	xyzPos[i].type = m;
	n++; if(n%100000 == 0) std::cerr <<".";
	l++;
      }
      std::cerr <<".";
      if(cubicShape){
	std::sort(xyzPos.begin(),xyzPos.end(),cmpCubic);
      }
      else {
	std::sort(xyzPos.begin(),xyzPos.end(),cmpSpherical);
      }
      std::cerr <<".";
      std::stable_sort(xyzPos.begin(),xyzPos.end(),cmpType);
      std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;
    }
  }

  double maxR = 0;
  for(long i=0;i<total;i++)
    maxR = max(maxR,xyzPos[i].x*xyzPos[i].x+xyzPos[i].y*xyzPos[i].y+xyzPos[i].z*xyzPos[i].z);
  maxR = sqrt(maxR);
  if(maxRadius > 0.0 && maxR > 0.0){
    double scale = maxRadius/maxR;
    for(long i=0;i<total;i++){
      xyzPos[i].x *=scale;
      xyzPos[i].y *=scale;
      xyzPos[i].z *=scale;
    }
    distance *= scale;
    maxR *= scale;
  }
  //
  // VEL
  //
  if(!dim3d){
    std::cerr <<"Writing velocities ";

    std::string outputnameVEL = outputname+".vel.xyz";
    std::ofstream vel;
    vel.open(outputnameVEL.c_str(), std::ios::out); 
    vel << total << std::endl;
    vel << copyright << std::endl;
    vel << std::setprecision(14);

    double kbT = temperature*0.001987191*3.0/2.0;
    long m = 0;
    long l = 0;
    for(long i=0;i<total;i++){
      if (l >= atomTypes[m].number){
	l = 0;
	m++;
      }
      double kbToverM = sqrt(kbT/atomTypes[m].mass);
      vel << atomTypes[m].name << " " 
	  << kbToverM*generate_driver(6.0)<< " " 
	  << kbToverM*generate_driver(6.0) << " " 
	  << 0.0 <<"\n";
      if((i+1)%100000 == 0)
	std::cerr <<".";
      
    }

    vel.close();
    std::cerr << " ("<<total<<")."<<std::endl;
    std::cerr << "Wrote \'"<<outputnameVEL<<"\'."<<std::endl;
  }

  //
  // XYZ
  //
  std::cerr <<"Writing positions ";

  std::string outputnameXYZ = outputname+".pos.xyz";
  std::ofstream xyz;
  xyz.open(outputnameXYZ.c_str(), std::ios::out); 
  xyz << total << std::endl;
  xyz << copyright <<" "<<cmdLine<<std::endl;
  xyz << std::setprecision(14);

  for(long i=0;i<total;i++){
    xyz << atomTypes[xyzPos[i].type].name << " " 
	<< xyzPos[i].x << " " 
	<< xyzPos[i].y << " " 
	<< xyzPos[i].z << std::endl;      
    if((i+1)%100000 == 0)
      std::cerr <<".";
  }

  xyz.close();
  std::cerr <<" ("<<xyzPos.size()<<")."<<std::endl;
  std::cerr << "Wrote \'"<<outputnameXYZ<<"\'."<<std::endl;
  std::cerr <<"Max radius : "<<maxR<<std::endl;
  std::cerr <<"Distance   : "<<distance<<std::endl;
  

  //
  // PSF
  //
  std::cerr <<"Writing PSF ";

  std::string outputnamePSF = outputname+".psf";
  std::ofstream psf;
  psf.open(outputnamePSF.c_str(), std::ios::out); 
  psf << "PSF"<<std::endl;
  psf <<std::endl;
  psf <<"       1 !NTITLE"<<std::endl;
  psf <<" REMARKS "<<copyright<<std::endl;
  psf <<std::endl;
  psf << formated(total,8)<<" !NATOM"<<std::endl;
  long m = 0;
  long l = 0;
  std::string s5 = " XXXX    1 IONE X    "+ formated(atomTypes[m].name,6)+ " "+ formated(atomTypes[m].charge,8) 
    + "       "+formated(atomTypes[m].mass,8)+ "          0\n";
  for(long i=1;i<=total;i++){
    if (l >= atomTypes[m].number){
      l = 0;
      m++;
      s5 = " XXXX    1 IONE X    "+ formated(atomTypes[m].name,6)+ " "+ formated(atomTypes[m].charge,8) 
	+ "       "+formated(atomTypes[m].mass,8)+ "          0\n";
    }//formated(i,8)
    psf << std::setw(8) << (long)(i%100000000) << s5;
    l++;
    if((i+1)%100000 == 0)
      std::cerr <<".";

  }

  psf.close();
  std::cerr <<" ("<<total<<")."<<std::endl;  
  std::cerr << "Wrote \'"<<outputnamePSF<<"\'."<<std::endl;

  //
  // PAR
  //
  std::string outputnamePAR = outputname+".par";
  std::ofstream par;
  par.open(outputnamePAR.c_str(), std::ios::out); 
  par << " REMARKS Charmm parameter set for water and ions v22 b4"<<std::endl
      << " REMARKS FILENAME=\""<<outputname<<".par\""<<std::endl
      << ""<<std::endl
      << " SET ECHO=FALSE END"<<std::endl
      << ""<<std::endl
      << " {>>>>>>>>>> Developmental Parameter File for Proteins <<<<<<<<<<"<<std::endl
      << "  >>>>>>>>>>>>>>>>> Using All Hydrogens (ALLH) <<<<<<<<<<<<<<<<<<"<<std::endl
      << "  >>>>>>>>>>>>>>>>>>>>>>> Jan 1993 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl
      << "  >>>>>>> Direct comments to Alexander D. MacKerell Jr. <<<<<<<<<"<<std::endl
      << "  >>>>>> 410-706-7442 or bitnet: alex,tammy.harvard.edu <<<<<<<<<"<<std::endl
      << "  These files are a beta release; additional parameter development"<<std::endl
      << "  and testing may lead to alteration of the contents.}"<<std::endl
      << ""<<std::endl
      << "!-----------------------------------------------------------------------------------------------"<<std::endl
      << "!	VARIOUS BOND PARAMETERS:  Force Constant, Equilibrium Radius"<<std::endl
      << "!-----------------------------------------------------------------------------------------------"<<std::endl
      << ""<<std::endl
      << "!------------------------------------------------------------------------------------------------------"<<std::endl
      << "!	VARIOUS ANGLE PARAMETERS:  Force Constant, Equilibrium Angle, Urie-Bradley Force Const., "<<std::endl
      << "!									U.-B. equilibrium  (if any)"<<std::endl
      << "!------------------------------------------------------------------------------------------------------"<<std::endl
      << ""<<std::endl
      << "!-------------------------------------------------------------------------------------------------"<<std::endl
      << "!	VAN-DER-VAALS PARAMETERS: Energy Well Depth, Distance of Minimum(div. by 2**1/6)"<<std::endl
      << "!   for \"atoms itself\", i.e. for X-X interactions; for X-Y interactions the parameters"<<std::endl
      << "!   values are taken as arithmetic mean of appropriate \"atomic\" parameters."<<std::endl
      << "!   \"(1:4)\" values refer to special intramolecular interactions between atoms connected through "<<std::endl
      << "!   three bonds. These interactions are taken into consideration only when NBXMOD switch is set"<<std::endl
      << "!   to 5 in NBONDS...END statement, which is taken for default (see 10 lines above).	"<<std::endl
      << "!-------------------------------------------------------------------------------------------------"<<std::endl
      << ""<<std::endl;
  for(int i=0;i<atomTypeSize;i++)
    par << " NONBONDED  "<<formated(atomTypes[i].name,4)<<"     .0       .0       .0       .0    !"<<std::endl;
  par<< ""<<std::endl
     << ""<<std::endl
     << ""<<std::endl
     << "!---------------------------------------------------------------------------------------------"<<std::endl
     << "!---------------------------------------------------------------------------------------------"<<std::endl
     << " "<<std::endl
     << " SET ECHO=TRUE END"<<std::endl
     << ""<<std::endl;
  par.close();
  std::cerr << "Wrote \'"<<outputnamePAR<<"\'."<<std::endl;

  //
  // CONFIG
  //
  std::string s3 = (multigrid?"#":"");
  std::string s4 = (multigrid?"":"#");
  std::string outputnameCONFIG = outputname+".config";
  std::ofstream config;
  config.open(outputnameCONFIG.c_str(), std::ios::out); 
  config.precision(13);
  config <<"#"<<copyright<<std::endl
	 <<"#Autogenerated ProtoMol config file"<<std::endl
	 <<"#"<<cmdLine<<std::endl;
  config << "numsteps "<<numsteps<<std::endl
	 << "firststep 0"<<std::endl
	 << ""<<std::endl
	 << "boundaryConditions Normal"<<std::endl
	 << "cellManager Cubic"<<std::endl
	 << "cellsize "<<cellsize<<std::endl
	 << "exclude none"<<std::endl
	 << ""<<std::endl;
  if(dim3d)
    config << "temperature "<<temperature<<std::endl;
  else
    config << "#temperature "<<temperature<<std::endl;
  config << "posfile "<<outputname<<".pos.xyz"<<std::endl;
  if(dim3d)
    config << "#velfile "<<outputname<<".vel.xyz "<<std::endl;
  else
    config << "velfile "<<outputname<<".vel.xyz "<<std::endl;
  config << ""<<std::endl
	 << "psffile "<<outputname<<".psf"<<std::endl
	 << "parfile "<<outputname<<".par"<<std::endl
	 << ""<<std::endl
	 << "outputfreq "<<outputfreq<<std::endl
	 << "FINXYZPOSFILE     "<<outputname<<".out.fin.pos.xyz"<<std::endl
	 << "#FINXYZVELFILE     "<<outputname<<".out.fin.vel.xyz"<<std::endl
	 << "PAULOUTPUTFREQ 1"<<std::endl
	 << "PAULHISTSIZE 0"<<std::endl
	 << "PAULFILE          "<<outputname<<".out.paul"<<std::endl
	 << "PAULLOWFILE       "<<outputname<<".out.paul.low.xyz"<<std::endl
	 << "PAULOMEGA         "<<omegar<<std::endl
	 << "PAULOMEGAZ        "<<omegaz<<std::endl
	 << "#XYZPOSFILE "<<outputname<<".out.trajectory.pos.xyz"<<std::endl
	 << "#XYZPOSOUTPUTFREQ "<<(long)(1e3)<<std::endl
	 << "#ALLENERGIESFILE   "<<outputname<<".out.energies"<<std::endl
	 << ""<<std::endl
	 << "usecharmm28parfile no"<<std::endl
	 << "seed 1234"<<std::endl
	 << ""<<std::endl
	 <<     "Integrator {"<<std::endl
	 <<     "  level 0 NoseNVTLeapfrog {"<<std::endl
	 <<     "      timestep "<<timestep<<std::endl
	 <<     "      temperature "<<finalTemperature<<std::endl
	 <<     "      thermal "<<thermal<<std::endl
         <<s4<< "    force Coulomb -algorithm MultiGrid -interpolation Hermite -kernel C2 "<<std::endl
         <<s4<< "          -s "<<cutoff<<" -levels "<<levels<<" -h "<<gridsize<<" "<<gridsize<<" "<<gridsize<<"  -order "<<order<<" -ratio 2"<<std::endl
         <<s3<< "    force Coulomb -algorithm NonbondedSimpleFull"<<std::endl;

  if(dim3d)
    config  << "    force  PaulTrap -omegaR "<<omegar<<" -omegaZ "<<omegaz<<" -u0 0.0 -k "<<kx<<" "<<ky<<" "<<kz<<std::endl;
  else
    config  << "    force  PaulTrap -omegaR "<<omegar<<" -omegaZ 0.0 -u0 0.0 -k "<<kx<<" "<<ky<<" "<<kz<<" "<<std::endl;
  config << "  }"<<std::endl
	 << "}"<<std::endl;
  config.close(); 
  std::cerr << "Wrote \'"<<outputnameCONFIG<<"\'."<<std::endl;

  return 0;
}

//_____________________________________________________________________ usage
void usage(char* name){
  std::cerr << "Usage: "<<name<<std::endl
	    <<"         outputname              : output base name"<<std::endl
	    <<"         seed <n>                : (optional) random seed "<<std::endl
	    <<"         -3d|-2d                 : (optional, default 3D) dimension "<<std::endl
	    <<"         -temp <temperature>     : (optional) temperature [K]"<<std::endl
	    <<"         -ftemp <temperature>    : (optional) final temperature [K]"<<std::endl
	    <<"         -termal <factor>        : (optional) thermal factor"<<std::endl
	    <<"         -d <distance>           : average distance between atoms; [AA]"<<std::endl
	    <<"         -w <omega>              : Paul Trap Omega; [fs-1]"<<std::endl
	    <<"                                 : distance and/or omega required"<<std::endl
	    <<"         -wz <omega>             : (optional) Paul Trap Omega-z; [fs-1]"<<std::endl
	    <<"         -maxradius <distance>   : (optional) maximal radius; [AA]"<<std::endl
            <<"         -p <perturbation>       : (optional) perturbation of the grid structure"<<std::endl
            <<"         -c|-s                   : (optional, default spherical) cubic or sphere shape"<<std::endl
            <<"         -sc                     : (optional, default) Simple Cubic"<<std::endl
            <<"         -fcc                    : (optional) Face Centered Cubic"<<std::endl
            <<"         -bcc                    : (optional) Body Centered Cubic"<<std::endl
            <<"         -hcp                    : (optional) Hexagonal Closest Packing"<<std::endl
            <<"         -strict                 : (optional, default random assignment) strict atom type assignment"<<std::endl
            <<"         -mg                     : (optional, default direct method) Multigrid algorithm"<<std::endl
            <<"         -levels <n>             : (optional) number of multigrid levels"<<std::endl
            <<"         -cutoff <distance>      : (optional) softening/cutoff distance; [AA]"<<std::endl
            <<"         -gridsize <distance>    : (optional) meshsize; [AA]"<<std::endl
            <<"         -cellsize <distance>    : (optional) cell size; [AA]"<<std::endl
            <<"         -numsteps <n>           : (optional) number of steps"<<std::endl
            <<"         -timestep <time>        : (optional) time step [fs]"<<std::endl
            <<"         -outputfreq <n>         : (optional) output frequency"<<std::endl
	    <<"         <number of atoms> <atom name> <mass> <charge> : atom definition; [amu], [e]"<<std::endl
	    <<"         [<number of atoms> <atom name> <mass> <charge> [ ...]]"<<std::endl;
}

//_____________________________________________________________________ formated
std::string formated(long i, int n){
  std::stringstream ss;
  ss << i;
  std::string tmp = std::string(n,' ')+ss.str();
  return tmp.substr(tmp.size()-n,tmp.size());
}

//_____________________________________________________________________ formated
std::string formated(float x, int n){
  std::stringstream ss;
  ss << x;
  std::string tmp = ss.str()+std::string(n,' ');
  return tmp.substr(0,n);
} 

//_____________________________________________________________________ formated
std::string formated(std::string s, int n){
  std::stringstream ss;
  ss << s;
  std::string tmp = ss.str()+std::string(n,' ');
  return tmp.substr(0,n);
} 


//_____________________________________________________________________ equalNocase
bool equalNocase(const std::string& s1, const std::string& s2){
  if(s1.size() != s2.size())
    return false;
  for(int i=0;i<s1.size();i++)
    if(::toupper(s1[i]) != ::toupper(s2[i]))
      return false;
  return true;
}

//_____________________________________________________________________ globalTime
double globalTime(){
#if WIN32
  return (double)timeGetTime() / 1000.0;
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec+tv.tv_usec*.000001;
#endif
}

//_____________________________________________________________________ equal
bool equal(double a, double b){
  if(a == b)
    return true;
  if(fabs(a-b) <= (fabs(a)+fabs(b))*EPSILON)
    return true;
  return false;
}

//_____________________________________________________________________ generate_driver
double generate_driver (double sdv){
  double rnd = 0.0;
  for (int i = 0; i < 2*sdv; i++)
    rnd += drand48();
  rnd -= sdv;
  return rnd;
}

//_________________________________________________________________ parseCmdLine
void parseCmdLine(int argc, char* argv[],
		  std::string& outputname,
		  long& seed,
		  bool& dim3d,
		  double& temperature,
		  double& finalTemperature,
		  double& thermal,
		  double& distance,
		  double& omega,
		  double& omegaz,
		  double& perturbation,
		  Lattice& lattice,
		  bool& cubicShape,
		  bool& strictAssign,
		  std::vector<AtomType>& atomTypes,
		  bool& multigrid,
		  long& levels,
                  long& order,
		  double& cutoff,
		  double& gridsize,
		  double& cellsize,
		  long& numsteps,
		  double& timestep,
		  long& outputfreq,
                  double& maxRadius){
  
  if (argc < 2 || (argc >= 2 && (std::string(argv[1]) =="-h" ||
                                 std::string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  outputname = std::string(argv[cur++]);

  while (cur<(argc-1) && argv[cur][0]=='-') {

    std::string str(argv[cur]);


    if (str == "-strict"){
      strictAssign = true;
      cur++;
      continue;
    }
    if (str == "-mg"){
      multigrid = true;
      cur++;
      continue;
    }

    if (str == "-sc"){
      lattice = SC;
      cur++;
      continue;
    }

    if (str == "-fcc"){
      lattice = FCC;
      cur++;
      continue;
    }

    if (str == "-bcc"){
      lattice = BCC;
      cur++;
      continue;
    }

    if (str == "-hcp"){
      lattice = HCP;
      cur++;
      continue;
    }

    if (str == "-seed") {
      seed = atoi(argv[++cur]);
      cur++;
      continue;
    }
    if (str == "-levels") {
      levels = atoi(argv[++cur]);
      if(levels < 1){
	std::cerr << "Levels >= 1."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }


    if (str == "-order") {
      order = atoi(argv[++cur]);
      if(order < 1){
	std::cerr << "Order >= 1."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-3d"){
      dim3d = true;
      cur++;
      continue;
    }

    if (str == "-2d"){
      dim3d = false;
      cur++;
      continue;
    }

    if (str == "-temp") {
      temperature = atof(argv[++cur]);
      if(temperature<0.0){
	std::cerr << "Temerature >= 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-ftemp") {
      finalTemperature = atof(argv[++cur]);
      if(finalTemperature<0.0){
	std::cerr << "Final temerature >= 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-thermal") {
      thermal = atof(argv[++cur]);
      if(temperature<0.0){
	std::cerr << "Thermal factor >= 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-d") {
      distance = atof(argv[++cur]);
      if(distance<=0.0){
	std::cerr << "Distance > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-maxradius") {
      maxRadius = atof(argv[++cur]);
      if(maxRadius<=0.0){
	std::cerr << "Maximal radius > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-w") {
      omega = atof(argv[++cur]);
      if(omega<=0.0){
	std::cerr << "Omega > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-wz") {
      omegaz = atof(argv[++cur]);
      if(omegaz<=0.0){
	std::cerr << "Omega > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-p") {
      perturbation = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-c"){
      cubicShape = true;
      cur++;
      continue;
    }

    if (str == "-s"){
      cubicShape = false;
      cur++;
      continue;
    }
    if (str == "-cutoff") {
      cutoff = atof(argv[++cur]);
      if(cutoff<=0.0){
	std::cerr << "Cutoff > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }
    if (str == "-gridsize") {
      gridsize = atof(argv[++cur]);
      if(gridsize<=0.0){
	std::cerr << "Gridsize > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    if (str == "-cellsize") {
      cellsize = atof(argv[++cur]);
      if(cellsize<=0.0){
	std::cerr << "Cellsize > 0.0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }
    if (str == "-numsteps") {
      numsteps = atoi(argv[++cur]);
      if(numsteps<=0){
	std::cerr << "Number of steps  > 0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }
    if (str == "-timestep") {
      timestep = atof(argv[++cur]);
      if(timestep<=0.0){
	std::cerr << "Integration time step > 0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }
    if (str == "-outputfreq") {
      outputfreq = atoi(argv[++cur]);
      if(outputfreq<=0){
	std::cerr << "Output frequency > 0."<<std::endl;
	usage(argv[0]);
	exit(0);
      }
      cur++;
      continue;
    }

    break;

  }
  if(omegaz == 0 && omega > 0.0)
    omegaz = omega;


  // Check if the atom defintion is complete
  if(((argc-cur) % 4 != 0) || argc-cur < 4){
    std::cerr << "Atom definition(s) incomplete, remaining "<<argc-cur<<" argument(s)."<<std::endl;
    usage(argv[0]);
    exit(0);
  }

  //
  // Atom definitions
  //
  for(int i=cur;i<argc;i+=4)
    atomTypes.push_back(AtomType(atoi(argv[i]),std::string(argv[i+1]),atof(argv[i+2]),atof(argv[i+3])));

}

//_____________________________________________________________________ cmpSpherical
bool cmpSpherical (const Atom& a1, const Atom& a2){
    double r1 = a1.x*a1.x+a1.y*a1.y+a1.z*a1.z;
    double r2 = a2.x*a2.x+a2.y*a2.y+a2.z*a2.z;
    if (r1 < r2) return true;
    if (r1 > r2) return false;
    if (a1.x<a2.x) return true;
    if (a1.x>a2.x) return false;
    if (a1.y<a2.y) return true;
    if (a1.y>a2.y) return false;
    if (a1.z<a2.z) return true;
    return false;
}

//_____________________________________________________________________ cmpCubic
bool cmpCubic (const Atom& a1, const Atom& a2){
    double r1 = max(max(fabs(a1.x),fabs(a1.y)), fabs(a1.z));
    double r2 = max(max(fabs(a2.x),fabs(a2.y)), fabs(a2.z));
    if (r1 < r2) return true;
    if (r1 > r2) return false;
    if (a1.x<a2.x) return true;
    if (a1.x>a2.x) return false;
    if (a1.y<a2.y) return true;
    if (a1.y>a2.y) return false;
    if (a1.z<a2.z) return true;
    return false;
}

//_____________________________________________________________________ cmpIndex
bool cmpIndex (const Atom& a1, const Atom& a2){
  return (a1.i < a2.i);
}

//_____________________________________________________________________ cmpTag
bool cmpTag (const Atom& a1, const Atom& a2){
  return (a1.tag < a2.tag);
}

//_____________________________________________________________________ cmpType
bool cmpType (const Atom& a1, const Atom& a2){
  return (a1.type < a2.type);
}
