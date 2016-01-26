/*

g++ xyzread.cpp -O2 -o xyzread
xlC xyzread.cpp -q64 -O3 -qmaxmem=-1 -qstrict -qarch=pwr4 -qtune=pwr4 -o xyzread
*/

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <algorithm>

#include <math.h>

const std::string copyright="XYZ/bin Reader/writer, 2004, matthey@ii.uib.no, "+std::string(__DATE__)+".";

//_________________________________________________________________ Real
//typedef float Real;
typedef double Real;

//_________________________________________________________________ Atom
struct Atom {
  Atom():x(0.0),y(0.0),z(0.0),r(0.0),name(""){};
  Atom(double a, double b, double c, std::string n):x(a),y(b),z(c),r(sqrt(a*a+b*b+c*c)),name(n){};
  Atom operator-(const Atom& a) const  {
    return Atom(x-a.x,y-a.y,z-a.z,a.name+"-"+name);
  }
  double x,y,z,r;
  std::string name;
};


//_____________________________________________________________________ Forward declaration
static void parseCmdLine(int argc, char* argv[], std::string& inputname, std::string& outputname, bool& rb, bool& wb);
static void usage(char* name);
static bool cmpSpherical (const Atom& a1, const Atom& a2);
static bool cmpCubic (const Atom& a1, const Atom& a2);
static void writeXYZ(std::ostream& outfile, std::vector<Atom> atoms, std::string str="");
static std::vector<Atom> readXYZ(std::istream& infile);
static std::vector<Atom> readBin(std::istream& infile);

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
    double r1 = std::max(std::max(fabs(a1.x),fabs(a1.y)), fabs(a1.z));
    double r2 = std::max(std::max(fabs(a2.x),fabs(a2.y)), fabs(a2.z));
    if (r1 < r2) return true;
    if (r1 > r2) return false;
    if (a1.x<a2.x) return true;
    if (a1.x>a2.x) return false;
    if (a1.y<a2.y) return true;
    if (a1.y>a2.y) return false;
    if (a1.z<a2.z) return true;
    return false;
}



//_____________________________________________________________________ readXYZ
std::vector<Atom> readXYZ(std::istream& infile){
  
  std::string tempstr;
  int n;
  infile >> n; // number of atoms
  
  std::getline (infile, tempstr); // Next line
  std::getline (infile, tempstr); // comment
  double x,y,z;

  std::vector<Atom> atoms;
  atoms.resize(n);
  
  for(int i = 0; i < n; i++){ // now loop through the atoms	
    infile >> tempstr;
    if(infile.eof())
      std::cerr << n << " elements (atomname x y z), read " << i << ", reached end of file.\n";
    infile >> x;
    infile >> y;
    infile >> z;
    atoms[i] = Atom(x,y,z,tempstr);
  }
  return atoms;
}
//_____________________________________________________________________ readBin
std::vector<Atom> readBin(std::istream& infile){
  infile.seekg (0, std::ios::end);
  long size = infile.tellg();
  size -= 4;
  infile.seekg (0, std::ios::beg);
  int n = 0;

  if(sizeof(int) == 4){
    infile.read(reinterpret_cast<char*>(&n),sizeof(int));
    std::cerr <<"Using 32-bit int, ";
  }
  else if(sizeof(short) == 4){
    short l = 0;
    infile.read(reinterpret_cast<char*>(&l),sizeof(short));    
    n = static_cast<int>(l);
    std::cerr <<"Using 32-bit short, ";
  }
  else if(sizeof(long) == 4){
    long l = 0;
    infile.read(reinterpret_cast<char*>(&l),sizeof(long));    
    n = static_cast<long>(l);
    std::cerr <<"Using 32-bit long, ";
  }
  else {
    std::cerr <<"No adequate type found.\n";    
    exit(-1);
  }

  std::vector<Atom> atoms(n);

  if(size/(3*sizeof(float)) == n){
    std::cerr <<sizeof(float)<<"-byte float.\n";    

    float* vec = new float[n*3];
    infile.read(reinterpret_cast<char*>(vec),n*3*sizeof(float));    

    for(unsigned int i=0;i<n;i++){
      atoms[i].x = static_cast<double>(vec[i*3+0]);
      atoms[i].y = static_cast<double>(vec[i*3+1]);
      atoms[i].z = static_cast<double>(vec[i*3+2]);
      atoms[i].name = "A";
    }


    delete [] vec;
  }
  else if(size/(3*sizeof(double)) == n){
    std::cerr <<sizeof(double)<<"-byte double.\n";    

    double* vec = new double[n*3];
    infile.read(reinterpret_cast<char*>(vec),n*3*sizeof(double));

    for(unsigned int i=0;i<n;i++){
      atoms[i].x = static_cast<double>(vec[i*3+0]);
      atoms[i].y = static_cast<double>(vec[i*3+1]);
      atoms[i].z = static_cast<double>(vec[i*3+2]);
      atoms[i].name = "A";
    }


    delete [] vec;
  }
  else {
    std::cerr <<"no adequate type found.\n";    
    exit(-1);
  }

  return atoms;
}
//_____________________________________________________________________ writeXYZ
void writeXYZ(std::ostream& outfile, std::vector<Atom> atoms, std::string str){
  
  outfile << atoms.size() <<std::endl << "xyzread";
  if(str != "")
    outfile << " ("<< str << ")";
  outfile <<std::endl;
  outfile << std::setprecision(14);

  for(int i = 0; i < atoms.size(); i++)
    outfile << atoms[i].name <<" "<< atoms[i].x <<" "<< atoms[i].y <<" "<< atoms[i].z <<std::endl;
}

//_____________________________________________________________________ readBin
void writeBin(std::ostream& outfile, std::vector<Atom> atoms, std::string str){

  if(sizeof(int) == 4){
    int n = static_cast<int>(atoms.size());
    std::cerr <<"Using 32-bit int, ";
    outfile.write(reinterpret_cast<char*>(&n), 4);
  }
  else if(sizeof(short) == 4){
    short n = static_cast<short>(atoms.size());
    outfile.write(reinterpret_cast<char*>(&n), 4);
    std::cerr <<"Using 32-bit short, ";
  }
  else if(sizeof(long) == 4){
    long n = static_cast<long>(atoms.size());
    outfile.write(reinterpret_cast<char*>(&n), 4);
    std::cerr <<"Using 32-bit long, ";
  }
  else {
    std::cerr <<"No adequate type found.\n";    
    exit(-1);
  }

  if(sizeof(Real) == sizeof(float)){
    std::cerr <<sizeof(float)<<"-byte float.\n";    

    float* vec = new float[atoms.size()*3];
    for(unsigned int i=0;i<atoms.size();i++){
      vec[i*3+0] = static_cast<float>(atoms[i].x);
      vec[i*3+1] = static_cast<float>(atoms[i].y);
      vec[i*3+2] = static_cast<float>(atoms[i].z);
    }
    outfile.write(reinterpret_cast<char*>(vec),atoms.size()*3*sizeof(float));    

    delete [] vec;
  }
  else   if(sizeof(Real) == sizeof(double)){
    std::cerr <<sizeof(double)<<"-byte double.\n";    

    double* vec = new double[atoms.size()*3];
    for(unsigned int i=0;i<atoms.size();i++){
      vec[i*3+0] = static_cast<double>(atoms[i].x);
      vec[i*3+1] = static_cast<double>(atoms[i].y);
      vec[i*3+2] = static_cast<double>(atoms[i].z);
    }
    outfile.write(reinterpret_cast<char*>(vec),atoms.size()*3*sizeof(double));    

    delete [] vec;
  }
  else {
    std::cerr <<"no adequate type found.\n";    
    exit(-1);
  }

}
//_____________________________________________________________________ parseCmdLine
void parseCmdLine(int argc, char* argv[], std::string& inputname, std::string& outputname, bool& rb, bool& wb){
  if (argc < 1 || (argc >= 2 && (std::string(argv[1]) =="-h" ||
                                 std::string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  while (cur<argc && argv[cur][0]=='-') {

    std::string str(argv[cur]);
    if (str == "-rb"){
      rb = true;
      cur++;
      continue;
    }
    if (str == "-wb"){
      wb = true;
      cur++;
      continue;
    }

    break;

  }

  inputname = "";
  outputname = "";

  if (cur < argc)
    inputname = std::string(argv[cur++]);
  if (cur < argc)
    outputname = std::string(argv[cur]);
  
}

//_____________________________________________________________________ usage
void usage(char* name){
  std::cerr <<copyright<<std::endl;
  std::cerr << "Usage: " << name << " [-h][--help][-rb][-wb] [inputname [outputname]]\n"
            << " where:\n"
            << "    -h or --help          : This help output\n"
            << "    inputname             : (optional) filename of input\n"
            << "    outputname            : (optional) filename output\n";
}

//_____________________________________________________________________ main
int main (int argc, char** argv) {

  std::cerr <<copyright<<std::endl<<std::endl;

  std::vector<Atom> atoms;
  std::string inputname;
  std::string outputname;
  bool rb = false;
  bool wb = false;


  // Command line
  parseCmdLine(argc,argv,inputname,outputname, rb, wb);
  

  // Input
  if(inputname != ""){
    std::ifstream infile;
    infile.open(inputname.c_str(), std::ios::in); // open file for reading

    if (!infile){
      std::cerr << "Cannot open \'" << inputname << "\'.\n";
      exit(0);
    }
    if(rb)
      atoms = readBin(infile);      
    else
      atoms = readXYZ(infile);
    infile.close();
  }
  else {
    if(rb)
      atoms = readBin(std::cin);      
    else
      atoms = readXYZ(std::cin);
  }
  std::cerr <<"Read "<<atoms.size()<<" position(s) from "<<(inputname != ""?"\'"+inputname+"\'":"cin")<<"."<<std::endl;


  
  // Output
  if(outputname != ""){
    std::ofstream outfile;
    outfile.open(outputname.c_str(), std::ios::out); // open file for writing
    if (!outfile){
      std::cerr << "Cannot open \'" << outputname << "\'.\n";
      exit(0);
    }
    if(wb)
      writeBin(outfile, atoms, inputname);
    else
      writeXYZ(outfile, atoms, inputname);
    outfile.close();
  }
  else {
    if(wb)
      writeBin(std::cout, atoms, inputname);
    else
      writeXYZ(std::cout, atoms, inputname);
  }
  std::cerr <<"Wrote "<<atoms.size()<<" position(s) to "<<(outputname != ""?"\'"+outputname+"\'":"cout")<<"."<<std::endl;

  // Done

  return 0;
}
