/*  -*- c++ -*-  */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <set>
#include <map>
#include <cctype>
#include <sstream>
#include <iomanip>

/*
  CC -O2 -LANG:std -LANG:restrict=ON -LANG:ansi-for-init-scope=ON ioncrystalSelection.C -o ioncrystalSelection
  g++ -O9 -ffast-math -finline-functions -funroll-loops ioncrystalSelection.C -o ioncrystalSelection -lm
  xlC -q64 -O3 -qmaxmem=-1 -qstrict -qarch=pwr4 -qtune=pwr4 ioncrystalSelection.C -o ioncrystalSelection -lm
  xlC -q64  -qmaxmem=-1 -qstrict -qarch=pwr4 -qtune=pwr4 ioncrystalSelection.C -o ioncrystalSelection
*/

//_________________________________________________________________ static const
static const std::string copyright="Coulomb Crystal Selection, 2002, matthey@ii.uib.no, "+std::string(__DATE__)+".";

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

//_________________________________________________________________ Interval
struct Interval {
  Interval():d(0.0),a(0.0),b(0.0){};
  Interval(double x, double y):d(fabs(y-x)),a(std::min(x,y)),b(std::max(x,y)){};
  bool inside(double x) const{return (a <= x && x <= b);};
  bool outside(double x) const{return (x < a || b < x);};
  double d,a,b;
  bool operator<(const Interval &i)  const{return (d<i.d);}
};

//_________________________________________________________________ Cell
class Cell {
public:
  Cell():x(0),y(0),z(0){}
  Cell(int a, int b, int c): x(a),y(b),z(c){}
  Cell(Atom a, double r): x((int)floor(a.x/r)),y((int)floor(a.y/r)),z((int)floor(a.z/r)){}
public:
  bool operator<(const Cell & c) const{
    if (x<c.x) return true;
    if (x>c.x) return false;
    if (y<c.y) return true;
    if (y>c.y) return false;
    if (z<c.z) return true;
    return false;
  }
  bool operator==(const Cell & c) const{
    if (x == c.x && y == c.y && z == c.z) return true;
    return false;
  };
public:
  int x,y,z;
};

//_________________________________________________________________ usage
static void usage(char* name);
//_________________________________________________________________ lowercase
static std::string lowercase(const std::string& word);
//_________________________________________________________________ formated
static std::string formated(int i, int n);
//_________________________________________________________________ formated
static std::string formated(std::string s, int n);
//_________________________________________________________________ formated
static std::string formated(float x, int n);
//_________________________________________________________________ parseCmdLine
static void parseCmdLine(int argc, char* argv[],
                         std::string& inputname,
                         std::string& inputnamePSF,
                         std::string& outputname,
                         double& scale,
                         std::vector<double>& intervalSplitList,
                         std::vector<double>& scaleList,
                         bool& innerCharge);

//_________________________________________________________________ main
int main (int argc, char* argv[]){
  std::cerr <<copyright<<std::endl<<std::endl;

  std::string inputname = "";          // Inputname
  std::string inputnamePSF = "";          // Inputname
  std::string outputname = "";         // Outputname
  std::vector<double> intervalSplitList;
  std::vector<double> scaleList;
  double scale = 1e-4;                 // Scaling factor of positions
  bool innerCharge = true;

  int counter = 0;

  //
  // Parse input
  //
  parseCmdLine(argc,argv,
               inputname,
               inputnamePSF,
               outputname,
               scale,
               intervalSplitList,
               scaleList,
               innerCharge);


  //
  // Open the input files
  //
  std::ifstream input;
  input.open(inputname.c_str(), std::ios::in);
  if(!input){
    std::cerr<< "Could not open \'"<<inputname<<"\'.\n";
    usage(argv[0]);
    exit(0);
  }
  std::ifstream inputPSF;
  inputPSF.open(inputnamePSF.c_str(), std::ios::in);
  if(!inputPSF){
    std::cerr<< "Could not open \'"<<inputnamePSF<<"\'.\n";
    usage(argv[0]);
    exit(0);
  }

  //
  // Read the header
  //
  bool dim2d = true;
  int n = 0;
  input >> n;
  std::string str;
  std::getline(input,str);
  std::getline(input,str);
  std::set<std::string> atomTypes;
  std::map<std::string,int> atomTypesMap;
  std::vector<Atom> atoms;
  double maxRadius = 0.0;

  //
  // Read the atoms
  //
  for(int i=0; i<n;i++){
    double x,y,z;
    input >> str >> x >> y >> z;
    Atom atom = Atom(x,y,z,str);
    maxRadius = std::max(atom.r,maxRadius);
    atoms.push_back(atom);
    atomTypes.insert(str);
    if(fabs(z) > 1e-20)
      dim2d = false;
  }
  input.close();
  if(n != atoms.size()){
    std::cerr<< "Input file corrupt, number of atoms incorrect.\n";
    exit(0);
  }

  //
  // Map for atom types
  //
  int index=0;
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++)
    atomTypesMap[(*j)]=index++;


  //
  // Read the PSF header
  //
  std::vector<double> charges(atomTypes.size());
  std::vector<double> mass(atomTypes.size());
  bool start = false;
  counter = 0;
  while(inputPSF >> str){
    if(!start){
      if(lowercase(str) == "!natom"){
        start = true;
        std::getline(inputPSF,str);
      }
      continue;
    }
    std::string s0,s1,s2;
    inputPSF >> str >> str >> str >> str >> s0 >> s1 >> s2;
    std::getline(inputPSF,str);
    
    //std::cerr << s0<<","<< s1<<","<<s2<<"\n";
    if(atomTypes.find(s0) == atomTypes.end()){
      std::cerr<< "PSF Input file corrupt, found atom \'"<<s0<<"\'.\n";
      exit(0);
    }
    charges[atomTypesMap[s0]] = atof(s1.c_str());
    mass[atomTypesMap[s0]] = atof(s2.c_str());
    counter++;
  }
  if(n != counter){
    std::cerr<< "PSF Input file corrupt, number of atoms incorrect.\n";
    exit(0);
  }

  //
  // Generate interval list
  //
  std::vector<Interval> intervalList;
  if(intervalSplitList.size()> 0){
    intervalList.push_back(Interval(0,intervalSplitList[0]/scale));
    for(int i=1;i<intervalSplitList.size();i++)
      intervalList.push_back(Interval(intervalSplitList[i-1]/scale,intervalSplitList[i]/scale));
    intervalList.push_back(Interval(intervalSplitList[intervalSplitList.size()-1]/scale,maxRadius));
  }
  else {
    intervalList.push_back(Interval(0,maxRadius));
  }
  scaleList.resize(intervalList.size()*atomTypes.size());
  for(int i=0;i<scaleList.size();i++)
    if(scaleList[i] <= 0.0)
      scaleList[i] = 1.0;
  
  //
  // Print info
  //
  std::cerr<< "Input:  \'"<<inputname<<"\' ("<<n<<")\n";
  std::cerr<< "Input PSF:  \'"<<inputnamePSF<<"\' ("<<n<<")\n";
  std::cerr<< (dim2d?std::string("2D"):std::string("3D"))<<".\n";
  std::cerr<< "Output: \'"<<outputname<<"\'\n";
  std::cerr<< "Max radius: "<<maxRadius*scale<<"\n";
  std::cerr<< "Atom types: "<<atomTypes.size()<<"\n";
  std::cerr<< "Interval list: ";
  for(int i=0;i<intervalList.size();i++)
    std::cerr<< "["<<intervalList[i].a*scale<<","<<intervalList[i].b*scale<<"]";
  std::cerr<< "\n";
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++){
    std::cerr<< "Scale list \'"<<(*j)<<"\':";
    for(int i=0;i<intervalList.size();i++)
      std::cerr << scaleList[i*atomTypes.size()+atomTypesMap[(*j)]]<<" ";
    std::cerr<< "\n";
  }
  //
  // Split
  //
  std::vector<std::vector<Atom> > splitAtoms(intervalList.size());
  std::vector<std::vector<int> > splitAtomCount(intervalList.size(),std::vector<int>(atomTypes.size()));
  for(int i=0;i<n;i++){
    for(int j=0;j<intervalList.size();j++){
      if(intervalList[j].inside(atoms[i].r)){
        splitAtoms[j].push_back(atoms[i]);
        splitAtomCount[j][atomTypesMap[atoms[i].name]]++;
      }
    }
  }

  for(int i=0;i<intervalList.size();i++){
    char digit[10];
    sprintf(digit,"%03d",i);

    //
    // XYZ
    // 
    std::string name = outputname + "." + std::string(digit)+".pos.xyz";
    std::ofstream output;
    output.open(name.c_str(), std::ios::out);
    if(!output.good()){
      std::cerr<< "Could not open \'"<<name<<"\'.\n";
      exit(0);       
    }
    output << splitAtoms[i].size()+(innerCharge&&i>0?1:0) << std::endl;
    output << copyright <<" ";
    for(int j=0;j<argc;j++)
      output << std::string(argv[j])<<" ";
    output << std::endl;
    output << std::setprecision(14);
    for(int j=0;j<splitAtoms[i].size();j++){
      double f = scaleList[i*atomTypes.size()+atomTypesMap[splitAtoms[i][j].name]];
      output << splitAtoms[i][j].name << " " 
             << splitAtoms[i][j].x*f << " " 
             << splitAtoms[i][j].y*f << " " 
             << splitAtoms[i][j].z*f << std::endl;
    }
    if(innerCharge&&i>0)
      output << "XINC 0.0 0.0 0.0\n"; 
    output.close();
    std::cerr << "Wrote \'"<<name<<"\' ("<<splitAtoms[i].size()+(innerCharge&&i>0?1:0)<<")."<<std::endl;

    name = outputname + "." + std::string(digit) + ".psf";
    std::ofstream outputPSF;
    outputPSF.open(name.c_str(), std::ios::out);
    if(!outputPSF.good()){
      std::cerr<< "Could not open \'"<<name<<"\'.\n";
      exit(0);       
    }
    outputPSF << "PSF"<<std::endl;
    outputPSF <<std::endl;
    outputPSF <<"       1 !NTITLE"<<std::endl;
    outputPSF <<" REMARKS "<<copyright<<std::endl;
    outputPSF <<std::endl;
    outputPSF << formated((int)splitAtoms[i].size()+(innerCharge&&i>0?1:0),8)<<" !NATOM"<<std::endl;
    for(int j=0;j<splitAtoms[i].size();j++){
      int index = atomTypesMap[splitAtoms[i][j].name];
      outputPSF << formated((int)j+1,8) 
                << " XXXX " 
                << "   1" 
                << " IONE " 
                << formated("X",4)<< " "
                << formated(splitAtoms[i][j].name,6)<< " "
                << formated((float)charges[index],8) << "       "
                << formated((float)mass[index],8)
                << "          0\n";
    }
    if(innerCharge&&i>0){
      double q = 0.0;
      double m = 0.0;
      for(int k=0;k<i;k++){
        for(int l=0;l<atomTypes.size();l++){
          q += splitAtomCount[k][l]*charges[l];
          m += splitAtomCount[k][l]*mass[l];
        }
      }
      outputPSF << formated((int)splitAtoms[i].size()+1,8) 
                << " XXXX " 
                << "   1" 
                << " IONE " 
                << formated("X",4)<< " "
                << formated("XINC",6)<< " "
                << formated((float)q,8) << "       "
                << formated((float)m,8)
                << "          0\n";
    }
    outputPSF.close();
    std::cerr << "Wrote \'"<<name<<"\' ("<<splitAtoms[i].size()+(innerCharge&&i>0?1:0)<<")."<<std::endl;
  }

  //
  // XYZ
  // 
  std::string name = outputname + ".scaled.pos.xyz";
  std::ofstream output;
  output.open(name.c_str(), std::ios::out);
  if(!output.good()){
    std::cerr<< "Could not open \'"<<name<<"\'.\n";
    exit(0);       
  }
  output << atoms.size() << std::endl;
  output << copyright <<" ";
  for(int i=0;i<argc;i++)
    output << std::string(argv[i])<<" ";
  output << std::endl;
  output << std::setprecision(14);
  for(int i=0;i<atoms.size();i++){
    for(int j=0;j<intervalList.size();j++){
      if(intervalList[j].inside(atoms[i].r)){  
        double f = scaleList[j*atomTypes.size()+atomTypesMap[atoms[i].name]];
        output << atoms[i].name << " " 
               << atoms[i].x*f << " " 
               << atoms[i].y*f << " " 
               << atoms[i].z*f << std::endl;
        break;
      }
    }
  }
  output.close();
  std::cerr << "Wrote \'"<<name<<"\' ("<<atoms.size()<<")."<<std::endl;

  return 0;
}

//_________________________________________________________________ lowercase
std::string lowercase(const std::string& word){
  std::string ret_str = "";

  for (unsigned int index = 0; index < word.size(); index++)
    ret_str += std::tolower(word[index]);

  return ret_str;
}

//_________________________________________________________________ parseCmdLine
 void parseCmdLine(int argc, char* argv[],
                         std::string& inputname,
                         std::string& inputnamePSF,
                         std::string& outputname,
                         double& scale,
                         std::vector<double>& intervalSplitList,
                         std::vector<double>& scaleList,
                         bool& innerCharge){

  if (argc < 2 || (argc >= 2 && (std::string(argv[1]) =="-h" ||
                                 std::string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  while (cur<(argc-1) && argv[cur][0]=='-') {

    std::string str(argv[cur]);


    if (str == "-il") {
      while(std::isdigit(argv[++cur][0])){
        intervalSplitList.push_back(atof(argv[cur]));
      }
      std::sort(intervalSplitList.begin(),intervalSplitList.end());
      if(std::string(argv[cur]) == "-li"){
        cur++;
      }
      continue;
    }

    if (str == "-nocharge"){
      innerCharge = false;
      cur++;
      continue;
    }

    if (str == "-charge"){
      innerCharge = true;
      cur++;
      continue;
    }

    if (str == "-scale") {
      scale = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-sl") {
      while(std::isdigit(argv[++cur][0])){
        scaleList.push_back(atof(argv[cur]));
      }
      if(std::string(argv[cur]) == "-ls"){
        cur++;
      }
      continue;
    }


    break;

  }

  std::string str = std::string(argv[cur++]);
  if (cur < argc-1){
    inputname = str;
    inputnamePSF = std::string(argv[cur++]);
    outputname = std::string(argv[cur]);
  }
  else if(cur < argc){
    inputname = str;
    outputname = std::string(argv[cur]);
    inputnamePSF = outputname+".psf";
  }
  else {
    outputname =  str;
    inputnamePSF = str+".psf";
    inputname  = str+".out.fin.pos.xyz";
  }
}

//_________________________________________________________________ usage
void usage(char* name){
  std::cerr <<"Usage: "<<name<<" [-h][--help] -il f0 [f1 [f2 ..]-li][-sl f0 [f1 [.]]][-ls][-scale f][-nocharge][-charge][inputname][inputnamePSF] outputname\n"
            <<" where:\n"
            <<"    -h or --help         : This help output\n"
            <<"    -il f0 [f1 [.]][-li] : interval list\n"
            <<"    -sl f0 [f1 [.]][-ls] : (optional) scaling list of positions\n"
            <<"                           (default: 1)\n"
            <<"    -nocharge            : (optional) no additional inner charge\n"
            <<"    -charge              : (optional, default) additional inner charge\n"
            <<"    -scale f             : (optional) scaling factor\n"
            <<"                           (default 1e-4)\n"
            <<"    inputname            : (optional) filename of input\n"
            <<"    inputnamePSF         : (optional) PSF filename of input\n"
            <<"                           (default: outputname.out.fin.pos.xyz)\n"
            <<"    outputname           : filename of outputs\n";
}
//_____________________________________________________________________ formated
std::string formated(int i, int n){
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

