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
#include <iomanip>
#ifdef WIN32
#define	M_PI	3.14159265358979323846
#define M_PI_2  1.57079632679489661923	// pi/2
#define	M_2PI	6.283185307179586476926 
#endif

/*
  CC -O2 -LANG:std -LANG:restrict=ON -LANG:ansi-for-init-scope=ON ioncrystalAnalyzer.cpp -o ioncrystalAnalyzer
  g++ -O9 -ffast-math -finline-functions -funroll-loops ioncrystalAnalyzer.cpp -o ioncrystalAnalyzer -lm
  xlC -q64 -O3 -qmaxmem=-1 -qstrict -qarch=pwr4 -qtune=pwr4 ioncrystalAnalyzer.cpp -o ioncrystalAnalyzer -lm
  xlC -q64  -qmaxmem=-1 -qstrict -qarch=pwr4 -qtune=pwr4 ioncrystalAnalyzer.cpp -o ioncrystalAnalyzer
*/

//_________________________________________________________________ static const
static const std::string copyright="Coulomb Crystal Distribution Analyzer, 2002, matthey@ii.uib.no, "+std::string(__DATE__)+".";
static const std::string xFigHeader1 = "#FIG 3.2\nLandscape\nCenter\nMetric\nLetter  \n100.00\nSingle\n-2\n# "+copyright+"\n";
static const std::string xFigHeader2 = "1200 2\n";
static const std::string xFigHeader = xFigHeader1+xFigHeader2;
static const int xFigXMax = 12500;
static const int xFigYMax = 9700;
static const int xFigXYMin = (xFigXMax < xFigYMax ? xFigXMax : xFigYMax);
static const int xFigColors[32] = {0,20,13,10,31,1,2,3,4,5,6,7,8,9,11,12,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30};
static const std::string xFigBW1[4] ={"1 3 0 1 7 0",
                                      "1 3 0 2 0 0",
                                      "1 3 1 2 0 0",
                                      "1 3 2 2 0 0"};
static const std::string xFigBW2[4] ={"0 20 0.000 1 0.0000",
                                      "0 -1 0.000 1 0.0000",
                                      "0 -1 4.500 1 0.0000",
                                      "0 -1 4.500 1 0.0000"};

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
//_________________________________________________________________ parseCmdLine
static void parseCmdLine(int argc, char* argv[],
                         std::string& inputname,
                         std::string& outputname,
                         int& subdivisions,
                         double& scale,
                         double& maxRadiusHistogram,
                         Interval& zplane,
                         Interval& hemisphere,
                         bool& color,
                         double& minRadiusInterval,
                         double& correlationMaxDistance,
                         Interval& correlationInterval,
                         double& cluster,
                         std::vector<double>& intervalSplitList,
                         double& xfigScaleRadius,
                         bool& correlationSpherical,
			 bool& xplane,
			 bool& yplane,
			 double& maxr,
                         double& alpha,
                         double& beta,
                         double& gamma);
//_________________________________________________________________ writeXFigTrailer
static void writeXFigTrailer(std::ofstream& xfig,
                             const std::set<std::string>& atomTypes,
                             bool color);
//_________________________________________________________________ writeXFigAtom
static void writeXFigAtom(std::ofstream& xfig,
                          int x, int y, int r, int depth, int col, bool color);

//_________________________________________________________________ main
int main (int argc, char* argv[]){
  std::cerr <<copyright<<std::endl<<std::endl;

  std::string inputname = "";          // Inputname
  std::string outputname = "";         // Outputname
  int subdivisions = 0;                // Number of subdivions for the histogram
  double scale = 1e-4;                 // Scaling factor of positions
  double maxRadiusHistogram = 0.0;     // Max radius for the histogram;
  Interval zplane;                     // Z-plane
  bool color = true;                   // XFig output
  Interval hemisphere;                 // Half interval hemisphere
  double correlationMaxDistance = 0.0;
  Interval correlationInterval;
  double minRadiusInterval = -1;
  bool correlationSpherical = false;
  double cluster = 0;
  std::vector<double> intervalSplitList;
  double xfigScaleRadius = 0;
  bool xplane = false;
  bool yplane = false;
  double maxRadius = 0.0;
  double alpha = 0.0;
  double beta = 0.0;
  double gamma = 0.0;

  int written = 0;
  std::string cmdLine = "";
  for(int j=0;j<argc;j++)
    cmdLine += std::string(argv[j])+(j<argc-1?" ":"");

  //
  // Parse input
  //
  parseCmdLine(argc,argv,
               inputname,
               outputname,
               subdivisions,
               scale,
               maxRadiusHistogram,
               zplane,
               hemisphere,
               color,
               minRadiusInterval,
               correlationMaxDistance,
               correlationInterval,
               cluster,
               intervalSplitList,
               xfigScaleRadius,
               correlationSpherical,
	       xplane,
	       yplane,
	       maxRadius,
               alpha,
               beta,
               gamma);


  //
  // Open the input file
  //
  std::ifstream input;
  input.open(inputname.c_str(), std::ios::in);
  if(!input){
    std::cerr<< "Could not open \'"<<inputname<<"\'.\n";
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
  double maxRadiusInput = 0.0;

  //
  // Read the atoms
  //
  int count = 0;
  for(int i=0; i<n;i++, count++){
    double x,y,z;
    double u,v;
    input >> str >> x >> y >> z;
    Atom atom;
    u = x*cos(alpha) - y*sin(alpha);
    v = y*cos(alpha) + x*sin(alpha);
    x = u;
    y = v;
    u = x*cos(beta) - z*sin(beta);
    v = z*cos(beta) + x*sin(beta);
    x = u;
    z = v;
    u = y*cos(gamma) - z*sin(gamma);
    v = z*cos(gamma) + y*sin(gamma);
    y = u;
    z = v;

    if(xplane)
      atom = Atom(z*scale,y*scale,x*scale,str);
    else if(yplane)
      atom = Atom(x*scale,z*scale,y*scale,str);
    else
      atom = Atom(x*scale,y*scale,z*scale,str);

    maxRadiusInput = std::max(atom.r,maxRadiusInput);
    if(maxRadius <=   0.0 || (maxRadius > 0.0 && atom.r<= maxRadius)){
      atoms.push_back(atom);
      atomTypes.insert(str);
      if(fabs(z) > 1e-20)
	dim2d = false;
    }
  }
  input.close();
  if(n != count){
    std::cerr<< "Input file corrupt, number of atoms incorrect.\n";
    exit(0);
  }

  if(maxRadius <= 0.0)
    maxRadius = maxRadiusInput;
  //
  // Map for atom types
  //
  int index=0;
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++)
    atomTypesMap[(*j)]=index++;
  if(atomTypes.size() > 4)
    color = true;

  //
  // Compute max radius for the histogram
  //
  if(maxRadius*1.01>maxRadiusHistogram)
    maxRadiusHistogram = maxRadius*1.01;

  //
  // XFig  radius scaling factor
  //
  if(xfigScaleRadius <=0.0){
    if(n > 20000){
      xfigScaleRadius = 1.0/(log(10*sqrt(n/20000.0))/log(10.0));
    }
    else{
      xfigScaleRadius = 1.0;
    }
  }
  //
  // Generate interval list
  //
  std::vector<Interval> intervalList;
  if(intervalSplitList.size()> 0){
    intervalList.push_back(Interval(0,intervalSplitList[0]));
    for(int i=1;i<intervalSplitList.size();i++)
      intervalList.push_back(Interval(intervalSplitList[i-1],intervalSplitList[i]));
    intervalList.push_back(Interval(intervalSplitList[intervalSplitList.size()-1],maxRadius));
  }

  //
  // Print info
  //
  std::cerr<< "Input:                 \'"<<inputname<<"\' ("<<n<<")\n";
  std::cerr<< (dim2d?std::string("2D"):std::string("3D"))<<".\n";
  std::cerr<< "Output:                \'"<<outputname<<"\'\n";
  std::cerr<< "Max radius:             "<<maxRadius<<"\n";
  std::cerr<< "Max radius input:       "<<maxRadiusInput<<"\n";
  std::cerr<< "Scaling factor:         "<<scale<<"\n";
  std::cerr<< "Max radius histogram:   "<<maxRadiusHistogram<<"\n";
  std::cerr<< "XFig scaling:           "<<xfigScaleRadius<<"\n";
  if(zplane.d > 0)
    std::cerr<< "Plane:                 ["<<zplane.a<<","<<zplane.b<<"]\n";
  if(hemisphere.d > 0)
    std::cerr<< "Hemisphere:            ["<<hemisphere.a<<","<<hemisphere.b<<"]\n";
  std::cerr<< (color?std::string("XFig color"):std::string("XFig B&W"))<<".\n";
  std::cerr<< "Atom types:             "<<atomTypes.size()<<"\n";
  if(correlationMaxDistance > 0.0){
    std::cerr<< "Radius correlation :    "<<correlationMaxDistance<<", [";
    if(correlationInterval.d > 0)
      std::cerr<<correlationInterval.a<<","<<correlationInterval.b;
    else
      std::cerr<<"0,"<<maxRadius;
    std::cerr<<"].\n";
  }
  if(minRadiusInterval >= 0.0)
    std::cerr<< "Minimum radius interval:  "<<minRadiusInterval<<"\n";
  if(intervalSplitList.size()> 0){
    std::cerr<< "Interval list:            ";
    for(int i=0;i<intervalList.size();i++)
      std::cerr<< "["<<intervalList[i].a<<","<<intervalList[i].b<<"]";
    std::cerr<< "\n";
  }


  //
  // Split the output for each atom type
  //
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++){

    //
    // Local copy of actual atom type
    //
    std::vector<Atom> atomsOneType;
    for(int i=0; i<n;i++)
      if((*j) == atoms[i].name)
        atomsOneType.push_back(atoms[i]);

    //
    // Write the histogram
    //
    if(subdivisions > 0){
      std::vector<int> histogram;
      histogram.resize(subdivisions);
      double a = (double)subdivisions/maxRadiusHistogram;
      for(int i=0; i<atomsOneType.size();i++)
        histogram[int(atomsOneType[i].r*a)] += 1;

      std::string outputnameHist = outputname+".hist."+(*j);
      std::ofstream hist;
      hist.open(outputnameHist.c_str(), std::ios::out);
      for(int i=0; i<histogram.size();i++)
        hist<<histogram[i]<<std::endl;
      hist.close();
      std::cerr << "Wrote \'"<<outputnameHist<<"\' ("<<histogram.size()<<").\n";
    }

    std::vector<double> radius;
    for(int i=0; i<atomsOneType.size();i++)
      radius.push_back(atomsOneType[i].r);

    //
    // Write radii
    //
    std::string outputnameRadii = outputname+".r."+(*j);
    std::ofstream radii;
    radii.open(outputnameRadii.c_str(), std::ios::out);
    radii << std::setprecision(8);
    for(int i=0; i<radius.size();i++)
      radii<<radius[i]<<std::endl;
    radii.close();
    std::cerr << "Wrote \'"<<outputnameRadii<<"\' ("<<radius.size()<<").\n";




    //
    // Write sorted radii
    //
    std::string outputnameRadiiSorted = outputname+".rs."+(*j);
    std::ofstream radiiSorted;
    radiiSorted.open(outputnameRadiiSorted.c_str(), std::ios::out);
    radiiSorted << std::setprecision(8);
    std::sort(radius.begin(),radius.end());
    for(int i=0; i<radius.size();i++)
      radiiSorted<<radius[i]<<std::endl;
    radiiSorted.close();
    std::cerr << "Wrote \'"<<outputnameRadiiSorted<<"\' ("<<radius.size()<<").\n";

    //
    // Write interval counts
    //
    if(intervalSplitList.size() > 0 && atomTypes.size()>1){
      std::vector<int> intervalListN(intervalList.size());
      std::vector<Interval> intervalReal = intervalList;
      int k=0;
      int l=0;
      while(k<radius.size()){
        if(intervalList[l].inside(radius[k])){
          intervalListN[l]++;
          if(radius[k] < intervalReal[l].a || intervalReal[l].a == intervalList[l].a)
            intervalReal[l] = Interval(radius[k],intervalReal[l].b);
          if(radius[k] > intervalReal[l].b || intervalReal[l].b == intervalList[l].b)
            intervalReal[l] = Interval(intervalReal[l].a,radius[k]);
          k++;
        }
        else {
          l++;
        }
      }

      std::string outputnameIntervalList = outputname+".il."+(*j);
      std::ofstream intervalListOut;
      intervalListOut.open(outputnameIntervalList.c_str(), std::ios::out);
      intervalListOut << std::setprecision(8);
      for(int i=0; i<intervalList.size();i++){
        intervalListOut.width(7);
        intervalListOut << intervalListN[i];
        intervalListOut.width(14);
        intervalListOut << intervalList[i].d;
        intervalListOut.width(14);
        intervalListOut << intervalList[i].a;
        intervalListOut.width(14);
        intervalListOut <<intervalList[i].b;
        intervalListOut.width(14);
        intervalListOut << intervalReal[i].d;
        intervalListOut.width(14);
        intervalListOut << intervalReal[i].a;
        intervalListOut.width(14);
        intervalListOut <<intervalReal[i].b <<"\n";
      }
      intervalListOut.close();
      std::cerr << "Wrote \'"<<outputnameIntervalList<<"\' ("<<intervalList.size()<<").\n";
    }

    //
    // Write clusters
    // Cluster together all radii <= cluster
    //
    if(cluster > 0.0 && atomTypes.size()>1){
      std::vector<Interval> clusterIntervals;
      std::vector<int> clusterIntervalsN(1);
      std::vector<double> clusterIntervalsRavg(1);
      clusterIntervalsN[clusterIntervalsN.size()-1]++;
      clusterIntervalsRavg[clusterIntervalsRavg.size()-1] +=radius[0];
      Interval interval(radius[0],radius[0]);
      for(int i=1;i<radius.size();i++){
        if(interval.b>= radius[i]-cluster){
          interval = Interval(interval.a,radius[i]);
        }
        else {
          clusterIntervals.push_back(interval);
          interval = Interval(radius[i],radius[i]);
          clusterIntervalsN.resize(clusterIntervalsN.size()+1);
          clusterIntervalsRavg.resize(clusterIntervalsRavg.size()+1);
        }
        clusterIntervalsN[clusterIntervalsN.size()-1]++;
        clusterIntervalsRavg[clusterIntervalsRavg.size()-1] +=radius[i];
      }
      clusterIntervals.push_back(interval);
      std::string outputnameClusterInterval = outputname+".cluster."+(*j);
      std::ofstream clusterInterval;
      clusterInterval.open(outputnameClusterInterval.c_str(), std::ios::out);
      clusterInterval << std::setprecision(8);
      for(int i=0;i<clusterIntervals.size();i++){
        clusterInterval.width(7);
        clusterInterval << clusterIntervalsN[i];
        clusterInterval.width(14);
        clusterInterval << clusterIntervals[i].d;
        clusterInterval.width(14);
        clusterInterval << clusterIntervals[i].a;
        clusterInterval.width(14);
        clusterInterval <<clusterIntervals[i].b;
        clusterInterval.width(14);
        clusterInterval <<clusterIntervalsRavg[i]/clusterIntervalsN[i];
        if(i>0){
          clusterInterval.width(14);
          clusterInterval << clusterIntervals[i].a-clusterIntervals[i-1].b;
        }
        clusterInterval <<"\n";
      }
      clusterInterval.close();
      std::cerr << "Wrote \'"<<outputnameClusterInterval<<"\' ("<<clusterIntervals.size()<<").\n";
    }

    //
    // Write the spherical projection of positions
    // gnuplot input
    //
    written = 0;
    std::string outputnameSphere = outputname+"."+(*j);
    std::ofstream sphere;
    sphere.open(outputnameSphere.c_str(), std::ios::out);
    sphere << std::setprecision(8);
    for(int i=0; i<atomsOneType.size();i++){
      if(dim2d){
        sphere <<atomsOneType[i].x<<" "<<atomsOneType[i].y<<std::endl;
        written++;
      }
      else{
        if(zplane.d > 0.0 && zplane.outside(atomsOneType[i].z))
          continue;
        sphere <<atomsOneType[i].x<<" "
               <<sqrt(atomsOneType[i].z*atomsOneType[i].z+atomsOneType[i].y*atomsOneType[i].y)
               <<std::endl;
        written++;
      }
    }
    sphere.close();
    std::cerr << "Wrote \'"<<outputnameSphere<<"\' ("<<written<<").\n";
  }


  //
  // Write XFig plane projection
  //
  written = 0;
  std::string outputnameXfigProjection = outputname+".fig";
  std::ofstream xfigProjection;
  xfigProjection.open(outputnameXfigProjection.c_str(), std::ios::out);
  xfigProjection << xFigHeader1<<"#"<<cmdLine<<"\n"<<xFigHeader2;
  for(int i=0; i<n;i++){
    int x,y,depth;
    if(dim2d){
      x = atoms[i].x/(2.1*maxRadius)*xFigXYMin+xFigXMax/2;
      y = xFigYMax/2-atoms[i].y/(2.1*maxRadius)*xFigXYMin;
      depth = 0;
    }
    else {
      if(zplane.d > 0.0 && zplane.outside(atoms[i].z))
        continue;

      x = atoms[i].x/(2.1*maxRadius)*xFigXYMin+xFigXMax/2;
      y = xFigYMax/2-sqrt(atoms[i].z*atoms[i].z+atoms[i].y*atoms[i].y)/(2.1*maxRadius)*xFigXYMin;
      if(zplane.d > 0.0){
        depth = 100-(atoms[i].z-zplane.a)/(zplane.d)*50.0;
      }
      else {
        depth = 100-(atoms[i].z+maxRadius)/(2.0*maxRadius)*50.0;
      }
    }
    writeXFigAtom(xfigProjection,x,y,75*xfigScaleRadius,depth,atomTypesMap[atoms[i].name],color);
    written++;
  }
  writeXFigTrailer(xfigProjection,atomTypes,color);
  xfigProjection.close();
  std::cerr << "Wrote \'"<<outputnameXfigProjection<<"\' ("<<written<<").\n";


  //
  // Write XFig hemisphere
  //
  if(hemisphere.d > 0.0){
    std::string outputnameXfig = outputname+".sphere.fig";
    std::ofstream xfig;
    xfig.open(outputnameXfig.c_str(), std::ios::out);
    xfig << xFigHeader1<<"#"<<cmdLine<<"\n"<<xFigHeader2;
    written = 0;
    for(int i=0; i<n;i++){
      if(atoms[i].z >= 0 &&  hemisphere.inside(atoms[i].r)){
        int x = atoms[i].x/(2.1*maxRadius)*xFigXYMin+xFigXMax/2;
        int y = xFigYMax/2-atoms[i].y/(2.1*maxRadius)*xFigXYMin;
        int r = atoms[i].z/maxRadius*50+50;
        int depth = 100 - (atoms[i].r-hemisphere.a)/(hemisphere.d)*50;
        writeXFigAtom(xfig,x,y,r*xfigScaleRadius,depth,atomTypesMap[atoms[i].name],color);
        written++;
      }
    }
    writeXFigTrailer(xfig,atomTypes,color);
    xfig.close();
    std::cerr << "Wrote \'"<<outputnameXfig<<"\' ("<<written<<").\n";
  }

  //
  // Write XFig z-plane
  //
  if(zplane.d > 0.0){
    std::string outputnameXfig = outputname+".plane.fig";
    std::ofstream xfig;
    xfig.open(outputnameXfig.c_str(), std::ios::out);
    xfig << xFigHeader1<<"#"<<cmdLine<<"\n"<<xFigHeader2;
    xfig << "4 0 0 100 0 0 12 0.0000 4 105 300 9090 495 "<<outputname<<"\\001\n";
    written = 0;
    for(int i=0; i<n;i++){
      if(zplane.inside(atoms[i].z)){
        int x = atoms[i].x/(2.1*maxRadius)*xFigXYMin+xFigXMax/2;
        int y = xFigYMax/2-atoms[i].y/(2.1*maxRadius)*xFigXYMin;
        int depth = 100-(atoms[i].z-zplane.a)/zplane.d*50.0;
        writeXFigAtom(xfig,x,y,75*xfigScaleRadius,depth,atomTypesMap[atoms[i].name],color);
        written++;
      }
    }
    writeXFigTrailer(xfig,atomTypes,color);
    xfig.close();
    std::cerr << "Wrote \'"<<outputnameXfig<<"\' ("<<written<<").\n";
  }

  //
  // Precompute radii and sort them
  //
  std::vector<double> radius(n);
  if(intervalSplitList.size() > 0 ||
     cluster > 0.0 ||
     minRadiusInterval >= 0.0){
    for(int i=0;i<n;i++)
      radius[i] = atoms[i].r;
    std::sort(radius.begin(),radius.end());
  }

  //
  // Write interval counts
  //
  if(intervalSplitList.size() > 0){
    std::vector<int> intervalListN(intervalList.size());
    std::vector<Interval> intervalReal = intervalList;
    int k=0;
    int l=0;
    while(k<radius.size()){
      if(intervalList[l].inside(radius[k])){
        intervalListN[l]++;
        if(radius[k] < intervalReal[l].a || intervalReal[l].a == intervalList[l].a)
          intervalReal[l] = Interval(radius[k],intervalReal[l].b);
        if(radius[k] > intervalReal[l].b || intervalReal[l].b == intervalList[l].b)
          intervalReal[l] = Interval(intervalReal[l].a,radius[k]);
        k++;
      }
      else {
        l++;
      }
    }

    std::string outputnameIntervalList = outputname+".il";
    std::ofstream intervalListOut;
    intervalListOut.open(outputnameIntervalList.c_str(), std::ios::out);
    intervalListOut << std::setprecision(8);
    for(int i=0; i<intervalList.size();i++){
      intervalListOut.width(7);
      intervalListOut << intervalListN[i];
      intervalListOut.width(14);
      intervalListOut << intervalList[i].d;
      intervalListOut.width(14);
      intervalListOut << intervalList[i].a;
      intervalListOut.width(14);
      intervalListOut <<intervalList[i].b;
      intervalListOut.width(14);
      intervalListOut << intervalReal[i].d;
      intervalListOut.width(14);
      intervalListOut << intervalReal[i].a;
      intervalListOut.width(14);
      intervalListOut <<intervalReal[i].b <<"\n";
    }
    intervalListOut.close();
    std::cerr << "Wrote \'"<<outputnameIntervalList<<"\' ("<<intervalList.size()<<").\n";
  }

  //
  // Write cluster list of all atoms
  //
  if(cluster > 0.0){
    std::vector<Interval> clusterIntervals;
    std::vector<int> clusterIntervalsN(1);
    std::vector<double> clusterIntervalsRavg(1);
    clusterIntervalsN[clusterIntervalsN.size()-1]++;
    clusterIntervalsRavg[clusterIntervalsRavg.size()-1] +=radius[0];
    Interval interval(radius[0],radius[0]);
    for(int i=1;i<radius.size();i++){
      if(interval.b>= radius[i]-cluster){
        interval = Interval(interval.a,radius[i]);
      }
      else {
        clusterIntervals.push_back(interval);
        interval = Interval(radius[i],radius[i]);
        clusterIntervalsN.resize(clusterIntervalsN.size()+1);
        clusterIntervalsRavg.resize(clusterIntervalsRavg.size()+1);
      }
      clusterIntervalsN[clusterIntervalsN.size()-1]++;
      clusterIntervalsRavg[clusterIntervalsRavg.size()-1] +=radius[i];
    }
    clusterIntervals.push_back(interval);
    std::string outputnameClusterInterval = outputname+".cluster";
    std::ofstream clusterInterval;
    clusterInterval.open(outputnameClusterInterval.c_str(), std::ios::out);
    clusterInterval << std::setprecision(8);
    for(int i=0;i<clusterIntervals.size();i++){
      clusterInterval.width(7);
      clusterInterval << clusterIntervalsN[i];
      clusterInterval.width(14);
      clusterInterval << clusterIntervals[i].d;
      clusterInterval.width(14);
      clusterInterval << clusterIntervals[i].a;
      clusterInterval.width(14);
      clusterInterval <<clusterIntervals[i].b;
      clusterInterval.width(14);
      clusterInterval <<clusterIntervalsRavg[i]/clusterIntervalsN[i];
      if(i>0){
        clusterInterval.width(14);
        clusterInterval << clusterIntervals[i].a-clusterIntervals[i-1].b;
      }
      clusterInterval  <<"\n";
    }
    clusterInterval.close();
    std::cerr << "Wrote \'"<<outputnameClusterInterval<<"\' ("<<clusterIntervals.size()<<").\n";
  }


  //
  // Write max distance of radii
  //
  if(minRadiusInterval >= 0.0){

    std::vector<Interval> radiusIntervals;
    for(int i=0;i<radius.size()-1;i++){
      Interval interval = Interval(radius[i],radius[i+1]);
      if(interval.d>= minRadiusInterval)
        radiusIntervals.push_back(interval);
    }
    std::sort(radiusIntervals.begin(),radiusIntervals.end());

    std::string outputnameRadiusInterval = outputname+".ir";
    std::ofstream rintervals;
    rintervals.open(outputnameRadiusInterval.c_str(), std::ios::out);
    rintervals << std::setprecision(8);
    for(int i=0;i<radiusIntervals.size();i++){
      rintervals.width(14);
      rintervals << radiusIntervals[i].d;
      rintervals.width(14);
      rintervals << radiusIntervals[i].a;
      rintervals.width(14);
      rintervals <<radiusIntervals[i].b <<"\n";
    }
    rintervals.close();
    std::cerr << "Wrote \'"<<outputnameRadiusInterval<<"\' ("<<radiusIntervals.size()<<").\n";
  }

  //
  // Write radius correlation
  //
  if(correlationMaxDistance > 0.0){
    std::map<Cell,std::vector<Atom> > cellList;
    for(int i=0; i<n;i++){
      if(correlationInterval.d == 0.0 || correlationInterval.inside(atoms[i].r))
        cellList[Cell(atoms[i],correlationMaxDistance)].push_back(atoms[i]);
    }
    std::vector<double> radii;
    std::map<std::string,std::vector<double> > radiiSplit;
    for(std::map<Cell,std::vector<Atom> >::const_iterator i = cellList.begin();i!=cellList.end();i++){
      for(std::map<Cell,std::vector<Atom> >::const_iterator j = i;j!=cellList.end();j++){
        if(fabs((float)(i->first.x-j->first.x))< 2 &&
           fabs((float)(i->first.y-j->first.y))< 2 &&
           fabs((float)(i->first.z-j->first.z))< 2){
          //std::cerr << "["<<i->first.x <<","<< i->first.y<<","<< i->first.z<<"] "<<i->second.size()
          //          << "["<<j->first.x <<","<< j->first.y<<","<< j->first.z<<"] "<<j->second.size()
          //          <<"\n";
          for(int k = 0;k < i->second.size();k++){
            for(int l = (i!=j?0:k+1);l < j->second.size();l++){
              Atom atom(i->second[k]-j->second[l]);
              double r = atom.r;
              if(correlationSpherical && i->second[k].r > 0.0 && j->second[l].r > 0.0){
                double xi = i->second[k].x / i->second[k].r;
                double yi = i->second[k].y / i->second[k].r;
                double zi = i->second[k].z / i->second[k].r;
                double xj = j->second[l].x / j->second[l].r;
                double yj = j->second[l].y / j->second[l].r;
                double zj = j->second[l].z / j->second[l].r;
                double a = fabs(i->second[k].r-j->second[l].r);
                double b = fabs(acos(xi*xj + yi*yj + zi*zj)*(i->second[k].r+j->second[l].r)/2.0);
                r = sqrt(a*a + b*b);                
                //std::cerr << "||r||="<<atom.r <<",r="<<r<<",a="<<a<<",b="<<b<<"alpha="<<180.0/(atan(1.0)*4.0)*acos(xi*xj + yi*yj + zi*zj)<<"\n";
              }
              if(r <= correlationMaxDistance){
                //std::cerr << r <<"\n";
                Atom tmp(j->second[l]-i->second[k]);
                atom = (atom.name < tmp.name ? atom : tmp);
                radii.push_back(r);
                radiiSplit[atom.name].push_back(r);
              }
            }
          }
        }
      }
    }
    //
    // Write radius correlation for all atoms
    //
    std::string outputnameRadiiSorted = outputname+".c";
    std::ofstream radiiSorted;
    radiiSorted.open(outputnameRadiiSorted.c_str(), std::ios::out);
    radiiSorted << std::setprecision(8);
    std::sort(radii.begin(),radii.end());
    for(int i=0; i<radii.size();i++)
      radiiSorted<<radii[i]<<std::endl;
    radiiSorted.close();
    std::cerr << "Wrote \'"<<outputnameRadiiSorted<<"\' ("<<radii.size()<<").\n";

    //
    // Write radius correlation for each interaction type
    //
    if(atomTypes.size() > 1){
      for(std::map<std::string,std::vector<double> >::iterator j = radiiSplit.begin();
          j!= radiiSplit.end();j++){

        written = 0;
        std::string outputnameRadiiSortedOne = outputname+".c."+j->first;
        std::ofstream radiiSortedOne;
        radiiSortedOne.open(outputnameRadiiSortedOne.c_str(), std::ios::out);
        radiiSortedOne << std::setprecision(8);
        std::sort(j->second.begin(),j->second.end());
        std::cerr << "Wrote \'"<<outputnameRadiiSortedOne<<"\' ("<<j->second.size()<<").\n";
        for(int i=0; i<j->second.size();i++)
          radiiSortedOne<<j->second[i]<<std::endl;
        radiiSortedOne.close();
      }
    }
    //     // Write the sorted radii, direct method
    //     std::string outputnameRadiiSortedDirect = outputname+".cd";
    //     std::ofstream radiiSortedDirect;
    //     radiiSortedDirect.open(outputnameRadiiSortedDirect.c_str(), std::ios::out);
    //     radii.clear();
    //     for(int i=0; i<atoms.size();i++){
    //       for(int j=i+1; j<atoms.size();j++){
    //         Atom atom(atoms[i]-atoms[j]);
    //         if(atom.r <= correlationMaxDistance){
    //           radii.push_back(atom.r);
    //         }
    //       }
    //     }
    //     std::sort(radii.begin(),radii.end());
    //     for(int i=0; i<radii.size();i++)
    //       radiiSortedDirect<<radii[i]<<std::endl;
    //     radiiSortedDirect.close();
    //     std::cerr << "Wrote \'"<<outputnameRadiiSortedDirect<<"\' ("<<radii.size()<<").\n";

  }

  //
  // Write gnuplot command
  //
  std::cerr << "gnuplot command:\n";
  std::cerr << "   gnuplot\n";
  std::cerr << "set terminal postscript\n";
  std::cerr << "set output \""<<outputname<<".ps\"\n";
  std::cerr << "plot";
  bool first = true;
  index=6;
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++){
    if(!first)
      std::cerr <<",";
    std::cerr << " \""<<outputname<<"."<<(*j)<<"\" title \""<<(*j)<<"\"  with points pt "<<index;
    index++;
    first = false;
  }
  std::cerr <<"\n";
  std::cerr <<"show output\n";
  std::cerr <<"quit\n";

  return 0;
}

//_________________________________________________________________ writeXFigTrailer
void writeXFigTrailer(std::ofstream& xfig, const std::set<std::string>& atomTypes, bool color){
  int i = 0;
  for(std::set<std::string>::const_iterator j = atomTypes.begin();j!= atomTypes.end();j++,i++){
    int x = 540;
    int y = 540+i*400;
    int r = 100;
    writeXFigAtom(xfig,x,y,r,0,i,color);
    int x2 = x+220;
    int y2 = y+90;
    xfig <<"4 0 0 0 0 0 18 0.0000 4 195 10000 "<<x2<<" "<<y2 <<" "<<(*j)<<"\\001\n";
  }
}

//_________________________________________________________________ writeXFigAtom
void writeXFigAtom(std::ofstream& xfig, int x, int y, int r, int depth, int col, bool color){
  if(color){
    xfig <<"1 3 0 1 "<<xFigColors[col]<<" "<<xFigColors[col]<<" "<<depth<<" 0 41 0.000 1 0.0000 "
         << x <<" "
         << y <<" "
         << r <<" "
         << r <<" "
         << x <<" "
         << y <<" "
         << x+r <<" "
         << y << std::endl;
  }
  else {
    xfig <<xFigBW1[col]<<" "<<depth<<" "<<xFigBW2[col]<<" "
         << x <<" "
         << y <<" "
         << r <<" "
         << r <<" "
         << x <<" "
         << y <<" "
         << x+r <<" "
         << y << std::endl;
  }
}

//_________________________________________________________________ parseCmdLine
void parseCmdLine(int argc, char* argv[],
                  std::string& inputname,
                  std::string& outputname,
                  int& subdivisions,
                  double& scale,
                  double& maxRadiusHistogram,
                  Interval& zplane,
                  Interval& hemisphere,
                  bool& color,
                  double& minRadiusInterval,
                  double& correlationMaxDistance,
                  Interval& correlationInterval,
                  double& cluster,
                  std::vector<double>& intervalSplitList,
                  double& xfigScaleRadius,
                  bool& correlationSpherical,
		  bool& xplane,
		  bool& yplane,
		  double& maxr,
                  double& alpha,
                  double& beta,
                  double& gamma){

  if (argc < 2 || (argc >= 2 && (std::string(argv[1]) =="-h" ||
                                 std::string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  while (cur<(argc-1) && argv[cur][0]=='-') {

    std::string str(argv[cur]);

    if (str == "-plane" && cur < (argc-2)) {
      double a = atof(argv[++cur]);
      double b = atof(argv[++cur]);
      zplane = Interval(a,b);
      cur++;
      continue;
    }

    if (str == "-sphere" && cur < (argc-2)) {
      double a = atof(argv[++cur]);
      double b = atof(argv[++cur]);
      hemisphere = Interval(a,b);
      cur++;
      continue;
    }

    if (str == "-hs") {
      subdivisions = atoi(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-hr") {
      maxRadiusHistogram = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-scale") {
      scale = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-maxr") {
      maxr = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-a") {
      alpha = atof(argv[++cur])/180.0*M_PI;
      cur++;
      continue;
    }

    if (str == "-b") {
      beta = atof(argv[++cur])/180.0*M_PI;
      cur++;
      continue;
    }

    if (str == "-c") {
      gamma = atof(argv[++cur])/180.0*M_PI;
      cur++;
      continue;
    }

    if (str == "-xfigscale") {
      xfigScaleRadius = atof(argv[++cur]);
      cur++;
      continue;
    }

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

    if (str == "-cluster") {
      cluster = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-ir") {
      minRadiusInterval = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-cr") {
      correlationMaxDistance = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-ci" && cur < (argc-2)) {
      double a = atof(argv[++cur]);
      double b = atof(argv[++cur]);
      correlationInterval = Interval(a,b);
      cur++;
      continue;
    }

    if (str == "-color"){
      color = true;
      cur++;
      continue;
    }
    if (str == "-x"){
      xplane = true;
      cur++;
      continue;
    }
    if (str == "-y"){
      yplane = true;
      cur++;
      continue;
    }

    if (str == "-bw"){
      color = false;
      cur++;
      continue;
    }

    if (str == "-scor"){
      correlationSpherical = true;
      cur++;
      continue;
    }

    if (str == "-xy"){
      alpha = 45.0/180*M_PI;
      beta  = 0.0;
      gamma = 0.0;
      cur++;
      continue;
    }

    if (str == "-xz"){
      alpha = 0.0;
      beta  = 45.0/180*M_PI;
      gamma = 0.0;
      cur++;
      continue;
    }

    if (str == "-yz"){
      alpha = 0.0;
      beta  = 0.0;
      gamma = 45.0/180*M_PI;
      cur++;
      continue;
    }

    if (str == "-xyz"){
      alpha = 45.0/180*M_PI;
      beta  = acos(1.0/sqrt(3.0));
      gamma = 0.0;
      cur++;
      continue;
    }

    break;

  }

  inputname = std::string(argv[cur++]);
  if (cur < argc){
    outputname = std::string(argv[cur]);
  }
  else {
    outputname =  inputname;
    inputname  = inputname+".out.fin.pos.xyz";
  }
}

//_________________________________________________________________ usage
void usage(char* name){
  std::cerr <<"Usage: "<<name<<" [-h][--help][-scale]][-plane f0 f1][-sphere f0 f1]"
            <<"[-max r][-hs i][-hr f][-ir f][-scale f][-xfigscale f][-cluster f][-il f0 [f1 [f2 ..]][-li]][-cr f][-ci f0 f1][-color][-bw][-scor][inputname][-x][-y][-a f][-b f][-c f] outputname\n"
            <<" where:\n"
            <<"    -h or --help         : This help output\n"
            <<"    -plane f0 f1         : (optional) XFig, sperical projection and z-plane cut [f0 f1]]\n"
            <<"    -sphere f0 f1        : (optional) XFig, front hemisphere [f0,f1]\n"
            <<"    -hs i                : (optional) histogram, number subdivisons\n"
            <<"    -hr f                : (optional) histogram, maximum radius\n"
            <<"    -scale f             : (optional) scaling factor\n"
            <<"                           (default 1e-4)\n"
            <<"    -maxr f              : (optional) max radius of input and output\n"
            <<"    -a                   : (optional) z-rotation, first\n"
            <<"    -b                   : (optional) y-rotation, second\n"
            <<"    -c                   : (optional) x-rotation, third\n"
            <<"    -xfigscale f         : (optional) XFig radius scaling factor\n"
            <<"                           (default 1.0)\n"
            <<"    -il f0 [f1 [.]][-li] : (optional) interval list\n"
            <<"    -cluster f           : (optional) maximum cluster radius\n"
            <<"    -ir f                : (optional) minimum radius interval\n"
            <<"    -cr f                : (optional) maximum correlation radius\n"
            <<"    -ci f0 f1            : (optional) correlation hemisphere interval\n"
            <<"    -color               : (default, optional) XFig, color output \n"
            <<"    -x                   : (optional) x-plane of input\n"
            <<"    -y                   : (optional) y-plane of input\n"
            <<"    -xy                  : (optional) xy-plane of input\n"
            <<"    -xz                  : (optional) yz-plane of input\n"
            <<"    -yz                  : (optional) yz-plane of input\n"
            <<"    -xyz                 : (optional) xyz-plane of input\n"
            <<"    -bw                  : (optional) XFig, b&w output \n"
            <<"    -scor                : (optional) Use spherical norm for correaltion. \n"
            <<"                         : (default: Euclidean norm\n"
            <<"    inputname            : (optional) filename of input\n"
            <<"                           (default: outputname.out.fin.pos.xyz)\n"
            <<"    outputname           : filename of outputs\n";
}
