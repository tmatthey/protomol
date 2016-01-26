/*  ------------------------------------------------------------------------  *
 *  To compile, modify one of the following commands to reflect your local
 *  OpenGL installation.  Don't forget that you need to define HAVE_GLUT in your
 *  environment or in the Makefile.
 *  
 *  g++ -I/usr/X11R6/include -L/usr/X11R6/lib -I/Home/romhegg/matthey/packages/glut-3.7/include -L/Home/romhegg/matthey/packages/glut-3.7/lib/glut -Wall xyzviz.cpp -O9 -ffast-math -finline-functions -funroll-loops -DNDEBUG -lglut -lGLU -lGL -lXmu -lXt -lSM -lICE -lXext -lX11 -lXi -lXext -lX11 -lm  -o xyzviz

 *  g++ -Wall xyzviz.cpp -O9 -ffast-math -finline-functions -funroll-loops -DNDEBUG -I$HOME/include/ -L$HOME/lib/ -L/usr/X11R6/lib -lglut -lGLU -lGL -lm -lXmu -lXi -lX11 -o xyzviz

 *  xlC xyzviz.cpp -O2 -I/usr/lpp/OpenGL/include/ -L/usr/lpp/OpenGL/lib/ -lglut -lGLU -lGL -lm -lXmu -lXi -lX11 -o xyzviz
 *  
 *  ------------------------------------------------------------------------  */

#ifndef HAVE_GLUT

#include <iostream>
int main(int argc, char **argv) {

    std::cout << "\nWarning: HAVE_GLUT is undefined.  See source code for more "
              << "information.\n";

    return(1);

}

#else

#include "DCDTrajectoryReader.h"
#include "PPMWriter.h"
#include "PGMWriter.h"
#include "PDBReader.h"
#include "XYZBinReader.h"
#include "XYZReader.h"
#include "XYZTrajectoryReader.h"
#include "mathutilities.h"
#include "systemutilities.h"
#include "Report.h"
#include "Timer.h"
#include "openglutilities.h"
#include "Matrix3by3.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>

using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::vector;
using std::string;
using std::map;

//_____________________________________________________________________ Opaque
struct Opaque {
  Opaque():i(0),x(0.0),y(0.0),z(0.0),o(0.0){}
  Opaque(int a, Vector3D b, double c):i(a),x(b.x),y(b.y),z(b.z),o(c){}
  int i;
  double x,y,z,o;
};

//_____________________________________________________________________ Feedback3Dcolor
struct Feedback3Dcolor {
  GLfloat x;
  GLfloat y;
  GLfloat z;
  GLfloat red;
  GLfloat green;
  GLfloat blue;
  GLfloat alpha;
};

//_____________________________________________________________________ DepthIndex
struct DepthIndex {
  DepthIndex():ptr(NULL),depth(0.0){}
  DepthIndex(GLfloat *a, GLfloat b):ptr(a),depth(b){}

  GLfloat *ptr;
  GLfloat depth;
};


//_____________________________________________________________________ Enum
enum MENU {ANIM_IDLE=-1,
	   PRINT_BUFFER=0,
	   SAVE_FRAMEBUFFER,
	   SAVE_FRAMEBUFFER_EPS,
	   QUIT,
	   TOGGLE_VIEW_CUBE,
	   TOGGLE_PLANE_ROTATION,
	   POLYGON_MODE,
	   TOGGLE_LIGHTING,
	   ANIM_SAVE,
	   ANIM_ON,
	   ANIM_OFF};

enum POLYMODE {FILL=0,
	       LINE,
	       POINT,
	       LAST}; // Just for modulo

//_____________________________________________________________________ Forward declaration
static bool cmpOpaque (const Opaque& a1, const Opaque& a2);

static bool nextFrame();
static void animationCallback(int value);
static void bmovieCallback(int);
static void bmovieCallback(long value);
static void display(void);
static void doAnimation(bool flag);
static void frameRateAnimationCallback(int value);
static void getFrame(int i);
static void keyboardCallback(unsigned char key, int x, int y);
static void mainMenuCallback(int value);
static void motion(int x, int y);
static void mouse(int button, int state, int x, int y);
static void parseCmdLine(int argc, 
			 char* argv[], 
			 string& inputname, 
			 string& outputname, 
			 double& scale, 
			 double& size, 
			 int& d, 
			 double& from, 
			 double& to,
			 Vector3D& rot,
			 int& tracer,
			 int& startFieldOfView,
			 double& fieldOfView,
			 GLfloat& angp,
			 GLfloat& angp2,
			 bool& batch,
			 double& dfr,
			 Vector3D& bbmin,
			 Vector3D& bbmax);
static void readXYZ(const string& inputname, const Vector3D& rot, int tracer, Vector3D b1, Vector3D b2);
static void reshape(int w, int h);
static void saveFramebufferEPSToFile();
static void saveFramebufferToFile();
static void startAnimation();
static void stopAnimation();
static void timerCallback(int value);
static void timerCallback(long value);
static void timerIdleCallback();
static void updateLighting();
static void updatePolygonMode();
static void usage(char* name);
static void verbose();
static void visibleParticlesCallback(int value);
static double brightness();
static Timer timer;

//_____________________________________________________________________ Statics
static bool         moving=false;
static GLfloat      angle = 0;      // in degrees 
static GLfloat      angle2 = 0;     // in degrees 
static GLfloat      anglePlane = 0;      // in degrees 
static GLfloat      anglePlane2 = 0;     // in degrees 
static int          startx=0;
static int          starty=0;
static bool         movingFieldOfView=false;
static int          startFieldOfView=0;
static double       fieldOfView = 80;
static double       xL=0;
static double       yL=0;
static double       zL=0;
static double       xT=0;
static double       yT=0;
static double       zT=0;
static double       maxL=0;
static double       scale = 1.0;
static double       opaque = 0.6;
static bool         viewCube = true;
static bool         planeRotation = false;
static GLsizei      winX = 600;
static GLsizei      winY = 600;
static double       radius = 0;
static int          currentFrame = 0;
static POLYMODE     polygonMode = FILL;
static bool         lighting = true;
static GLint        bufferSize = 0;
static bool         batch = false;
static double       dFocusRadius = 0.0;
static int subdivisions = 10;
static GLfloat lineWidth = 1.0;
static GLfloat pointSize = 1.0;

static bool saveAnimFrame = false;
static int timeout        = -1;
static int frame          = 0;
static int frameRate      = -1;
static int frameStep      = 1;
static int totalFrames    = 1;
static bool movie         = false;
static int animation;
static int mainMenu;
static const int EPSILON_TIME  = 5;  // [ms]
static int frameTimeN = 0;
static double from = 0.0;
static double to   = 0.0;
static string                    outputname = "Noname";
static vector<Vector3D>          r;
static vector<Vector3D>          rTmp;
static vector<Vector3D>          rAll;
static vector<string>       name;
static map<string,Vector3D> colorMap;
static map<string,int>      color2index;
static vector<string>       index2color;
static map<int,bool>             colorVisible;
static vector<Opaque>            opaqueSort;

//_____________________________________________________________________ brightness
double brightness(){
  GLint viewport[4]; // -> x, y, width, height
  glGetIntegerv(GL_VIEWPORT, viewport);
  GLint width  = viewport[2]-viewport[0];
  GLint height = viewport[3]-viewport[1];
  
  long size = width*height*3;
  GLfloat* pixels = new GLfloat[size];
  GLint drawBuffer;
  glGetIntegerv(GL_DRAW_BUFFER, &drawBuffer);
  glReadBuffer((GLenum)drawBuffer);  // read from correct buffer
  glReadPixels(viewport[0], viewport[1], width, height, GL_RGB, GL_FLOAT, pixels);
  double count = 0.0;
  for(long i=0;i<size;i++)
    count += std::max(std::min((double)pixels[i],1.0),0.0);
  delete [] pixels;
  return count/size;
}

//_____________________________________________________________________ cmpOpaque
bool cmpOpaque (const Opaque& a1, const Opaque& a2){
  if(a1.o > a2.o)
    return true;
  if(a1.o < a2.o)
    return false;
  if(a1.z < a2.z)
    return true;
  return false;
}

//_____________________________________________________________________ saveFramebufferEPSToFile
void saveFramebufferEPSToFile(){

  char postFix[20];
  sprintf(postFix,".%06d.eps",currentFrame);
  string filename = outputname + string(postFix);

  std::ofstream output(filename.c_str(), std::ios::out);
  if (!output){
    report << recoverable << "Attempt to open EPS file \'" << filename
	   << "\' for output failed."<<endr;
    return;
  }

  unsigned int count = openglToEPS(output,display);

  output.close();
  report << plain << "EPS \'" << filename << "\' created containing "<<count<<" OpenGL primitives."<<endr;
  currentFrame++;
}

//_____________________________________________________________________ updatePolygonMode
void updatePolygonMode(){
  switch (polygonMode) {
  case POINT:
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    break;
  case LINE:
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;
  case FILL:
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;
  case LAST:
  default:
    break;
  }
}

//_____________________________________________________________________ updateLighting
void updateLighting(){
  if (lighting) {
    glEnable(GL_LIGHTING);
  } else {
    glDisable(GL_LIGHTING);
  }
}


//_____________________________________________________________________ saveFramebufferToFile
void saveFramebufferToFile(){

  char postFix[20];
  sprintf(postFix,".%06d.ppm",currentFrame);
  string filename = outputname + string(postFix);

  if(dFocusRadius > 0.0){
    PGMWriter writer(filename);
    PGM pgm;
    openglToPGM(pgm);
    writer << pgm;
  }
  else {
    PPMWriter writer(filename);
    PPM ppm;
    openglToPPM(ppm);
    writer << ppm;
  }
  
  report << plain << "PPM \'" << filename << "\' created, frame "<<frame<<"."<<endr;
  
  currentFrame++;
}

//_____________________________________________________________________ readXYZ
void readXYZ(const string& inputname, const Vector3D& rot, int tracer, Vector3D bbmin, Vector3D bbmax){
  
  XYZ xyz;
  vector<XYZ> trajectory;

  // PDB
  if(PDBReader(inputname).tryFormat()){
    PDBReader reader(inputname);
    PDB tmp;
    if(!(reader >> tmp))
      report << error << "Could not read PDB file \'"<<inputname<<"\'."<<endr;
    xyz.clear();
    xyz.coords = tmp.coords;
    for(unsigned int j=0;j<tmp.atoms.size();++j)
      xyz.names.push_back(tmp.atoms[j].elementName);			
    trajectory.push_back(xyz);
  }
  // XYZ
  else if(XYZReader(inputname).tryFormat()){
    XYZReader reader(inputname);
    if(!(reader >> xyz))
      report << error << "Could not read XYZ file \'"<<inputname<<"\'."<<endr;
	    trajectory.push_back(xyz);
  }
  // Binary XYZ
  else if(XYZBinReader(inputname).tryFormat()){
    XYZBinReader reader(inputname);
    if(!(reader >> xyz))
      report << error << "Could not read binary XYZ file \'"<<inputname<<"\'."<<endr;
    trajectory.push_back(xyz);
  }
  // DCD
  else if(DCDTrajectoryReader(inputname).tryFormat()){
    DCDTrajectoryReader reader(inputname);
    while((reader >> xyz)){
      trajectory.push_back(xyz);
    }
    //const int n = (argc >1 ? atoi(argv[2]):0);
  }
  // Trajectory XYZ 
  else if(XYZTrajectoryReader(inputname).tryFormat()){
    XYZTrajectoryReader reader(inputname);
    while((reader >> xyz)){
      trajectory.push_back(xyz);
    }
  }
  else if(!isAccessible(inputname)){
    report << error<< "Can not open \'"<<inputname<<"\'."<<endr;
  }
  else {
    report << error<< "Can not figure out format of \'"<<inputname<<"\', skipping."<<endr;
  }

  movie = (trajectory.size() > 1);
  totalFrames = trajectory.size();

  vector<Vector3D> colors;
  if(from < to || dFocusRadius > 0){
    colors.push_back(Vector3D(0.0,0.0,0.0));
  }
  else {
    colors.push_back(Vector3D(1.0,0.0,0.0));
    colors.push_back(Vector3D(0.0,1.0,0.0));
    colors.push_back(Vector3D(0.0,0.0,1.0));
    colors.push_back(Vector3D(0.0,1.0,1.0));
    colors.push_back(Vector3D(1.0,1.0,0.0));
    colors.push_back(Vector3D(1.0,0.0,1.0));
    colors.push_back(Vector3D(0.0,0.0,0.0));
    colors.push_back(Vector3D(0.7,0.7,0.7));
  }

  colorMap.clear();
  color2index.clear();
  index2color.clear();
  colorVisible.clear();
  r.clear();
  rAll.clear();
  Matrix3by3 mat;
  mat.rotate(rot,Vector3D(0.0,0.0,1.0));

  Vector3D a,b;
  for(unsigned int i=0;i<trajectory.size();++i){
    for(unsigned int j=0;j<trajectory[i].size();++j){
      trajectory[i].coords[j] = mat*trajectory[i].coords[j];
      Real x = trajectory[i].coords[j].x;
      Real y = trajectory[i].coords[j].y;
      Real z = trajectory[i].coords[j].z;
      string tempstr = trajectory[i].names[j];

      if(i==0&&(j==0||x<a.x)) a.x=x;
      if(i==0&&(j==0||y<a.y)) a.y=y;
      if(i==0&&(j==0||z<a.z)) a.z=z;
      if(i==0&&(j==0||x>b.x)) b.x=x;
      if(i==0&&(j==0||y>b.y)) b.y=y;
      if(i==0&&(j==0||z>b.z)) b.z=z;

      if(tracer > 1 && (i + (trajectory[0].size()%tracer)) % tracer == 0)
	tempstr += "-Tracer";

      name.push_back(tempstr);
      rAll.push_back(Vector3D(x,y,z));
      if(colorMap.find(tempstr) == colorMap.end()){
	int index = colorMap.size() %  colors.size();
	colorMap[tempstr]=colors[index];
	color2index[tempstr]=index;
	index2color.push_back(tempstr);
	colorVisible[index]=true;
      }
    }
  }
  Real count;
  if((bbmin-bbmax).normSquared() > Constant::EPSILON){
    report << plain <<"Org BB:" << a << "-"<< b<<endr;
    bbmin = mat * bbmin;
    bbmax = mat * bbmax;
    a.x = std::min(bbmin.x,bbmax.x);
    a.y = std::min(bbmin.y,bbmax.y);
    a.z = std::min(bbmin.z,bbmax.z);
    b.x = std::max(bbmin.x,bbmax.x);
    b.y = std::max(bbmin.y,bbmax.y);
    b.z = std::max(bbmin.z,bbmax.z);
    count = 50.0;
  }
  report << plain <<"Use BB:" << a << "-"<< b<<endr;

  xL = b.x-a.x;
  yL = b.y-a.y;
  zL = b.z-a.z;
  xT = a.x;
  yT = a.y;
  zT = a.z;
  maxL = std::max(xL,std::max(yL,zL));
  for(unsigned int i = 0; i < rAll.size(); i++){
    rAll[i].x += xT;
    rAll[i].y += yT;
    rAll[i].z += zT;
  }

  r.resize(trajectory[0].size());
  count = r.size();
  radius = pow(maxL*maxL*maxL/count,1.0/3.0)/10.0*scale;
  getFrame(0);
}

//_____________________________________________________________________ parseCmdLine
void parseCmdLine(int argc, 
		  char* argv[], 
		  string& inputname, 
		  string& outputname, 
		  double& scale, 
		  double& size, 
		  int& d, 
		  double& f, 
		  double& t,
		  Vector3D& rot,
		  int& tracer,
		  int& startview,
		  double& fov,
		  GLfloat& angp,
		  GLfloat& angp2,
		  bool& batch,
		  double& dfr,
		  Vector3D& bbmin,
		  Vector3D& bbmax){
  if (argc < 2 || (argc >= 2 && (string(argv[1]) =="-h" ||
                                 string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  while (cur<(argc-1) && argv[cur][0]=='-') {

    string str(argv[cur]);

    if (str == "-scale" && cur<(argc-1)) {
      scale = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-plane" && cur<(argc-2)) {
      f = atof(argv[++cur]);
      t = atof(argv[++cur]);      
      cur++;
      continue;
    }

    if (str == "-view" && cur<(argc-2)) {
      startview = atoi(argv[++cur]);
      fov = atof(argv[++cur]);      
      cur++;
      continue;
    }

    if (str == "-angplane" && cur<(argc-2)) {
      angp = atof(argv[++cur]);
      angp2 = atof(argv[++cur]);      
      cur++;
      continue;
    }

    if (str == "-dfradius" && cur<(argc-1)) {
      dfr = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-batch") {
      batch = true;
      cur++;
      continue;
    }

    if (str == "-size") {
      size = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-tracer" && cur<(argc-1)) {
      tracer = atoi(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-subdivisions" && cur<(argc-1)) {
      d = atoi(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-z" || str == "-001"){
      rot = Vector3D(0,0,1);
      cur++;
      continue;
    }

    if (str == "-y" || str == "-010"){
      rot = Vector3D(0,1,0);
      cur++;
      continue;
    }

    if (str == "-x" || str == "-100"){
      rot = Vector3D(1,0,0);
      cur++;
      continue;
    }

    if (str == "-xy" || str == "-110"){
      rot = Vector3D(1,1,0);
      cur++;
      continue;
    }

    if (str == "-xz" || str == "-101"){
      rot = Vector3D(1,0,1);
      cur++;
      continue;
    }

    if (str == "-yz" || str == "-011"){
      rot = Vector3D(0,1,1);
      cur++;
      continue;
    }

    if (str == "-xyz" || str == "-111"){
      rot = Vector3D(1,1,1);
      cur++;
      continue;
    }

    if (str == "-rot" && cur<(argc-3)) {
      rot.x = atof(argv[++cur]);
      rot.y = atof(argv[++cur]);
      rot.z = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-box" && cur<(argc-2)) {
      bbmin.x = atof(argv[++cur]);
      bbmin.y = atof(argv[++cur]);      
      bbmin.z = atof(argv[++cur]);      
      bbmax.x = atof(argv[++cur]);
      bbmax.y = atof(argv[++cur]);      
      bbmax.z = atof(argv[++cur]);      
      cur++;
      continue;
    }

    break;

  }

  inputname = string(argv[cur++]);
  if (cur < argc){
    outputname = string(argv[cur]);
  }
}

//_____________________________________________________________________ usage
void usage(char* name){
  report << plain << "Usage: " << name << " [-h][--help] <options> inputname [outputname]\n"
	 << " where:\n"
	 << "    -h or --help             : This help output\n"
	 << "    -scale <real>            : (optional) sphere scaling factor\n"
	 << "    -plane <real> <real>     : (optional) z-plane projection\n"
	 << "    -view <real> <real>      : (optional) view\n"
	 << "    -angplane <real> <real>  : (optional) view angle\n"
	 << "    -dfradius <real>         : (optional) plane radius factor\n"
	 << "    -batch                   : (optional) batch mode, save all frames and quit\n"
	 << "    -size <real>             : (optional) point size and line width\n"
	 << "    -tracer <int>            : (optional) color every n'th differently\n"
	 << "    -subdivisions <int>      : (optional) number of subdivisions around/anlong the Z axis for sphere\n"
	 << "    -z or -001               : (optional) projection axis\n"
	 << "    -y or -010               : (optional) projection axis\n"
	 << "    -x or -100               : (optional) projection axis\n"
	 << "    -xy or -110              : (optional) projection axis\n"
	 << "    -xz or -101              : (optional) projection axis\n"
	 << "    -yz or -011              : (optional) projection axis\n"
	 << "    -xyz or -111             : (optional) projection axis\n"
	 << "    -rot <vector>            : (optional) rotate to (0,0,1)\n"
	 << "    -box <r><r><r> <r><r><r> : (optional) boundary box\n"
	 << "    inputname                : filename of input\n"
	 << "    outputname               : (optional) filename of EPS/PPM output"<<endr;
}

//_____________________________________________________________________ display
void display() {
  double fact = dFocusRadius*dFocusRadius/radius/radius;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  glTranslatef(-xT, -yT, -5 * maxL);
  glTranslatef(xT, yT, zT);
  glRotatef(angle2, 1.0, 0.0, 0.0);
  glRotatef(angle, 0.0, 1.0, 0.0);
  glTranslatef(-xT, -yT, -zT);
  glPushMatrix();
  glTranslatef(xT, yT, zT);
  if(viewCube){
    glColor3f(0.0, 0.0, 0.0);
    glutWireCube(maxL);
  }
  glPopMatrix();
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
#ifdef NO_OPAQUE_SORT 
  for(unsigned int n = 0; n < r.size(); n++) {
    if(colorVisible[color2index[name[n]]]){
      //std::cerr << r[n].x << "," << r[n].y << "," << r[n].z << "," << radius << "\n";
      if(from >= to && dFocusRadius <= 0.0){
	glPushMatrix();
	glTranslatef(r[n].x, r[n].y, r[n].z);
	glColor4f(colorMap[name[n]].x, colorMap[name[n]].y, colorMap[name[n]].z, opaque);
	glutSolidSphere(radius, subdivisions, subdivisions);
	glPopMatrix();
      }
      else {
	double x = r[n].x-xT;
	double y = r[n].y-yT;
	double z = r[n].z-zT;
	double a1=cos(anglePlane/180.0*M_PI);
	double a2=sin(anglePlane/180.0*M_PI);
	double b1=cos(anglePlane2/180.0*M_PI);
	double b2=sin(anglePlane2/180.0*M_PI);
	double u =  a1*x   + a2*b1*y  -a2*b2*z;
	double v = -a2*x   + a1*b1*y  -a1*b2*z;
	double w =             -b2*y  +   b1*z;
	x = u+xT;
	y = v+yT;
	z = w+zT;
	if(dFocusRadius > 0.0){      
	  double op = 1.0/(1.0+fact*z*z);
	  glPushMatrix();
	  glTranslatef(x, y, z);
	  glColor4f(colorMap[name[n]].x, colorMap[name[n]].y, colorMap[name[n]].z, opaque*op);
	  glutSolidSphere(radius/sqrt(op), subdivisions, subdivisions);
	  glPopMatrix();	  
	}
	else if(z >= from && z <= to){
	  double op = z-from;
	  op = 1.0-fabs(2.0*op/(to-from)-1.0);
	  //std::cerr <<from << "," << to << "," << r[n].z << ","<<op<< "\n";
	  glPushMatrix();
	  glTranslatef(x, y, z);
	  glColor4f(colorMap[name[n]].x, colorMap[name[n]].y, colorMap[name[n]].z, opaque*op);
	  glutSolidSphere(radius, subdivisions, subdivisions);
	  glPopMatrix();
	}
      }
    }
  }
#else
  if(opaqueSort.size() != r.size())
    opaqueSort.resize(r.size());
  if(rTmp.size() != r.size())
    rTmp.resize(r.size());
  unsigned int count = 0;
  for(unsigned int n = 0; n < r.size(); n++) {
    if(colorVisible[color2index[name[n]]]){
      if(from >= to && dFocusRadius <= 0.0){
	opaqueSort[count]= Opaque(n,r[n],1.0);
	++count;
      }
      else {

	// Edit here
	double x = r[n].x-xT;
	double y = r[n].y-yT;
	double z = r[n].z-zT;
	double u,v;
	u = y*cos(anglePlane2/180.0*M_PI) - z*sin(anglePlane2/180.0*M_PI);
	v = z*cos(anglePlane2/180.0*M_PI) + y*sin(anglePlane2/180.0*M_PI);
	y = u;
	z = v;
	u = x*cos(anglePlane/180.0*M_PI) - z*sin(anglePlane/180.0*M_PI);
	v = z*cos(anglePlane/180.0*M_PI) + x*sin(anglePlane/180.0*M_PI);
	x = u;
	z = v;
// 	double a1=cos(anglePlane/180.0*M_PI);
// 	double a2=sin(anglePlane/180.0*M_PI);
// 	double b1=cos(anglePlane2/180.0*M_PI);
// 	double b2=sin(anglePlane2/180.0*M_PI);
// 	double u =     b1*x   +    b2*y  ;
// 	double v = -a1*b2*x   + a1*b1*y + a2*z;
// 	double w =  a2*b2*x   - a2*b1*y + a1*z;
// 	x = u;
// 	y = v;
// 	z = w;
	double z0 = z;
	x += xT;
	y += yT;
	z += zT;
	if(dFocusRadius > 0.0) {
	  double op = 1.0/(1.0+fact*z0*z0);
	  if(op > 1e-2){
	    opaqueSort[count]=Opaque(n,Vector3D(x,y,z),op);
	    ++count;
	  }
	}
	else if(z0 >= from && z0 <= to){
	  double op = z0-from;
	  op = 1.0-fabs(2.0*op/(to-from)-1.0);
	  opaqueSort[count]=Opaque(n,Vector3D(x,y,z),op);
	  ++count;
	}
      }
    }
  }
  std::sort(opaqueSort.begin(),opaqueSort.begin()+count,cmpOpaque);
  for(unsigned int i = 0; i < count; i++) {
    int n = opaqueSort[i].i;
    double op = opaqueSort[i].o*opaque;
    double rad = radius;
    if(dFocusRadius > 0.0){
      rad /= sqrt(opaqueSort[i].o);
    }
    glPushMatrix();
    glTranslatef(opaqueSort[i].x,opaqueSort[i].y,opaqueSort[i].z);
    glColor4f(colorMap[name[n]].x, colorMap[name[n]].y, colorMap[name[n]].z, op);
    glutSolidSphere(rad, subdivisions, subdivisions);
    glPopMatrix();
  }
#endif
  glDisable(GL_BLEND);
  glPopMatrix();
  glutSwapBuffers();
}


//_____________________________________________________________________ reshape
void reshape (int w, int h) {
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fieldOfView, w / static_cast<float>(h), 0.1*maxL,  10*maxL);
  winX = w;
  winY = h;
}

//_____________________________________________________________________ mouse
void mouse (int button, int state, int x, int y) {
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      movingFieldOfView = true;
      startFieldOfView = y;
    }
    if (state == GLUT_UP) {
      movingFieldOfView = false;
    }
    break;
  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) {
      moving = true;
      startx = x;
      starty = y;
    }
    if (state == GLUT_UP) {
      moving = false;
    }       
    break;

  case GLUT_RIGHT_BUTTON:
    break;
  default:
    break;
  }
}

//_____________________________________________________________________ motion
void motion(int x, int y) {
  if (moving) {
    if((from < to || dFocusRadius > 0) && planeRotation){
      anglePlane = anglePlane + (x - startx);
      anglePlane2 = anglePlane2 + (y - starty);
    }
    else {
      angle = angle + (x - startx);
      angle2 = angle2 + (y - starty);
    }
    startx = x;
    starty = y;
    glutPostRedisplay();
  }
  if (movingFieldOfView) {
    fieldOfView = fieldOfView + (y - startFieldOfView)/10.0;
    startFieldOfView = y;
    if(fieldOfView < 10){
      fieldOfView = 10;
    }
    if(fieldOfView > 170){
      fieldOfView = 170;
    }    
    //std::cerr << fieldOfView << "\n";
    reshape(winX,winY);
    glutPostRedisplay();
  }
}

void mainMenuCallback(int value){
  switch (value){
  case TOGGLE_VIEW_CUBE:
    viewCube = !viewCube;
    glutPostRedisplay();    
    break;
  case TOGGLE_PLANE_ROTATION:
    planeRotation = !planeRotation;
    glutSetMenu(mainMenu);
    glutChangeToMenuEntry(2,const_cast<char*>((string("Toggle Plane Rotation : ")+(planeRotation?string("Yes"):string("No"))).c_str()),TOGGLE_PLANE_ROTATION);
    break;
  case QUIT:
    exit(0);
    break;
  case SAVE_FRAMEBUFFER:
    saveFramebufferToFile();
    break;
  case PRINT_BUFFER:
    openglToPlain(std::cerr,display);
    break;
  case SAVE_FRAMEBUFFER_EPS:
    saveFramebufferEPSToFile();

    break;
  case POLYGON_MODE:
    polygonMode = static_cast<POLYMODE>((static_cast<int>(polygonMode) + 1) % static_cast<int>(LAST));
    updatePolygonMode();
    glutPostRedisplay();
    break;
  case TOGGLE_LIGHTING:
    lighting = !lighting;
    updateLighting();
    glutPostRedisplay();
    break;
  default:
    break;
  }
}

//_____________________________________________________________________ keyboardCallback
void animationCallback(int value) {
  switch (value) {
  case ANIM_OFF:        // Off
    stopAnimation();
    break;
  case ANIM_SAVE:
    if (saveAnimFrame){    // no-save frames
      saveAnimFrame = false;
      glutSetMenu(animation);
      glutChangeToMenuEntry(5,"Save frames : No",ANIM_SAVE);
    } 
    else {                   // save frames
      saveAnimFrame = true;
      glutSetMenu(animation);
      glutChangeToMenuEntry(5,"Save frames : Yes",ANIM_SAVE);
      glutPostRedisplay();
      saveFramebufferToFile();
    }
    break;
  case ANIM_ON:
    startAnimation();
    break;
  default:  
    break;
  }  
}
//_____________________________________________________________________ stopAnimation
void bmovieCallback(long value){
  bmovieCallback(static_cast<int>(value));
}
void bmovieCallback(int v){
  saveFramebufferToFile();
  startAnimation();
}
//_____________________________________________________________________ stopAnimation
void stopAnimation(){
  timer.stop();
  if(timer.getTime().getRealTime() > 0 && frameTimeN > 0)
    report << plain <<"Frames/sec: "<<static_cast<double>(frameTimeN)/timer.getTime().getRealTime()<<"."<<endr;
  timeout = -1;
  glutIdleFunc(0);
  if(batch)
    exit(0);
}
//_____________________________________________________________________ startAnimation
void startAnimation(){
  timer.stop();
  timer.reset();
  timer.start();
  frameTimeN = 0;
  if (frameRate > 0){                      // Timeout
    timeout  = (frameRate > EPSILON_TIME) ? frameRate - EPSILON_TIME : 1;
    glutIdleFunc(0);
    glutTimerFunc(timeout,&timerCallback,1);
  }
  else {   // Idle, as fast as possible
    timeout = 0;
    glutIdleFunc(timerIdleCallback);
  }
}

//_____________________________________________________________________ frameRateAnimationCallback
void frameRateAnimationCallback(int value){
  frameRate = value;
  if(timeout >= 0)
    startAnimation();
}

//_____________________________________________________________________ visibleParticlesCallback
void visibleParticlesCallback(int value){
  colorVisible[value] = !colorVisible[value];
  display();
}

//_____________________________________________________________________ timerCallback
void timerCallback(long value){
  timerCallback(static_cast<int>(value));
}
void timerCallback(int value){
  Timer t;
  t.start();
  doAnimation(false);   
  t.stop();
  int n = static_cast<int>(floor(t.getTime().getRealTime()*1000.0));
  if (timeout>0){
    if(timeout-n > 0)
      glutTimerFunc(timeout-n,timerCallback,value + 1);  
    else
      glutTimerFunc(1,timerCallback,value + 1);  
  }
}
//_____________________________________________________________________ timerIdleCallback
void timerIdleCallback(){
  doAnimation(false);  
}

//_____________________________________________________________________ doAnimation
void doAnimation(bool flag){
  if (timeout >= 0 || flag){
    if(nextFrame()){
      display();
      if (saveAnimFrame)
	saveFramebufferToFile();
    }
    else {
      if(timeout >= 0){
	report << plain << "Last frame ("<<frame<<")."<<endr;
      }
      stopAnimation();
    }
  }
  else {
    display();
  }
}

//_____________________________________________________________________ nextFrame
bool nextFrame(){
  int j = frame;
  frame += frameStep;
  int i = frame;

  if(i < 0)
    i = 0;
  else if(i >= totalFrames){
    i = totalFrames -1;
  }
  
  if(j < 0)
    j = 0;
  else if(j >= totalFrames){
    j = totalFrames -1;
  }

  bool update = (i != j);
  if(update){
    getFrame(i);
    frameTimeN++;
  }
  frame = i;
  return update;
}
//_____________________________________________________________________ getFrame
void getFrame(int i){
  std::copy(&(rAll[i*r.size()]), &(rAll[(i+1)*r.size()]), r.begin());
}

//_____________________________________________________________________ frameStepAnimationCallback
void frameStepAnimationCallback(int value){
  frameStep = value;
}
 
//_____________________________________________________________________ verbose
void verbose(){
  report << plain <<"Particle(s)         : " << r.size() << "\n";
  report  <<"Frame(s)            : " << totalFrames << "\n";
  report  <<"Box                 : [" << xL << "," << yL << "," << zL << "]\n";
  report  <<"Types(s)            : " << colorMap.size() << "\n";
  report  <<"Sphere radius       : " << radius << "\n";
  report  <<"Sphere subdivisions : " << subdivisions << "\n";
  report  <<"Scaling             : " << scale << "\n";
  report  <<"Opaque              : " << opaque << "\n";
  report  <<"Point/line size     : " << lineWidth << "\n";
  report  <<"View cube           : " << (viewCube?string("on"):string("off")) << "\n";
  report  <<"Lighting            : " << (lighting?string("on"):string("off")) << "\n";
  report  <<"Switch fill mode    : " << (polygonMode==FILL?string("fill"):(polygonMode==LINE?string("line"):string("point"))) << "\n";
  report  <<"Window              : " << winX<<"x"<<winY<<", angles ("<<angle<<", "<<angle2<<") angles plane ("<<anglePlane<<", "<<anglePlane2<<"), start ("<<startx<<", "<<starty<<"), "
	    << "view ("<<startFieldOfView<<", "<<fieldOfView<<"), buffer: "<<bufferSize<<".\n";
  report  <<"Picture number      : " << currentFrame << "\n";
  report  <<"Animation           : frame: " << frame <<", frame step: " << frameStep <<", frames/sec:";
  if(frameRate > 0)
    report  <<1000.0/static_cast<double>(frameRate);
  else
    report  <<"max";
  if(timer.getActualTime().getRealTime() > 0 && frameTimeN > 0)
    std::cerr <<" showing at frames/sec: "<<static_cast<double>(frameTimeN)/timer.getActualTime().getRealTime();
  std::cerr <<".\n";
  if(from < to) {
    report  <<"Projection plane    : " <<from<<" - "<<to<<endr;
  }
  else if(dFocusRadius > 0) {
    report  <<"Projection focus    : " <<dFocusRadius<<std::endl;
    report  <<"Max radius scaling  : " <<sqrt(1.0+dFocusRadius*dFocusRadius/radius/radius*maxL*maxL/4.0)<<std::endl;    
    report  <<"Brightness          : " <<brightness()<<endr;
  }
  else {
    report  <<"3D"<<endr;
  }
}
//_____________________________________________________________________ keyboardCallback
void keyboardCallback(unsigned char key, int, int){
  switch (key) {
  case 'z':
    {
      double b = 0.0;
      GLfloat a1 = anglePlane;
      GLfloat a2 = anglePlane2;
      GLfloat b1 = anglePlane;
      GLfloat b2 = anglePlane2;
      
      for(int i = -5;i <= 5;i++){
	for(int j = -5;j <= 5;j++){
	  anglePlane = b1+i;
	  anglePlane2 = b2+j;
	  display();
	  double a = brightness();
	  if(a > b){
	    a1 = anglePlane;
	    a2 = anglePlane2;
	    b = a;
	  }
	}
      }
      anglePlane = a1;
      anglePlane2 = a2;
      display();
    }
    break;
  case 'w':                   // Save framebuffer PPM
    saveFramebufferToFile();
    break;
  case 'p':                   // Save framebuffer EPS
    saveFramebufferEPSToFile();
    break;
  case 27:    // ESC and 'q' will quit.
  case 'q': 
    exit(1); 
    break;
  case '-':
    opaque -= 0.05;
    if(opaque < 0.0)
      opaque = 0.0;
    else
      glutPostRedisplay();
    break;
  case '+':
    opaque += 0.05;
    if(opaque > 1.0)
      opaque = 1.0;
    else
      glutPostRedisplay();
    break;
  case 'r':
    planeRotation = !planeRotation;
    glutPostRedisplay();
    break;
  case ' ':
    viewCube = !viewCube;
    glutPostRedisplay();
    break;
  case 'v':
    verbose();
    break;
  case 'a':                    // Animation statistic
    if(movie){
      report << plain <<"Frame: " << frame <<", frame step: " << frameStep <<", frames/sec:";
      if(frameRate > 0)
	report  <<1000.0/static_cast<double>(frameRate);
      else
	report  <<"max";
      report  <<"."<<endr;

    }
    break; 
  case 's':                    // Stop animation
    if(movie){
      stopAnimation();
    }
    break; 
  case 'i':                   // Idle animation (as fast as possible)
    if(movie){
      startAnimation();
    }
    break;
 
  case 'f':                  // Stop animation and one step forward
    if(movie){
      stopAnimation();
      doAnimation(true);
    }
    break;
    
  case 'b':                 // Stop animation and one step back
    if(movie){
      stopAnimation();
      frameStep *= -1;
      doAnimation(true);
      frameStep *= -1;
    }
    break;
  case 'c':                // Stop and reset animation
    if(movie){
      stopAnimation();
      frame = 0;            // Go back to frame 0
      frameStep = 1;         // frameStep = 1
      getFrame(frame);
      glutPostRedisplay();
    }
    break;
    
  case 9:                  // Change sign of frameStep
    if(movie){
      frameStep *=-1;
    }
    break;
  case 'h':
  case '?':
    report << plain <<"Summary of keyboard commands:" << std::endl
	   << "  q, ESC :  Quit." << std::endl
	   << "  w      :  Write PPM file." << std::endl
	   << "  p      :  Write EPS file." << std::endl
	   << "  SPACE  :  Toggle view cube." << std::endl
	   << "  +, -   :  Increase/decrease colors." << std::endl
	   << "  h, ?   :  This message." << std::endl
	   << "  v      :  Debug information." << std::endl;
    if(from < to || dFocusRadius > 0.0)
      report << "  r      :  Toggle plane rotation." << std::endl;
    if(movie){
      report << "  s      :  Stop animation." << std::endl
	     << "  i      :  Start animation." << std::endl
	     << "  f      :  Stop animation and one step forward." << std::endl
	     << "  b      :  Stop animation and one step back." << std::endl
	     << "  a      :  Animation statistic." << std::endl
	     << "  TAB    :  Change sign of frameStep." << std::endl
	     << "  c      :  Set actual frame = 0 and frameStep = 1." << endr;
    }
    
    break;
  default:
    break;
  }
}

//_____________________________________________________________________ main
int main (int argc, char** argv) {
  glutInit(&argc, argv);


  string inputname;
  double size = 0.0;
  int d = 0;
  Vector3D normal(1.0,0,0);
  Vector3D b1(1.0,0,0);
  Vector3D b2(1.0,0,0);
  int tracer = 0;
  parseCmdLine(argc,argv, inputname, outputname, scale, size, d, from, to, normal,tracer,startFieldOfView,fieldOfView,anglePlane,anglePlane2,batch,dFocusRadius,b1,b2);

  if(size > 0){
    lineWidth = size;
    pointSize = size;
  }
  if(d > 2){
    subdivisions = d;
  }

  if(batch || (from < to) || (dFocusRadius > 0.0))
    viewCube = false;    

  planeRotation = (from < to) || (dFocusRadius > 0.0);

  if(dFocusRadius > 0.0)
    opaque = 1.0;

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(winX, winY);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);

  GLfloat whiteLight[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat lightPosition1[] = {5, 5, 10.0, 0.0};

  glClearColor(1.0, 1.0, 1.0, 0.0);
  glShadeModel(GL_SMOOTH);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteLight);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition1);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);

  glLineWidth(lineWidth);
  glPointSize(pointSize);

  updateLighting();
  updatePolygonMode();

  glutMotionFunc(motion);
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutKeyboardFunc(keyboardCallback);

  readXYZ(inputname,normal,tracer,b1,b2);

  // Animation
  if(movie){
    int frameStepAnimation = glutCreateMenu(frameStepAnimationCallback);
    glutAddMenuEntry("-100",           -100);
    glutAddMenuEntry(" -50",            -50);
    glutAddMenuEntry(" -20",            -20);
    glutAddMenuEntry(" -10",            -10);
    glutAddMenuEntry("  -5",             -5);
    glutAddMenuEntry("  -2",             -2);
    glutAddMenuEntry("  -1",             -1);
    glutAddMenuEntry("   1",              1);
    glutAddMenuEntry("   2",              2);
    glutAddMenuEntry("   5",              5);
    glutAddMenuEntry("  10",             10);
    glutAddMenuEntry("  20",             20);
    glutAddMenuEntry("  50",             50);
    glutAddMenuEntry(" 100",            100);
    
    int frameRateAnimation = glutCreateMenu(frameRateAnimationCallback);
    glutAddMenuEntry("25",            40);
    glutAddMenuEntry("20",            50);
    glutAddMenuEntry("10",           100);
    glutAddMenuEntry(" 5",           200);
    glutAddMenuEntry(" 4",           250);
    glutAddMenuEntry(" 2",           500);
    glutAddMenuEntry(" 1",          1000);
    glutAddMenuEntry("as fast as possible",ANIM_IDLE);
  
    animation = glutCreateMenu(animationCallback);
    glutAddMenuEntry("Start animation", ANIM_ON);
    glutAddMenuEntry("Stop animation", ANIM_OFF);
    glutAddSubMenu  ("Frames/sec", frameRateAnimation);
    glutAddSubMenu  ("Framestep", frameStepAnimation);
    glutAddMenuEntry("Save frames : No",  ANIM_SAVE);
  }

  int visibleParticles = glutCreateMenu(visibleParticlesCallback);
  for(unsigned int i=0;i<index2color.size();i++)
    glutAddMenuEntry(const_cast<char*>(index2color[i].c_str()),i);    

  mainMenu = glutCreateMenu(mainMenuCallback);
  glutAddMenuEntry("Toggle view cube",TOGGLE_VIEW_CUBE);
  if(from < to || dFocusRadius > 0.0)
    glutAddMenuEntry(const_cast<char*>((string("Toggle Plane Rotation : ")+(planeRotation?string("Yes"):string("No"))).c_str()),TOGGLE_PLANE_ROTATION);
  glutAddMenuEntry("Toggle lighting", TOGGLE_LIGHTING);
  glutAddMenuEntry("Save framebuffer PPM", SAVE_FRAMEBUFFER);
  glutAddMenuEntry("Save framebuffer EPS", SAVE_FRAMEBUFFER_EPS);

  if(movie){
    glutAddSubMenu("Animation",animation);
  }
  glutAddSubMenu("Toggle particles",visibleParticles);
  glutAddMenuEntry("Print framebuffer", PRINT_BUFFER);
  glutAddMenuEntry("Switch fill mode (poly, line, point)", POLYGON_MODE);
  glutAddMenuEntry("Quit", QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  
  verbose();
  
  if(batch){
    saveAnimFrame = true;
    glutSetMenu(animation);
    glutChangeToMenuEntry(5,"Save frames : Yes",ANIM_SAVE);
    glutPostRedisplay();
    glutTimerFunc(1000,&bmovieCallback,1);
  }
  
  glutMainLoop();
  return 0;
}
#endif
