/*
  g++ -Wall xyzviz.cpp -O9 -ffast-math -finline-functions -funroll-loops -DNDEBUG -I$HOME/include/ -L$HOME/lib/ -L/usr/X11R6/lib -lglut -lGLU -lGL -lm -lXmu -lXi -lX11 -o xyzviz
  xlC xyzviz.cpp -O2 -I/usr/lpp/OpenGL/include/ -L/usr/lpp/OpenGL/lib/ -lglut -lGLU -lGL -lm -lXmu -lXi -lX11 -o xyzviz
*/

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>


//_____________________________________________________________________ Vector3D
struct Vector3D {
  Vector3D():x(0.0),y(0.0),z(0.0){}
  Vector3D(double a, double b, double c):x(a),y(b),z(c){}
  double x,y,z;
};
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
static bool cmpRadius (const Vector3D& a1, const Vector3D& a2);
static bool compare (const DepthIndex& a1, const DepthIndex& a2);
static bool nextFrame();
static double getTime();
static double getTime(double last);
static void animationCallback(int value);
static void bmovieCallback();
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
			 std::string& inputname, 
			 std::string& outputname, 
			 double& scale, 
			 double& size, 
			 int& d, 
			 bool& m, 
			 double& from, 
			 double& to,
                         double& alpha,
                         double& beta,
                         double& gamma,
			 int& tracer,
			 int& startFieldOfView,
			 double& fieldOfView,
			 GLfloat& angp,
			 GLfloat& angp2,
			 bool& bmovie,
			 double& dfr);
static void print3DcolorVertex(GLint size, GLint * count, GLfloat * buffer);
static void printBuffer();
static void readXYZ(const std::string& inputname, bool movie, double alpha, double beta, double gamma, int tracer);
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
static bool         bmovie = false;
static double       dFocusRadius = 0.0;
static int subdivisions = 10;
static GLfloat lineWidth = 1.0;
static GLfloat pointSize = 1.0;
static const double EPS_GOURAUD_THRESHOLD = 0.1;
static const double EPS_SMOOTH_LINE_FACTOR = 0.06;

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
static double frameTime = 0.0;
static int frameTimeN = 0;
static double from = 0.0;
static double to   = 0.0;
static std::string                    outputname = "Noname";
static std::vector<Vector3D>          r;
static std::vector<Vector3D>          rTmp;
static std::vector<Vector3D>          rAll;
static std::vector<std::string>       name;
static std::map<std::string,Vector3D> colorMap;
static std::map<std::string,int>      color2index;
static std::vector<std::string>       index2color;
static std::map<int,bool>             colorVisible;
static std::vector<Opaque>            opaqueSort;

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
//_____________________________________________________________________ cmpVector3D
bool cmpVector3D (const Vector3D& a1, const Vector3D& a2){
  return ((a1.x*a1.x+a1.y*a1.y+a1.z*a1.z)<(a2.x*a2.x+a2.y*a2.y+a2.z*a2.z));
}
//_____________________________________________________________________ getTime
double getTime(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (static_cast<double>(tv.tv_sec) + static_cast<double>(tv.tv_usec)*.000001);
}

//_____________________________________________________________________ getTime
double getTime(double last){
  return (getTime() - last);
}
//_____________________________________________________________________ compare
bool compare(const DepthIndex& a1, const DepthIndex& a2){
  if(a2.depth < a1.depth)
    return true;
  return false;
}

//_____________________________________________________________________ print3DcolorVertex
void print3DcolorVertex(GLint size, GLint * count, GLfloat * buffer){
  printf("  ");
  for (int i = 0; i < 7; i++) {
    printf("%4.2f ", buffer[size - (*count)]);
    *count = *count - 1;
  }
  printf("\n");
}

//_____________________________________________________________________ printBuffer
void printBuffer(){
  std::vector<GLfloat> feedbackBuffer(bufferSize);
  while(true){
    glFeedbackBuffer(feedbackBuffer.size(), GL_3D_COLOR, &(feedbackBuffer[0]));
    glRenderMode(GL_FEEDBACK);
    display();
    GLint returned = glRenderMode(GL_RENDER);
    if(returned > 0){
      feedbackBuffer.resize(returned);
      break;
    }
    bufferSize += 100*1024;
    feedbackBuffer.resize(bufferSize);
  }
  GLint size = feedbackBuffer.size();
  GLfloat* buffer = &(feedbackBuffer[0]);

  GLint count;
  int token, nvertices;

  count = size;
  while (count) {
    token = static_cast<int>(buffer[size - count]);
    count--;
    switch (token) {
    case GL_PASS_THROUGH_TOKEN:
      printf("GL_PASS_THROUGH_TOKEN\n");
      printf("  %4.2f\n", buffer[size - count]);
      count--;
      break;
    case GL_POINT_TOKEN:
      printf("GL_POINT_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_TOKEN:
      printf("GL_LINE_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_LINE_RESET_TOKEN:
      printf("GL_LINE_RESET_TOKEN\n");
      print3DcolorVertex(size, &count, buffer);
      print3DcolorVertex(size, &count, buffer);
      break;
    case GL_POLYGON_TOKEN:
      printf("GL_POLYGON_TOKEN\n");
      nvertices = static_cast<int>(buffer[size - count]);
      count--;
      for (; nvertices > 0; nvertices--) {
        print3DcolorVertex(size, &count, buffer);
      }
    default:
      break;
    }
  }
}

//_____________________________________________________________________ saveFramebufferEPSToFile
void saveFramebufferEPSToFile(){

  char postFix[20];
  sprintf(postFix,".%06d.eps",currentFrame);
  std::string filename = outputname + std::string(postFix);

  std::ofstream output(filename.c_str(), std::ios::out);
  if (!output){
    std::cerr << "Attempt to open EPS file \'" << filename
	      << "\' for output failed.\n";
    return;
  }


  //
  //
  //
  std::vector<GLfloat> feedbackBuffer(bufferSize);
  while(true){
    glFeedbackBuffer(feedbackBuffer.size(), GL_3D_COLOR, &(feedbackBuffer[0]));
    glRenderMode(GL_FEEDBACK);
    display();
    GLint returned = glRenderMode(GL_RENDER);
    if(returned > 0){
      feedbackBuffer.resize(returned);
      break;
    }
    bufferSize += 100*1024;
    feedbackBuffer.resize(bufferSize);
  }


  GLfloat clearColor[4], viewport[4];

  // Read back a bunch of OpenGL state to help make the EPS
  // consistent with the OpenGL clear color, line width, point
  // size, and viewport.
  glGetFloatv(GL_VIEWPORT, viewport);
  glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
  glGetFloatv(GL_LINE_WIDTH, &lineWidth);
  glGetFloatv(GL_POINT_SIZE, &pointSize);

  // Emit EPS header. 
  output << "%!PS-Adobe-2.0 EPSF-2.0\n";
  // Notice %% for a single % in the fprintf calls. 
  output << "%%Creator: " << filename << " (using OpenGL feedback)\n";
  output << "%%BoundingBox: " << viewport[0] << " " << viewport[1] << " " << viewport[2] << " " << viewport[3] << "\n";
  output << "%%EndComments\n";
  output << "\n";
  output << "gsave\n";
  output << "\n";

  // Output Frederic Delhoume's "gouraudtriangle" PostScript fragment. 
  output << "% the gouraudtriangle PostScript fragement below is free\n";
  output << "% written by Frederic Delhoume (delhoume@ilog.fr)\n";
  output << "/threshold " << EPS_GOURAUD_THRESHOLD << " def\n"
	 << "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3\n"
	 << "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd\n"
	 << "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2\n"
	 << "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse\n"
	 << "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get\n"
	 << "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get\n"
	 << "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get\n"
	 << "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll\n"
	 << "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2\n"
	 << "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4\n"
	 << "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add\n"
	 << "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1\n"
	 << "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3\n"
	 << "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll\n"
	 << "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll\n"
	 << "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array\n"
	 << "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17\n"
	 << "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index\n"
	 << "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6\n"
	 << "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index\n"
	 << "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14\n"
	 << "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3\n"
	 << "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7\n"
	 << "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add\n"
	 << "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd";



  output << "\n" << lineWidth << " setlinewidth\n";

  // Clear the background like OpenGL had it.
  output << clearColor[0] << " " << clearColor[1] << " " << clearColor[2] << " setrgbcolor\n";
  output << viewport[0] << " " << viewport[1] << " " << viewport[2] << " " << viewport[3] << " rectfill\n\n";


  // Collect all primitives
  Feedback3Dcolor *vertex;
  std::vector<DepthIndex> prims;
  GLfloat* loc = &(feedbackBuffer[0]);
  GLfloat* end = loc + feedbackBuffer.size();
  int count = 0;
  while (loc < end) {
    DepthIndex item;
    double depthSum;
    int nvertices;
    double area;

    item.ptr = loc;  // Save this primitive's location.
    int token = static_cast<int>(*loc);
    loc++;
    count++;
    switch (token) {
    case GL_LINE_TOKEN:
    case GL_LINE_RESET_TOKEN:
      vertex = (Feedback3Dcolor*) loc;
      depthSum = vertex[0].z + vertex[1].z;
      item.depth = depthSum / 2.0;
      prims.push_back(item);
      loc += 14;
      break;
    case GL_POLYGON_TOKEN:
      nvertices = static_cast<int>(*loc);
      loc++;
      vertex = (Feedback3Dcolor*) loc;
      depthSum = vertex[0].z;
      area = 0.0;
      for(int i = 1; i < nvertices; i++) {
        depthSum += vertex[i].z;
	area += (vertex[0].y- vertex[i-1].y)*(vertex[i].x-vertex[i-1].x)-(vertex[0].x- vertex[i-1].x)*(vertex[i].y-vertex[i-1].y);
      }
      item.depth = depthSum / nvertices;
      // Only add if triangle is visible and not degenerated
      if(area > 0.04)
	prims.push_back(item);
      loc += (7 * nvertices);
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor*) loc;
      item.depth = vertex[0].z;
      prims.push_back(item);
      loc += 7;
      break;
    default:
      std::cerr << token << " - opps a GL token I do not understand.\n";

      break;
    }
  }

  // Sort them
  std::sort(prims.begin(),prims.end(),compare);

  // Print them
  for(unsigned int j=0;j<prims.size();j++){
    GLfloat* loc = prims[j].ptr;
    //std::cerr <<j << ".\n";

    GLfloat red = 0.0, green = 0.0, blue = 0.0;
    GLfloat dx = 0.0, dy = 0.0, dr = 0.0, dg = 0.0, db = 0.0, absR = 0.0, absG = 0.0, absB = 0.0, colormax = 0.0;
    GLfloat xstep = 0.0, ystep = 0.0, rstep = 0.0, gstep = 0.0, bstep = 0.0;
    GLfloat xnext = 0.0, ynext = 0.0, rnext = 0.0, gnext = 0.0, bnext = 0.0, distance = 0.0;
    int steps,nvertices;

    int token = static_cast<int>(*loc);
    loc++;
    switch (token) {
    case GL_LINE_RESET_TOKEN:
    case GL_LINE_TOKEN:
      vertex = (Feedback3Dcolor*) loc;
      
      dr = vertex[1].red - vertex[0].red;
      dg = vertex[1].green - vertex[0].green;
      db = vertex[1].blue - vertex[0].blue;
      
      if (dr != 0 || dg != 0 || db != 0) {
	// Smooth shaded line.
	dx = vertex[1].x - vertex[0].x;
	dy = vertex[1].y - vertex[0].y;
	distance = sqrt(dx * dx + dy * dy);
	
	absR = fabs(dr);
	absG = fabs(dg);
	absB = fabs(db);
	

	colormax = std::max(absR, std::max(absG, absB));
	steps = static_cast<int>(std::max(1.0, static_cast<double>(colormax * distance) * EPS_SMOOTH_LINE_FACTOR));

	xstep = dx / steps;
	ystep = dy / steps;

	rstep = dr / steps;
	gstep = dg / steps;
	bstep = db / steps;

	xnext = vertex[0].x;
	ynext = vertex[0].y;
	rnext = vertex[0].red;
	gnext = vertex[0].green;
	bnext = vertex[0].blue;

	// Back up half a step; we want the end points to be
	//   exactly the their endpoint colors. 
	xnext -= xstep / 2.0;
	ynext -= ystep / 2.0;
	rnext -= rstep / 2.0;
	gnext -= gstep / 2.0;
	bnext -= bstep / 2.0;
      } else {
	// Single color line. 
	steps = 0;
      }

      output << vertex[0].red << " " << vertex[0].green << " " << vertex[0].blue << " setrgbcolor\n";
      output << vertex[0].x << " " << vertex[0].y << " moveto\n";

      for(int i = 0; i < steps; i++) {
	xnext += xstep;
	ynext += ystep;
	rnext += rstep;
	gnext += gstep;
	bnext += bstep;
	output << xnext << " " << ynext << "  lineto stroke\n";
	output << rnext << " " << gnext << " " << bnext << "   setrgbcolor\n";
	output << xnext << " " << ynext << "  moveto\n";
      }
      output << vertex[1].x << " " << vertex[1].y << "  lineto stroke\n";
      
      loc += 14;          // Each vertex element in the feedback
      //   buffer is 7 GLfloats. 
      
      break;
    case GL_POLYGON_TOKEN:
      nvertices = static_cast<int>(*loc);
      loc++;
      
      vertex = (Feedback3Dcolor*) loc;
      
      
      if (nvertices > 0) {
	red = vertex[0].red;
	green = vertex[0].green;
	blue = vertex[0].blue;
	bool smooth = (vertex[0].alpha >= 1.0);
	for(int i = 1; i < nvertices; i++) {
	  if (vertex[0].alpha < 1.0 || red != vertex[i].red || green != vertex[i].green || blue != vertex[i].blue) {
	    smooth = true;
	    break;
	  }
	}

	if (smooth) {
	  // Smooth shaded polygon; varying colors at vetices. 
	  int triOffset;
	  
	  // Break polygon into "nvertices-2" triangle fans. 
	  for(int i = 0; i < nvertices - 2; i++) {
	    triOffset = i * 7;
	    output << "[ " << vertex[0].x << " " << vertex[i + 1].x << " " << vertex[i + 2].x << " "
		   << vertex[0].y << " " << vertex[i + 1].y << " " << vertex[i + 2].y << " ] ";
	    output << "[ " << vertex[0].red << " " << vertex[0].green << " " << vertex[0].blue << " ] "
		   << "[ " << vertex[i + 1].red << " " << vertex[i + 1].green << " " << vertex[i + 1].blue << " ] "
		   << "[ " << vertex[i + 2].red << " " << vertex[i + 2].green << " " << vertex[i + 2].blue << " ] gouraudtriangle\n";
	  }
	} else {
	  // Flat shaded polygon; all vertex colors the same. 
	  output << "newpath\n";
	  output << red << " " << green << " " << blue << " setrgbcolor\n";

	  // Draw a filled triangle. 
	  output << vertex[0].x << " " << vertex[0].y << " moveto\n";
	  for(int i = 1; i < nvertices; i++) {
	    output << vertex[i].x << " " << vertex[i].y << " lineto\n";
	  }
	  output << "closepath fill\n\n";
	}
      }
      loc += nvertices * 7;  // Each vertex element in the
      //   feedback buffer is 7 GLfloats. 
      break;
    case GL_POINT_TOKEN:
      vertex = (Feedback3Dcolor*) loc;
      output << vertex[0].red << " " << vertex[0].green << " " << vertex[0].blue << " setrgbcolor\n";
      output << vertex[0].x << " " << vertex[0].y << " " << pointSize / 2.0 << " 0 360 arc fill\n\n";
      loc += 7;           // Each vertex element in the feedback
      //   buffer is 7 GLfloats.
      break;
    default:
      std::cerr << token << " - opps a GL token I do not understand.\n";
      break;
    }


  }



  // Emit EPS trailer.
  output << "grestore\n\n";
  //output << "showpage\n%% stop using temporary dictionary\nend\n\n%% restore original state\norigstate restore\n\n%%Trailer\n\n";
  //output << "%Add `showpage' to the end of this file to be able to print to a printer.\n";
  output << "showpage\n\n";
  
  output.close();
  std::cerr << "EPS \'" << filename << "\' created containing "<<prims.size()<<" OpenGL primitives, removed "<<count-prims.size()<<".\n";
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
  std::string filename = outputname + std::string(postFix);

  std::ofstream output(filename.c_str(), std::ios::out);
  if (!output){
    std::cerr << "Attempt to open PPM file \'" << filename
	      << "\' for output failed.\n";
    return;
  }

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
  
  char* ppm = new char[size];
  int i=0;
  for(int j=0;j<height;j++){
    for(int k=0;k<width*3;k++){
      int n = (height-1-j)*width*3+k;
      if(pixels[i]<=0.0){
	ppm[n] = 0;
      }
      else if (pixels[i]>=1.0){
	ppm[n] = 255;
      }
      else {
	ppm[n] = (char)(255*pixels[i]);
      }
      i++;
    }
  }
  delete [] pixels;

  output << "P6\n# " << filename << "\n" << width << " " << height << "\n255\n";
  output.write(ppm,size);
  output.close();

  std::cerr << "PPM \'" << filename << "\' created, frame "<<frame<<".\n";

  delete [] ppm;
  currentFrame++;
}

//_____________________________________________________________________ readXYZ
void readXYZ(const std::string& inputname, bool mov, double alpha, double beta, double gamma, int tracer){
  
  std::vector<Vector3D> colors;
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
  r.resize(0);
  std::ifstream infile;
  infile.open(inputname.c_str(), std::ios::in); // open file for reading

  if (!infile){
    std::cerr << "Cannot open \'" << inputname << "\'.\n";
    exit(0);
  }

  std::string tempstr;
  int n;
  totalFrames = 1;
  if(mov){
    infile >> totalFrames; 
    std::getline (infile, tempstr); // Next line
    infile >> n; 
    std::getline (infile, tempstr); // comment
  }
  else {
    infile >> n; // number of atoms first
  
    std::getline (infile, tempstr); // Next line
    std::getline (infile, tempstr); // comment
  }
  double x,y,z;

  colorMap.clear();
  color2index.clear();
  index2color.clear();
  colorVisible.clear();
  r.clear();
  rAll.clear();

  Vector3D a,b;
  for(int j = 0; j < totalFrames; j++){ // now loop through the frames
    for(int i = 0; i < n; i++){ // now loop through the atoms	
      infile >> tempstr;
      if(infile.eof()){
	std::cerr << n << " elements (atomname x y z), read " << i << ", reached end of file.\n";
	i = n;
	j = totalFrames;
      }
      infile >> x;
      infile >> y;
      infile >> z;
      double u,v;
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


      if(j==0&&(i==0||x<a.x)) a.x=x;
      if(j==0&&(i==0||y<a.y)) a.y=y;
      if(j==0&&(i==0||z<a.z)) a.z=z;
      if(j==0&&(i==0||x>b.x)) b.x=x;
      if(j==0&&(i==0||y>b.y)) b.y=y;
      if(j==0&&(i==0||z>b.z)) b.z=z;

      if(tracer > 1 && (i + (n%tracer)) % tracer == 0)
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
    if(mov){
      infile >> tempstr;
      std::getline (infile, tempstr); // Next line
    }
  }

  r.resize(n);
  getFrame(0);
  std::sort(r.begin(),r.end(),cmpVector3D);
  Vector3D next(r[0]);
  
  for(unsigned int i = 0; i < r.size(); i++){
    r[i].x -= next.x;
    r[i].y -= next.y;
    r[i].z -= next.z;
  }
  std::sort(r.begin(),r.end(),cmpVector3D);

  for(unsigned int i = 0; i < 9; i++){
    std::cerr << r[i].x<<","<<r[i].y<<","<<r[i].z<<","<<sqrt(pow(r[i].x-r[0].x,2)+pow(r[i].y-r[0].y,2)+pow(r[i].z-r[0].z,2))<<std::endl;
  }

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
  radius = pow(maxL*maxL*maxL/r.size(),1.0/3.0)/10.0*scale;
  getFrame(0);
}

//_____________________________________________________________________ parseCmdLine
void parseCmdLine(int argc, 
		  char* argv[], 
		  std::string& inputname, 
		  std::string& outputname, 
		  double& scale, 
		  double& size, 
		  int& d, 
		  bool& m, 
		  double& f, 
		  double& t,
		  double& alpha,
		  double& beta,
		  double& gamma,
		  int& tracer,
		  int& startview,
		  double& fov,
		  GLfloat& angp,
		  GLfloat& angp2,
		  bool& bm,
		  double& dfr){
  if (argc < 2 || (argc >= 2 && (std::string(argv[1]) =="-h" ||
                                 std::string(argv[1]) =="--help" ))){
    usage(argv[0]);
    exit(0);
  }

  int cur = 1;
  while (cur<(argc-1) && argv[cur][0]=='-') {

    std::string str(argv[cur]);

    if (str == "-scale") {
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

    if (str == "-movie") {
      m = true;
      cur++;
      continue;
    }

    if (str == "-bmovie") {
      bm = true;
      m = true;
      cur++;
      continue;
    }

    if (str == "-size") {
      size = atof(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-tracer") {
      tracer = atoi(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-subdivisions") {
      d = atoi(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-z"){
      alpha = 0.0;
      beta  = 0.0;
      gamma = 0.0;
      cur++;
      continue;
    }

    if (str == "-y"){
      alpha = 0.0;
      beta  = 0.0;
      gamma = 90.0/180*M_PI;
      cur++;
      continue;
    }

    if (str == "-x"){
      alpha = 0.0;
      beta  = 90.0/180*M_PI;
      gamma = 0.0;
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



    break;

  }

  inputname = std::string(argv[cur++]);
  if (cur < argc){
    outputname = std::string(argv[cur]);
  }
}

//_____________________________________________________________________ usage
void usage(char* name){
  std::cerr << "Usage: " << name << " [-h][--help] [-scale <real>] [-plane <real> <real> [-movie] [-subdivisions <int>] [-size <real>] [-x][-y][-z][-xy][-xz][-yz][-xyz][-a <real>][-b <real>][-c <real>][-tracer <int>] inputname [outputname]\n"
            << " where:\n"
            << "    -h or --help          : This help output\n"
            << "    -scale <real>         : (optional) sphere scaling factor\n"
            << "    -subdivisions <int>   : (optional) number of subdivisions around/anlong the Z axis for sphere\n"
            << "    -size <real>          : (optional) point size and line width\n"
            << "    -movie                : (optional) for XYZ trajectory input\n"
            << "    -plane <real> <real>  : (optional) z-plane projection\n"
            << "    -x                    : (optional) projection axis\n"
            << "    -y                    : (optional) projection axis\n"
            << "    -z                    : (optional) projection axis\n"
            << "    -xy                   : (optional) projection axis\n"
            << "    -xz                   : (optional) projection axis\n"
            << "    -yz                   : (optional) projection axis\n"
            << "    -xyz                  : (optional) projection axis\n"
            << "    -a <real>             : (optional) z-rotation, first\n"
            << "    -b <real>             : (optional) y-rotation, second\n"
            << "    -c <real>             : (optional) x-rotation, third\n"
	    << "    -tracer <int>         : (optional) color every n'th differently\n"
            << "    inputname             : filename of input\n"
            << "    outputname            : (optional) filename of EPS/PPM output\n";
}

//_____________________________________________________________________ display
void display() {
  double fact = dFocusRadius*dFocusRadius/radius/radius;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  glTranslatef(-xT, -yT, -5 * zL);
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
    glutChangeToMenuEntry(2,const_cast<char*>((std::string("Toggle Plane Rotation : ")+(planeRotation?std::string("Yes"):std::string("No"))).c_str()),TOGGLE_PLANE_ROTATION);
    break;
  case QUIT:
    exit(0);
    break;
  case SAVE_FRAMEBUFFER:
    saveFramebufferToFile();
    break;
  case PRINT_BUFFER:
    printBuffer();
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
void bmovieCallback(int v){
  saveFramebufferToFile();
  startAnimation();
}
//_____________________________________________________________________ stopAnimation
void stopAnimation(){
  frameTime = getTime(frameTime);
  if(frameTime > 0 && frameTimeN > 0)
    std::cerr <<"Frames/sec: "<<static_cast<double>(frameTimeN)/frameTime<<".\n";
  timeout = -1;
  glutIdleFunc(0);
  if(bmovie)
    exit(0);
}
//_____________________________________________________________________ startAnimation
void startAnimation(){
  frameTime = getTime();
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
  double t = getTime();
  doAnimation(false);   
  t = getTime(t);
  int n = static_cast<int>(floor(t*1000.0));
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
	std::cerr << "Last frame ("<<frame<<").\n";
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
  std::cerr << "Particle(s)         : " << r.size() << "\n";
  std::cerr << "Frame(s)            : " << totalFrames << "\n";
  std::cerr << "Box                 : [" << xL << "," << yL << "," << zL << "]\n";
  std::cerr << "Types(s)            : " << colorMap.size() << "\n";
  std::cerr << "Sphere radius       : " << radius << "\n";
  std::cerr << "Sphere subdivisions : " << subdivisions << "\n";
  std::cerr << "Scaling             : " << scale << "\n";
  std::cerr << "Opaque              : " << opaque << "\n";
  std::cerr << "Point/line size     : " << lineWidth << "\n";
  std::cerr << "View cube           : " << (viewCube?std::string("on"):std::string("off")) << "\n";
  std::cerr << "Lighting            : " << (lighting?std::string("on"):std::string("off")) << "\n";
  std::cerr << "Switch fill mode    : " << (polygonMode==FILL?std::string("fill"):(polygonMode==LINE?std::string("line"):std::string("point"))) << "\n";
  std::cerr << "Window              : " << winX<<"x"<<winY<<", angles ("<<angle<<", "<<angle2<<") angles plane ("<<anglePlane<<", "<<anglePlane2<<"), start ("<<startx<<", "<<starty<<"), "
	    << "view ("<<startFieldOfView<<", "<<fieldOfView<<"), buffer: "<<bufferSize<<".\n";
  std::cerr << "Picture number      : " << currentFrame << "\n";
  std::cerr << "Animation           : frame: " << frame <<", frame step: " << frameStep <<", frames/sec:";
  if(frameRate > 0)
    std::cerr << 1000.0/static_cast<double>(frameRate);
  else
    std::cerr << "max";
  if(frameTime > 0 && frameTimeN > 0)
    std::cerr <<" showing at frames/sec: "<<static_cast<double>(frameTimeN)/getTime(frameTime)<<".\n";
  std::cerr <<".\n";
  if(from < to) {
    std::cerr << "Projection plane    : " <<from<<" - "<<to<<std::endl;
  }
  else if(dFocusRadius > 0) {
    std::cerr << "Projection focus    : " <<dFocusRadius<<std::endl;
    std::cerr << "Max radius scaling  : " <<sqrt(1.0+dFocusRadius*dFocusRadius/radius/radius*maxL*maxL/4.0)<<std::endl;    
    std::cerr << "Brightness          : " <<brightness()<<std::endl;
  }
  else {
    std::cerr << "3D"<<std::endl;
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
      std::cerr << "Frame: " << frame <<", frame step: " << frameStep <<", frames/sec:";
      if(frameRate > 0)
	std::cerr << 1000.0/static_cast<double>(frameRate);
      else
	std::cerr << "max";
      std::cerr << ".\n";

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
    std::cerr << "Summary of keyboard commands:" << std::endl
	      << "  q, ESC :  Quit." << std::endl
	      << "  w      :  Write PPM file." << std::endl
	      << "  p      :  Write EPS file." << std::endl
	      << "  SPACE  :  Toggle view cube." << std::endl
	      << "  +, -   :  Increase/decrease colors." << std::endl
	      << "  h, ?   :  This message." << std::endl
	      << "  v      :  Debug information." << std::endl;
    if(from < to || dFocusRadius > 0.0)
      std::cerr << "  r      :  Toggle plane rotation." << std::endl;
    if(movie){
      std::cerr << "  s      :  Stop animation." << std::endl
		<< "  i      :  Start animation." << std::endl
		<< "  f      :  Stop animation and one step forward." << std::endl
		<< "  b      :  Stop animation and one step back." << std::endl
		<< "  a      :  Animation statistic." << std::endl
		<< "  TAB    :  Change sign of frameStep." << std::endl
		<< "  c      :  Set actual frame = 0 and frameStep = 1." << std::endl;
    }
    
    break;
  default:
    break;
  }
}

//_____________________________________________________________________ main
int main (int argc, char** argv) {
  glutInit(&argc, argv);

  std::cerr.precision(16);

  std::string inputname;
  double size = 0.0;
  int d = 0;
  double alpha = 0.0;
  double beta = 0.0;
  double gamma = 0.0;
  int tracer = 0;
  parseCmdLine(argc,argv, inputname, outputname, scale, size, d, movie, from, to, alpha, beta, gamma,tracer,startFieldOfView,fieldOfView,anglePlane,anglePlane2,bmovie,dFocusRadius);

  if(size > 0){
    lineWidth = size;
    pointSize = size;
  }
  if(d > 2){
    subdivisions = d;
  }

  if(bmovie || (from < to) || (dFocusRadius > 0.0))
    viewCube = false;    

  planeRotation = (from < to) || (dFocusRadius > 0.0);

  readXYZ(inputname,movie,alpha,beta,gamma,tracer);

  //  if(from < to) {
  //  from -= zT;
  //  to   -= zT;
  // }
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
    glutAddMenuEntry(const_cast<char*>((std::string("Toggle Plane Rotation : ")+(planeRotation?std::string("Yes"):std::string("No"))).c_str()),TOGGLE_PLANE_ROTATION);
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
  

  if(bmovie){
    saveAnimFrame = true;
    glutSetMenu(animation);
    glutChangeToMenuEntry(5,"Save frames : Yes",ANIM_SAVE);
    glutPostRedisplay();
    glutTimerFunc(1000,&bmovieCallback,1);
  }
  
  glutMainLoop();
  return 0;
}
