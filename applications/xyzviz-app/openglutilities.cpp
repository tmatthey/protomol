#include "openglutilities.h"
#include "mathutilities.h"
#include "stringutilities.h"
#include "pmconstants.h"
#include "Report.h"
#include <vector>
#include <algorithm>

#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)
#include <GL/gl.h>
#endif

using std::ostream;
using std::vector;
using std::string;
using std::sort;
using namespace ProtoMol;
using namespace ProtoMol::Report;

namespace ProtoMol {

#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)

  static GLint bufferSize = 0;

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
    bool operator<(const DepthIndex& a2) const{
      return (a2.depth < depth);
    }    
  };

  //_____________________________________________________________________ print3DcolorVertex
  static string print3DcolorVertex(GLint size, GLint * count, GLfloat * buffer){
    string res = " ";
    for (unsigned int i = 0; i < 7; i++) {
      res += " "+toString(buffer[size - (*count)],4,2);
      *count = *count - 1;
    }
    return res;
  }
#endif

  //_____________________________________________________________________ openglToEPS
  unsigned int openglToPlain(ostream& output,void (*display)()){
#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)
    vector<GLfloat> feedbackBuffer(bufferSize);
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

    int token, nvertices;
    GLint count = size;
    while (count) {
      token = static_cast<int>(buffer[size - count]);
      count--;
      switch (token) {
      case GL_PASS_THROUGH_TOKEN:
	output << "GL_PASS_THROUGH_TOKEN\n";
	output << "  "<<toString(buffer[size - count],4,2) <<"\n";
	count--;
	break;
      case GL_POINT_TOKEN:
	output << "GL_POINT_TOKEN\n";
	output << print3DcolorVertex(size, &count, buffer) <<"\n";
	break;
      case GL_LINE_TOKEN:
	output << "GL_LINE_TOKEN\n";
	output << print3DcolorVertex(size, &count, buffer) <<"\n";
	output << print3DcolorVertex(size, &count, buffer) <<"\n";
	break;
      case GL_LINE_RESET_TOKEN:
	output << "GL_LINE_RESET_TOKEN\n";
	output << print3DcolorVertex(size, &count, buffer) <<"\n";
	output << print3DcolorVertex(size, &count, buffer) <<"\n";
	break;
      case GL_POLYGON_TOKEN:
	output << "GL_POLYGON_TOKEN\n";
	nvertices = static_cast<int>(buffer[size - count]);
	count--;
	for (; nvertices > 0; nvertices--) {
	  output << print3DcolorVertex(size, &count, buffer) <<"\n";
	}
      default:
	break;
      }
    }
    count = size;
#else
    unsigned int count = 0;

#endif
    return count;
  }

  //_____________________________________________________________________ openglToEPS
  unsigned int openglToEPS(ostream& output,void (*display)()){
    unsigned int count = 0;
#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)

    //
    vector<GLfloat> feedbackBuffer(bufferSize);
    while(true){
      glFeedbackBuffer(feedbackBuffer.size(), GL_3D_COLOR, &(feedbackBuffer[0]));
      glRenderMode(GL_FEEDBACK);
      (*display)();
      GLint returned = glRenderMode(GL_RENDER);
      if(returned > 0){
	feedbackBuffer.resize(returned);
	break;
      }
      bufferSize += 100*1024;
      feedbackBuffer.resize(bufferSize);
    }


    GLfloat clearColor[4], viewport[4];
    GLfloat lineWidth;
    GLfloat pointSize;

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
    output << "%%Creator: ProtoMol (using OpenGL feedback)\n";
    output << "%%BoundingBox: " << viewport[0] << " " << viewport[1] << " " << viewport[2] << " " << viewport[3] << "\n";
    output << "%%EndComments\n";
    output << "\n";
    output << "gsave\n";
    output << "\n";

    // Output Frederic Delhoume's "gouraudtriangle" PostScript fragment. 
    output << "% the gouraudtriangle PostScript fragement below is free\n";
    output << "% written by Frederic Delhoume (delhoume@ilog.fr)\n";
    output << "/threshold " << Constant::EPS_GOURAUD_THRESHOLD << " def\n"
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
    vector<DepthIndex> prims;
    GLfloat* loc = &(feedbackBuffer[0]);
    GLfloat* end = loc + feedbackBuffer.size();
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
	report << recoverable << token << " - opps a GL token I do not understand."<<endr;

	break;
      }
    }

    // Sort them
    sort(prims.begin(),prims.end());

    // Print them
    for(unsigned int j=0;j<prims.size();j++){
      GLfloat* loc = prims[j].ptr;

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
	  steps = static_cast<int>(std::max(1.0, static_cast<double>(colormax * distance) * Constant::EPS_SMOOTH_LINE_FACTOR));

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
	    //int triOffset;
	  
	    // Break polygon into "nvertices-2" triangle fans. 
	    for(int i = 0; i < nvertices - 2; i++) {
	      //triOffset = i * 7;
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
	report << recoverable << token << " - opps a GL token I do not understand"<<endr;
	break;
      }


    }



    // Emit EPS trailer.
    output << "grestore\n\n";
    //output << "showpage\n%% stop using temporary dictionary\nend\n\n%% restore original state\norigstate restore\n\n%%Trailer\n\n";
    //output << "%Add `showpage' to the end of this file to be able to print to a printer.\n";
    output << "showpage\n\n";    
#else
    output << "%!PS-Adobe-2.0 EPSF-2.0\n";
    output << "%%Title: nosupport.eps\n";
    output << "%%Creator: ProtoMol\n";
    output << "%%CreationDate: Jan  1 00:00:00 1970\n";
    output << "%%For: ProtoMol\n";
    output << "%%BoundingBox: 0 0 165 12\n";
    output << "%%Magnification: 1.0000\n";
    output << "%%EndComments\n";
    output << "/$F2psDict 200 dict def\n";
    output << "$F2psDict begin\n";
    output << "$F2psDict /mtrx matrix put\n";
    output << "end\n";
    output << "save\n";
    output << "newpath 0 12 moveto 0 0 lineto 165 0 lineto 165 12 lineto closepath clip newpath\n";
    output << "-63.0 88.7 translate\n";
    output << "1 -1 scale\n";
    output << "10 setmiterlimit\n";
    output << " 0.06000 0.06000 scale\n";
    output << "/Times-Roman findfont 180.00 scalefont setfont\n";
    output << "1050 1425 moveto\n";
    output << "gsave 1 -1 scale (Compiled without OpenGL support!) 0.000 0.000 0.000 setrgbcolor show grestore\n";
    output << "restore\n";
#endif
    return count;

  }

  //_____________________________________________________________________ openglToPPM
  void openglToPPM(PPM& ppm){

#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)
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

    ppm.resize(width,height);

    GLfloat* p = pixels;
    for(int j=0;j<height;j++){
      for(int i=0;i<width;i++){
	unsigned char c[3];
	for(unsigned int k=0;k<3;k++){
	  if((*p) <= 0.0){
	    c[k] = 0U;
	  }
	  else 	if((*p) >= 1.0){
	    c[k] = 255U;
	  }
	  else {
	    c[k] = static_cast<unsigned char>((*p)*255);
	  }
	  p++;
	}
	ppm.set(i,j,c[0],c[1],c[2]);
      }
    }
    delete [] pixels;
#else
    ppm.resize(0,0);
#endif


  }

  //_____________________________________________________________________ openglToPGM
  void openglToPGM(PGM& pgm){

#if defined(HAVE_OPENGL) || defined(HAVE_GLUT)
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

    pgm.resize(width,height);

    GLfloat* p = pixels;
    for(int j=0;j<height;j++){
      for(int i=0;i<width;i++){
	int c = 0;
	for(unsigned int k=0;k<3;k++){
	  if((*p) >= 1.0){
	    c += 255;
	  }
	  else 	if((*p) > 0.0){
	    c += (*p)*255;
	  }
	  p++;
	}
	pgm.set(i,j,static_cast<unsigned char>(c/3));
      }
    }
    delete [] pixels;
#else
    pgm.resize(0,0);
#endif


  }


}
