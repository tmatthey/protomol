#include "PPM.h"
#include "PGM.h"
#include <vector> // NULL definition

namespace ProtoMol {

  //_____________________________________________________________________ PPM
  PPM::~PPM(){delete [] p;}

  PPM::PPM():w(0),h(0),p(NULL){
  }

  PPM::PPM(int x, int y):w(x),h(y),p(new unsigned char[3*x*y]){
    clear();
  }
  
  void PPM::clear(){
    for(unsigned int i=0;i<3*w*h;i++) 
      p[i]= static_cast<unsigned char>(0);
  }

  void PPM::resize(unsigned int width, unsigned int height){
    if(p != NULL && width*height*3 != w*h*3){
      delete [] p;
      p = NULL;
    }
    if(p == NULL && width*height*3>0)
      p = new unsigned char[width*height*3];
    w = width;
    h = height;
    clear();
  }

  PPM& PPM::operator=(const PGM& pgm){
    resize(pgm.width(),pgm.height());
    unsigned char* p0 = begin();
    for(unsigned char* p1=pgm.begin();p1 != pgm.end();p1++){
      *p0++ = *p1;
      *p0++ = *p1;
      *p0++ = *p1;
    }
    return (*this);
  }

}
