/*  -*- c++ -*-  */
#ifndef PPM_H
#define PPM_H

namespace ProtoMol {
  class PGM;

  //_____________________________________________________________________ PPM
  /**
    Container for PPM binary image.@n @n

    PPM memory layout:@n @n

    h                @n
    |  1. ----->     @n
    |  2. ----->     @n
    |  ...           @n
    0  N. ----->     @n
    y                @n
       x  0-----w    @n
  */
  class PPM{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ~PPM();
    PPM();
    PPM(int x, int y);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PPM
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// set point at (x,y) with RGB (a,a,a)
    void set(unsigned int x, unsigned int y,unsigned char a){p[(h-y-1)*w*3+x*3] = a;p[(h-y-1)*w*3+x*3+1] = a;p[(h-y-1)*w*3+x*3+2] = a;}
    /// set point at (x,y) with RGB (r,g,b)
    void set(unsigned int x, unsigned int y,unsigned char r,unsigned char g,unsigned char b){p[(h-y-1)*w*3+x*3] = r;p[(h-y-1)*w*3+x*3+1] = g;p[(h-y-1)*w*3+x*3+2] = b;}
    /// get gray value, (r+g+b)/3
    unsigned char get(unsigned int x, unsigned int y) const {return (p[(h-y-1)*w*3+x*3]+p[(h-y-1)*w*3+x*3+1]+p[(h-y-1)*w*3+x*3+2])/3;}
    unsigned char getRed(unsigned int x, unsigned int y) const {return p[(h-y-1)*w*3+x*3];}
    unsigned char getGreen(unsigned int x, unsigned int y) const {return p[(h-y-1)*w*3+x*3+1];}
    unsigned char getBlue(unsigned int x, unsigned int y) const {return p[(h-y-1)*w*3+x*3+2];}

    unsigned int width() const{return w;}
    unsigned int height() const{return h;}
    unsigned int size() const{return 3*w*h;}
    void resize(unsigned int width, unsigned int height);
    unsigned char* begin() const{return &p[0];}
    unsigned char* end() const{return &p[3*w*h];}

    void clear();

    PPM& operator=(const PGM& pgm);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    unsigned int w,h;
    unsigned char* p;

  };
}
#endif /* PPM_H */
