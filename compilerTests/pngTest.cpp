#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#include <iostream>

/*
  icc -Wall -wd810,383,981,279,444 -O3 -ip -static -DNDEBUG  pngTest.cpp -o pngTest -lpng -lz;./pngTest ;xv pngTest.png
*/

class PPM{
public:
  ~PPM();
  PPM();
  PPM(int a, int b);
  void set(unsigned int i, unsigned int j,unsigned char a){p[(h-i-1)*w*3+j*3] = a;p[(h-i-1)*w*3+j*3+1] = a;p[(h-i-1)*w*3+j*3+2] = a;};
  void set(unsigned int i, unsigned int j,unsigned char r,unsigned char g,unsigned char b){p[(h-i-1)*w*3+j*3] = r;p[(h-i-1)*w*3+j*3+1] = g;p[(h-i-1)*w*3+j*3+2] = b;};
  unsigned char get(unsigned int i, unsigned int j) const {return (p[(h-i-1)*w*3+j*3]+p[(h-i-1)*w*3+j*3+1]+p[(h-i-1)*w*3+j*3+2])/3;}

  unsigned int width() const{return w;}
  unsigned int height() const{return h;}
  unsigned int size() const{return 3*w*h;}
  void resize(unsigned int width, unsigned int height);
  unsigned char* begin() const{return &p[0];}
  unsigned char* end() const{return &p[3*w*h];}

  void clear();
private:
  unsigned int w,h;
  unsigned char* p;

};

PPM::~PPM(){delete [] p;}

PPM::PPM():w(0),h(0),p(NULL){}

PPM::PPM(int a, int b):w(a),h(b),p(new unsigned char[3*a*b]){
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

//_____________________________________________________________________

int main(int argc, char **argv) {

  PPM ppm(256,256);
  int l = std::min(ppm.height(),ppm.width());
  for(int i=0;i<l;++i){
    ppm.set(i,i,i%256,(i+80)%256,(i+160)%256);
  }

  FILE *fp= fopen("pngTest.png", "wb");
  if (!fp){
    exit(-1);
  }
  // Create and initialize the png_struct with the default error handlers 
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    std::cerr <<"Could not initialize PNG library."<<std::endl;
    fclose(fp);
    exit(-1);
  }
 
  // Allocate/initialize the memory for image information.  REQUIRED.
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    std::cerr <<"Could not allocate/initialize the memory for PNG image information."<<std::endl;
    fclose(fp);
    exit(-1);
  }
 
  // Set error handling for setjmp/longjmp method of libpng error handling 
  if (setjmp(png_jmpbuf(png_ptr))) {
    // Free all of the memory associated with the png_ptr and info_ptr 
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    // If we get here, we had a problem writing the file 
    std::cerr <<"Could not open PNG image for writting."<<std::endl;
    fclose(fp);
    fclose(fp);
    exit(-1);
  }
 
  // Set up the input control if you are using standard C streams 
  png_init_io(png_ptr, fp);
 
  png_set_IHDR(png_ptr, info_ptr, ppm.width(), ppm.height(),
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
 
  png_set_gAMA(png_ptr, info_ptr, 1.0);
 
  png_textp text_ptr = (png_textp) png_malloc(png_ptr, (png_uint_32)sizeof(png_text) * 2);
 
  text_ptr[0].key = "Description";
  text_ptr[0].text = "ProtoMol PNG Writer";
  text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
#ifdef PNG_iTXt_SUPPORTED
  text_ptr[0].lang = NULL;
#endif
 
  text_ptr[1].key = "Software";
  text_ptr[1].text = "ProtoMol";
  text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
#ifdef PNG_iTXt_SUPPORTED
  text_ptr[1].lang = NULL;
#endif
 
  png_bytep* row_pointers = (png_bytep *) png_malloc(png_ptr, ppm.height()*sizeof(png_bytep));
  unsigned char* img = ppm.begin();
  for (unsigned int y=0; y<ppm.height(); y++) {
    row_pointers[y] = &img[y * ppm.width() * 3];
  }
 
  png_set_rows(png_ptr, info_ptr, row_pointers);
 
  // one-shot call to write the whole PNG file into memory 
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
 
  png_free(png_ptr, row_pointers);
  png_free(png_ptr, text_ptr);
 
  // clean up after the write and free any memory allocated - REQUIRED 
  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);

  fclose(fp);


  
  return 0;
}

