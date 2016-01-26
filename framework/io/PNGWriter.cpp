#include "protomol.h"
#ifdef HAVE_LIBPNG
#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#endif

#include "PNGWriter.h"

#include "systemutilities.h"
#include "Report.h"
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________PNGWriter

  PNGWriter::PNGWriter():Writer(){
#ifndef HAVE_LIBPNG
    report << recoverable <<"Compiled without PNG library."<<endr;    
#endif
  }

  PNGWriter::PNGWriter(const std::string& filename):Writer(filename){
#ifndef HAVE_LIBPNG
    report << recoverable <<"Compiled without PNG library."<<endr;    
#endif
  }

  bool PNGWriter::write(const PPM& ppm){

#ifndef HAVE_LIBPNG
    ppm.size();
    unsigned char png[] ={
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d, 0x49, 0x48, 0x44, 0x52, 
      0x00, 0x00, 0x00, 0xd1, 0x00, 0x00, 0x00, 0x0a, 0x08, 0x02, 0x00, 0x00, 0x00, 0xa5, 0xa1, 0x2c, 
      0x54, 0x00, 0x00, 0x00, 0x06, 0x62, 0x4b, 0x47, 0x44, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0xa0, 
      0xbd, 0xa7, 0x93, 0x00, 0x00, 0x00, 0x09, 0x70, 0x48, 0x59, 0x73, 0x00, 0x00, 0x0b, 0x12, 0x00, 
      0x00, 0x0b, 0x12, 0x01, 0xd2, 0xdd, 0x7e, 0xfc, 0x00, 0x00, 0x00, 0x07, 0x74, 0x49, 0x4d, 0x45, 
      0x07, 0xd5, 0x01, 0x07, 0x11, 0x19, 0x24, 0x94, 0x77, 0xb2, 0xbb, 0x00, 0x00, 0x01, 0x4b, 0x49, 
      0x44, 0x41, 0x54, 0x78, 0x9c, 0xed, 0x57, 0x5d, 0x0f, 0xc3, 0x20, 0x08, 0xd4, 0x65, 0xff, 0xff, 
      0x2f, 0x77, 0x0f, 0x26, 0x86, 0xc1, 0x71, 0xc2, 0xea, 0x3a, 0xb2, 0xf4, 0x1e, 0x16, 0x56, 0x29, 
      0x9c, 0xf2, 0x21, 0xed, 0xc7, 0x71, 0xb4, 0x77, 0xf4, 0xde, 0x87, 0x60, 0x97, 0xb6, 0xa3, 0xf7, 
      0x9e, 0xf2, 0x92, 0xd5, 0x4f, 0x99, 0x3a, 0x6f, 0x3c, 0x62, 0x61, 0x1e, 0xef, 0xc0, 0xd0, 0x1f, 
      0x0f, 0xa5, 0xec, 0xfd, 0x2d, 0x0b, 0x6f, 0xef, 0xf6, 0xf9, 0x03, 0x6a, 0x0c, 0xa8, 0xd3, 0xb9, 
      0xb1, 0x05, 0x23, 0x00, 0xe3, 0x84, 0x1b, 0xca, 0xb6, 0xa9, 0xd0, 0xfe, 0x34, 0x1c, 0xcf, 0xa0, 
      0x1e, 0x29, 0xbe, 0x99, 0xc8, 0x43, 0xb0, 0x75, 0x99, 0xaa, 0x54, 0xa8, 0x1c, 0x39, 0x6e, 0x45, 
      0x09, 0xca, 0x4a, 0x99, 0x2c, 0x41, 0xef, 0x3c, 0x4b, 0x88, 0x4d, 0x0f, 0x7b, 0x33, 0x29, 0xc2, 
      0x13, 0xc6, 0xc8, 0x0b, 0x9c, 0xb5, 0xe9, 0x29, 0xa7, 0xf6, 0xce, 0x72, 0x4e, 0xb9, 0x91, 0xb2, 
      0x64, 0x3f, 0x7f, 0x65, 0xfe, 0xc9, 0xb7, 0x54, 0x9a, 0x12, 0xc0, 0xd7, 0xd5, 0x43, 0xc2, 0x56, 
      0x7a, 0xe4, 0xde, 0x15, 0x55, 0x4b, 0xc0, 0xf3, 0xae, 0x5a, 0x8e, 0x95, 0xb3, 0xa3, 0x42, 0x84, 
      0xd5, 0x5c, 0xe5, 0xa6, 0x22, 0x3c, 0x61, 0x8c, 0xbc, 0xc0, 0x79, 0xd1, 0xb4, 0x24, 0x53, 0x7b, 
      0x8f, 0xf6, 0xb9, 0xcb, 0x70, 0xbe, 0xee, 0xcf, 0x8c, 0x65, 0x97, 0xcd, 0x4c, 0xd9, 0x29, 0xcd, 
      0x76, 0x5f, 0x05, 0xaf, 0x51, 0x15, 0x04, 0xcb, 0xb9, 0x8d, 0x03, 0x7b, 0x1c, 0xc5, 0xcf, 0x6b, 
      0x17, 0x48, 0x33, 0x83, 0x55, 0x17, 0xb9, 0x85, 0xe1, 0xbd, 0x54, 0x10, 0xfa, 0x1b, 0xa2, 0x0e, 
      0x3e, 0x68, 0x78, 0xea, 0xae, 0xff, 0xaa, 0xaf, 0x9f, 0xd8, 0xac, 0xe3, 0x8e, 0x43, 0x92, 0x01, 
      0x03, 0x0c, 0x9c, 0x69, 0xac, 0xb6, 0x37, 0x9f, 0xca, 0x89, 0x52, 0xd5, 0xa2, 0x37, 0x84, 0x79, 
      0xe3, 0x9d, 0xfd, 0x3a, 0x81, 0xde, 0x79, 0x3e, 0xd9, 0x9c, 0x83, 0xde, 0xa1, 0x59, 0x28, 0x67, 
      0x29, 0x2d, 0x2f, 0x38, 0x7e, 0x44, 0x64, 0x84, 0x5d, 0x56, 0x51, 0x90, 0xe7, 0x32, 0x46, 0xcb, 
      0x6f, 0x08, 0xa8, 0x49, 0x08, 0x68, 0xb5, 0xbd, 0x4d, 0xb8, 0x78, 0x57, 0xbf, 0xd1, 0x9c, 0x18, 
      0x5d, 0x19, 0xb8, 0xba, 0x77, 0xeb, 0x8d, 0x7f, 0xc5, 0x0b, 0xb2, 0xa6, 0xd6, 0xf6, 0x66, 0x30, 
      0x05, 0x29, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82    
    };
    myFile.write(reinterpret_cast<char*>(png),sizeof(png));
#endif

    close();



#ifdef HAVE_LIBPNG
    FILE *fp= fopen(myFilename.c_str(), "wb");
    if (!fp){
      myFile.setstate(std::ios::failbit);
      return false;
    }
    // Create and initialize the png_struct with the default error handlers 
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
      report << error <<"Could not initialize PNG library."<<endr;
      fclose(fp);
      myFile.setstate(std::ios::failbit);
      return false; 
    }
 
    // Allocate/initialize the memory for image information.  REQUIRED.
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      report << error <<"Could not allocate/initialize the memory for PNG image information."<<endr;
      fclose(fp);
      myFile.setstate(std::ios::failbit);     
      return false; 
    }
 
    // Set error handling for setjmp/longjmp method of libpng error handling 
    if (setjmp(png_jmpbuf(png_ptr))) {
      // Free all of the memory associated with the png_ptr and info_ptr 
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      // If we get here, we had a problem writing the file 
      report << error <<"Could not open PNG image for writting."<<endr;
      fclose(fp);
      myFile.setstate(std::ios::failbit);     
      return false; 
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
    return true;
#else
    report << recoverable <<"Compiled without PNG library."<<endr;    
    return false;
#endif
  }

  PNGWriter& operator<<(PNGWriter& PNGWriter, const PPM& ppm){
    PNGWriter.write(ppm);
    return PNGWriter;
  }

  bool PNGWriter::write(const PGM& pgm){
    PPM ppm;
    ppm = pgm;
    return write(ppm);   
  }

  PNGWriter& operator<<(PNGWriter& pngWriter, const PGM& pgm){
    pngWriter.write(pgm);
    return pngWriter;
  }


}