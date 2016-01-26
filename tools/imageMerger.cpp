#include "Report.h"
#include "stringutilities.h"
#include "PPMReader.h"
#include "PGMReader.h"
#include "PGMWriter.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::string;
using std::vector;

//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {


  // parse
  if (argc < 2 || (argc >= 2 && (string(argv[1]) =="-h" ||
                                 string(argv[1]) =="--help" ))){
    report << plain << quit << "usage: "<<argv[0]<<"[-h] [--help] [-auto] [-gamma <real>] [-i <range color from, 0-255> <range color to, 0-255>] -o <output image> <input> <intput> ..."<<endr;
  }


  string out;
  vector<string> in;
  int from = 0;
  int to = 255;
  Real gamma = 1.0;
  unsigned int cur = 1;
  bool autodist = false;

  while (cur<(argc-1) && argv[cur][0]=='-') {

    string str(argv[cur]);

    if (str == "-o") {
      out = argv[++cur];
      cur++;
      continue;
    }
    if (str == "-auto") {
      autodist = true;
      cur++;
      continue;
    }

    if (str == "-gamma") {
      gamma = toReal(argv[++cur]);
      cur++;
      continue;
    }

    if (str == "-i") {
      from = toInt(argv[++cur]);
      to   = toInt(argv[++cur]);
      cur++;
      continue;
    }

    break;
  }
  if(out.empty())
    report << error << "Missing output file."<<endr;

  for(unsigned int i=cur;i<argc;++i){
    if(PPMReader(argv[i]).tryFormat()){
      in.push_back(argv[i]);
    }
    else {
      report << hint <<"\'"<<argv[i]<<"\' neither PPM nor PGM."<<endr;
    }

  }
  if(in.empty())
    report << error << "Missing input file(s)."<<endr;
  

  // read
  long* image = NULL;
  unsigned int count = 0;
  unsigned int w=0,h=0;
  for(unsigned int i=0;i<in.size();++i){
    PGM pgm;
    PPMReader reader(in[i]);
    if(!(reader >> pgm)){
      report << hint <<"Could not read \'"<<in[i]<<"\'."<<endr;      
    }
    else {
      if(image == NULL){
	image = new long [pgm.size()];
	w = pgm.width();
	h = pgm.height();
      }
      else if(w != pgm.width() || h != pgm.height()){
	report << "Image \'"<<in[i]<<"\' "<<pgm.width()<<"x"<<pgm.height()<<" differs from "<<w<<"x"<<h<<"."<<endr;
      }
      count++;
      unsigned char* p1 = pgm.begin();
      long * p2         = image;
      for(unsigned int i=0;i<pgm.size();++i,p1++,p2++){
	(*p2) += (*p1);
      }
    }
  }

  // statistic
  vector<long> hist(256,0);
  int imin = 255;
  int imax = 0;
  for(unsigned int i=0;i<w*h;i++){
    int c = image[i]/count;
    hist[c]++;
    imin = std::min(c,imin);
    imax = std::max(c,imax);
  }
  report << hint << "Image range: ("<<imin <<","<<imax<<")"<<endr;
  if(autodist){
    from = imin;
    to = imax;
  }

  // merge
  PGM pgm(w,h);
  unsigned char* p1 = pgm.begin();
  long * p2         = image;
  imin = 255;
  imax = 0;
  for(unsigned int i=0;i<pgm.size();++i,p1++,p2++){
    int c = pow(((*p2)-from*count)*255.0/(count*(to-from))/255.0,gamma)*255.0+0.5;
    if( c < 0){
      c = 0;
    }
    else if(c > 255){
      c = 255;
    }
    
    (*p1) = 255-c;
    imin = std::min(255-c,imin);
    imax = std::max(255-c,imax);
  }
  report << hint << "Final image range: ("<<imin <<","<<imax<<")"<<endr;

  // write
  PGMWriter writer(out);
  if(!(writer << pgm))
    report << error << "Could not write \'"<<out<<"\'."<<endr;      


  // clean up
  if(image != NULL)
    delete [] image;

  return 0;
}
