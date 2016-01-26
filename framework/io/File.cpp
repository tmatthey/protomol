#include "File.h"

#include "systemutilities.h"

#if defined(__SUNPRO_CC)
// Sun WorkShop CC does not properly read more than one char ...
#define FSTREAM_READ_1BY1
#endif

using std::string;
using std::fstream;
namespace ProtoMol {
  //_________________________________________________________________ File
  File::File(std::ios::openmode mode):myMode(mode),myFilename(""),myComment(""){}

  File::File(std::ios::openmode mode,const string& filename):myMode(mode),myFilename(filename),myComment(""){}
  
  File::~File(){
    if(myFile.is_open())
      close();
	//std::cout<<"Open file is closed"<<endl;	
  }

  bool File::open(){
    if(myFile.is_open())
      close();
    myFile.clear();
    myFile.open(myFilename.c_str(),myMode);
    
    return !myFile.fail();
  }

  bool File::open(const string& filename){
    setFilename(filename);
    return open();
  }
  
  bool File::open(const char* filename){
    return open(string(filename));
  }
  

  void File::close(){
    if(myFile.is_open())
      myFile.close();
  }

  bool File::isAccessible(){
    return ProtoMol::isAccessible(myFilename);
  }

  File::operator void*() const{
    return (myFile.fail() ? 0 : const_cast<File*>(this));
  }

  bool File::operator!() const { 
    return myFile.fail(); 
  }

  std::fstream& File::read(char* c, std::streamsize count){
#ifdef FSTREAM_READ_1BY1
    for(std::streamsize i=0;i<count;++i)
      myFile.get(c[i]);
#else
    myFile.read(c,count);
#endif
    return myFile;
  }
  std::string File::getline(){
    string res;
    bool ok = (!myFile.fail()) && (!myFile.eof());
    std::getline(myFile,res);
    if(ok && !myFile.bad() && myFile.fail() && myFile.eof())
      myFile.clear(myFile.rdstate() & (~std::ios::failbit));
    return res;
  }


}
