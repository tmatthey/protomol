#include "XYZBinWriter.h"

#include <iomanip>

#include "Report.h"
#include "stringutilities.h"
#include "systemutilities.h"
#include "typeSelection.h"

using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________XYZBinWriter

  XYZBinWriter::XYZBinWriter(bool isLittleEndian, unsigned int size):
    Writer(std::ios::binary|std::ios::trunc),
    myIsLittleEndian(isLittleEndian),mySize(size)
  {
    setSize(size);
  }

  XYZBinWriter::XYZBinWriter(const std::string& filename, bool isLittleEndian, unsigned int size):
    Writer(std::ios::binary|std::ios::trunc,filename),
    myIsLittleEndian(isLittleEndian),mySize(size)
  {
    setSize(size);
  }

  XYZBinWriter::XYZBinWriter(const char* filename, bool isLittleEndian, unsigned int size):
    Writer(std::ios::binary|std::ios::trunc,string(filename)),
    myIsLittleEndian(isLittleEndian),mySize(size)
  {
    setSize(size);
  }

  void XYZBinWriter::setLittleEndian(bool littleEndian){
    myIsLittleEndian = littleEndian;
  }

  bool XYZBinWriter::open(const char* filename, bool isLittleEndian, unsigned int size){
    return open(std::string(filename),isLittleEndian,size);
  }


  bool XYZBinWriter::open(bool isLittleEndian, unsigned int size){
    setLittleEndian(isLittleEndian);
    setSize(size);
    return open();
  }

  bool XYZBinWriter::open(const std::string& filename, bool isLittleEndian, unsigned int size){
    setLittleEndian(isLittleEndian);
    setSize(size);
    return open(filename);
  }

  bool XYZBinWriter::write(const XYZ& xyz, bool isLittleEndian, unsigned int size){
    setLittleEndian(isLittleEndian);
    setSize(size);
    return write(xyz.coords);
  }

  bool XYZBinWriter::write(const Vector3DBlock& coords, bool isLittleEndian, unsigned int size){
    setLittleEndian(isLittleEndian);
    setSize(size);
    return write(coords);
  }

  bool XYZBinWriter::write(const XYZ& xyz){
    return write(xyz.coords);
  }

  bool XYZBinWriter::write(const Vector3DBlock& coords){
    if (!open())
      return false;

    bool swapEndian = (myIsLittleEndian != ISLITTLEENDIAN);
    const unsigned int count = coords.size();
    typedef TypeSelection::Int<4>::type int32;
    int32 n = static_cast<int32>(count);
    if(swapEndian)
      swapBytes(n);
    myFile.write(reinterpret_cast<char*>(&n), 4);

    if(swapEndian)
      report << hint <<"[XYZBinWriter::write] Writing "<<(ISLITTLEENDIAN?"big":"little")
	       <<"endian output on "<<(ISLITTLEENDIAN?"little":"big")<<"endian machine."<<endr;

    if(mySize == sizeof(Real)){
      Real* vec = new Real[count*3];
      for(unsigned int i=0;i<count;++i){
	vec[i*3+0] = static_cast<Real>(coords[i].x);
	vec[i*3+1] = static_cast<Real>(coords[i].y);
	vec[i*3+2] = static_cast<Real>(coords[i].z);
	if(swapEndian){
	  swapBytes(vec[i*3+0]);
	  swapBytes(vec[i*3+1]);
	  swapBytes(vec[i*3+2]);
	}
      }
      myFile.write(reinterpret_cast<char*>(vec),count*3*sizeof(Real));    
      
      delete [] vec;
    }
    else if(mySize == sizeof(float)){
      float* vec = new float[count*3];
      for(unsigned int i=0;i<count;++i){
	vec[i*3+0] = static_cast<float>(coords[i].x);
	vec[i*3+1] = static_cast<float>(coords[i].y);
	vec[i*3+2] = static_cast<float>(coords[i].z);
	if(swapEndian){
	  swapBytes(vec[i*3+0]);
	  swapBytes(vec[i*3+1]);
	  swapBytes(vec[i*3+2]);
	}
      }
      myFile.write(reinterpret_cast<char*>(vec),count*3*sizeof(float));    
      
      delete [] vec;
      report << hint << "[XYZBinWriter::write] Conversion from "<<sizeof(Real)<<" to "<<sizeof(float)<<" bytes."<<endr;
    }
    else if(mySize == sizeof(double)){
      double* vec = new double[count*3];
      for(unsigned int i=0;i<count;++i){
	vec[i*3+0] = static_cast<double>(coords[i].x);
	vec[i*3+1] = static_cast<double>(coords[i].y);
	vec[i*3+2] = static_cast<double>(coords[i].z);
	if(swapEndian){
	  swapBytes(vec[i*3+0]);
	  swapBytes(vec[i*3+1]);
	  swapBytes(vec[i*3+2]);
	}
      }
      myFile.write(reinterpret_cast<char*>(vec),count*3*sizeof(double));    
      
      delete [] vec;
      report << hint << "[XYZBinWriter::write] Conversion from "<<sizeof(Real)<<" to "<<sizeof(double)<<" bytes."<<endr;
    }
    else {
      report << error << "[XYZBinWriter::write]"
	     << " XYZBin file \'" << myFilename
	     << "\' could find adequate float nor double type."<<endr;    
    }
    
    close();
    return !myFile.fail();
  }

  void XYZBinWriter::setSize(unsigned int size){
    if(size != sizeof(Real) && size != sizeof(float) && size != sizeof(double))
      report << error << "[XYZBinWriter::setSize] Size "<<size<<" not supported."<<endr;
    mySize = size;
  }

  XYZBinWriter& operator<<(XYZBinWriter& xyzWriter, const XYZ& xyz){
    xyzWriter.write(xyz.coords);
    return xyzWriter;
  }

  XYZBinWriter& operator<<(XYZBinWriter& xyzWriter, const Vector3DBlock& coords){
    xyzWriter.write(coords);
    return xyzWriter;
  }


}
