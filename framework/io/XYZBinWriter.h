/*  -*- c++ -*-  */
#ifndef XYZBINWRITER_H
#define XYZBINWRITER_H

#include "Writer.h"
#include "XYZ.h"
#include "systemutilities.h"

namespace ProtoMol {

  //_________________________________________________________________XYZBinWriter
  /**
   * Writes a XYZ binary file with 4 or 8 byte floats and litle or big endian.@n
   *
   * Format:
   * - 4 byte int    : N, number of coordinates / Vector3D
   * - N*(3 float's) : N * (x, y, z) coordinates/ Vector3D
   */
  class XYZBinWriter : public Writer {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    explicit XYZBinWriter(bool isLittleEndian=ISLITTLEENDIAN, unsigned int size=sizeof(Real));
    explicit XYZBinWriter(const std::string& filename, bool isLittleEndian=ISLITTLEENDIAN, unsigned int size=sizeof(Real));
    explicit XYZBinWriter(const char* filename, bool isLittleEndian=ISLITTLEENDIAN, unsigned int size=sizeof(Real));
    // Need this implementation, otherwise const char* will bee converted to bool or int ...

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open(){return File::open();};
    virtual bool open(const std::string& filename){return File::open(filename);};
    virtual bool open(const char* filename){return File::open(filename);}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZ
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool open(bool isLittleEndian, unsigned int size=sizeof(Real));
    bool open(const std::string& filename, bool isLittleEndian, unsigned int size=sizeof(Real));
    bool open(const char* filename, bool isLittleEndian, unsigned int size=sizeof(Real));

    bool write(const XYZ& xyz);
    bool write(const Vector3DBlock& coords);
    bool write(const XYZ& xyz, bool isLittleEndian, unsigned int size=sizeof(Real));
    bool write(const Vector3DBlock& coords, bool isLittleEndian, unsigned int size=sizeof(Real));

    void setLittleEndian(bool littleEndian);
    void setSize(unsigned int size);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend XYZBinWriter& operator<<(XYZBinWriter& xyzWriter, const XYZ& xyz);
    friend XYZBinWriter& operator<<(XYZBinWriter& xyzWriter, const Vector3DBlock& coords);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    bool myIsLittleEndian;
    unsigned int mySize;
  };

  //____________________________________________________________________________INLINES
}
#endif /* XYZBINWRITER_H */
