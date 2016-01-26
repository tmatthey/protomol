/*  -*- c++ -*-  */
#ifndef ISGPARREADER_H
#define ISGPARREADER_H

#include "Reader.h"
#include "iSGPAR.h"

namespace ProtoMol {

  //_________________________________________________________________iSGPARReader
  class iSGPARReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    enum PARRecordTypeEnum {
      UNDEFINED,
      BOND,
      ANGLE,
      DIHEDRAL,
      IMPROPER,
      NONBONDED,
      NBFIX,
      HBOND
    };
    // Supported and read Charmm/XPLOR types/records

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    explicit iSGPARReader(iSGPAR::CharmmTypeEnum charmmType=iSGPAR::UNDEFINED);
    explicit iSGPARReader(const std::string& filename, iSGPAR::CharmmTypeEnum charmmType=iSGPAR::UNDEFINED);
    explicit iSGPARReader(const char* filename, iSGPAR::CharmmTypeEnum charmmType=iSGPAR::UNDEFINED);
    // Need this implementation, otherwise const char* will be converted to bool or int ...

    virtual ~iSGPARReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open(){return File::open();};
    virtual bool open(const std::string& filename){return File::open(filename);};
    virtual bool open(const char* filename){return File::open(filename);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGPAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool open(const std::string& filename, iSGPAR::CharmmTypeEnum charmmType);
    bool open(const char* filename, iSGPAR::CharmmTypeEnum charmmType);
    bool open(iSGPAR::CharmmTypeEnum charmmType);

    void setCharmmType(iSGPAR::CharmmTypeEnum charmmType);
    iSGPAR::CharmmTypeEnum getCharmmTypeDetected() const;

    void setNumComps(const int C) {numComp = C;}
    bool read(iSGPAR& par);

    iSGPAR* orphanPAR();

  private:
    static bool isKeywordCharmm28(const std::string& word);
    static bool isKeywordCharmm19(const std::string& word);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // operator that puts the info in the PAR file into the iSGPAR contained
    friend iSGPARReader& operator>>(iSGPARReader& parReader, iSGPAR& par);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    iSGPAR* myPAR;
    iSGPAR::CharmmTypeEnum myCharmmType;
    iSGPAR::CharmmTypeEnum myCharmmTypeDetected;
    int numComp;
  };

  //____________________________________________________________________________INLINES
  inline void iSGPARReader::setCharmmType(iSGPAR::CharmmTypeEnum charmmType){
    myCharmmType = charmmType;
    myCharmmTypeDetected = iSGPAR::UNDEFINED;
  }

  inline iSGPAR::CharmmTypeEnum iSGPARReader::getCharmmTypeDetected() const{
    return myCharmmTypeDetected;
  }

  inline bool iSGPARReader::open(const char* filename, iSGPAR::CharmmTypeEnum charmmType){
    return open(std::string(filename),charmmType);
  }

}
#endif /* ISGPARREADER_H */
