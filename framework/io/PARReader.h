/*  -*- c++ -*-  */
#ifndef PARREADER_H
#define PARREADER_H

#include "Reader.h"
#include "PAR.h"

namespace ProtoMol {

  //_________________________________________________________________PARReader
  /**
   * Reads PAR-Charmm/XPLOR 19/27 files (ASCII)
   */
  class PARReader : public Reader {
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
    explicit PARReader(PAR::CharmmTypeEnum charmmType=PAR::UNDEFINED);
    explicit PARReader(const std::string& filename, PAR::CharmmTypeEnum charmmType=PAR::UNDEFINED);
    explicit PARReader(const char* filename, PAR::CharmmTypeEnum charmmType=PAR::UNDEFINED);
    // Need this implementation, otherwise const char* will bee converted to bool or int ...

    virtual ~PARReader();

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
    // New methods of class PAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool open(const std::string& filename, PAR::CharmmTypeEnum charmmType);
    bool open(const char* filename, PAR::CharmmTypeEnum charmmType);
    bool open(PAR::CharmmTypeEnum charmmType);

    void setCharmmType(PAR::CharmmTypeEnum charmmType);
    PAR::CharmmTypeEnum getCharmmTypeDetected() const;

    bool read(PAR& par);

    PAR* orphanPAR();

  private:
    static bool isKeywordCharmm28(const std::string& word);
    static bool isKeywordCharmm19(const std::string& word);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class PAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend PARReader& operator>>(PARReader& parReader, PAR& par);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    PAR* myPAR;
    PAR::CharmmTypeEnum myCharmmType;
    PAR::CharmmTypeEnum myCharmmTypeDetected;

  };

  //____________________________________________________________________________INLINES
  inline void PARReader::setCharmmType(PAR::CharmmTypeEnum charmmType){
    myCharmmType = charmmType;
    myCharmmTypeDetected = PAR::UNDEFINED;
  }

  inline PAR::CharmmTypeEnum PARReader::getCharmmTypeDetected() const{
    return myCharmmTypeDetected;
  }

  inline bool PARReader::open(const char* filename, PAR::CharmmTypeEnum charmmType){
    return open(std::string(filename),charmmType);
  }

}
#endif /* PARREADER_H */
