/*  -*- c++ -*-  */
#ifndef STAGEREADER_H
#define STAGEREADER_H

#include "Reader.h"
#include "STAGE.h"
#include "Array.h"

namespace ProtoMol {

  //_________________________________________________________________STAGEReader
  class STAGEReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    enum STAGERecordTypeEnum {
      UNDEFINED,
      OSMOTICS,
      ATOMTYPE,
      ATOMCHARGE
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    explicit STAGEReader();
    explicit STAGEReader(const std::string& filename);
    explicit STAGEReader(const char* filename);
    // Need this implementation, otherwise const char* will be converted to bool or int ...

    virtual ~STAGEReader();

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
    // New methods of class STAGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(STAGE& stage);

    STAGE* orphanSTAGE();

    // stream operator
    friend STAGEReader& operator>>(STAGEReader& stageReader, STAGE& stage);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    STAGE* mySTAGE;

    //  the # of osmotic components
    unsigned int numOsmotics;

    // the structure and cooridinate filenames for the different osmotic molecule types in the system
    std::vector<std::string> PSFnames;
    std::vector<std::string> Coordnames;

    // flag so we know if any alphaLJ parameters have been specified
    bool got_alphaLJ;
  };

  //____________________________________________________________________________INLINES

}
#endif /* STAGEREADER_H */
