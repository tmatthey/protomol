/*  -*- c++ -*-  */
#ifndef TRANSREADER_H
#define TRANSREADER_H

#include "Reader.h"
#include "TRANS.h"
#include "Array.h"

namespace ProtoMol {

  //_________________________________________________________________TRANSReader
  class TRANSReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    enum TRANSRecordTypeEnum {
      UNDEFINED,
      IDENTITIES,
      STAGES,
      IDEAL_GAS_DELTAMU,
      ATOMTYPE,
      ATOMCHARGE
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    explicit TRANSReader();
    explicit TRANSReader(const std::string& filename);
    explicit TRANSReader(const char* filename);
    // Need this implementation, otherwise const char* will be converted to bool or int ...

    virtual ~TRANSReader();

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
    // New methods of class TRANS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(TRANS& trans);

    TRANS* orphanTRANS();

    // stream operator
    friend TRANSReader& operator>>(TRANSReader& transReader, TRANS& trans);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TRANS* myTRANS;
    
    // the number of different identities for each atom type
    int numIdentities;

    // the number of stages to break the molecular transformation into
    int numStages;

    // the number of different atom types in the system
    int numTypes;

    // array containing the different atom type names
    std::vector<std::string> myTypes;

    // flag so we know if the ideal gas chemical potential difference
    // matrix has been sized
    bool sized;

    // flag so we know if any alphaLJ parameters have been specified
    bool got_alphaLJ;
  };

  //____________________________________________________________________________INLINES

}
#endif /* TRANSREADER_H */
