/*  -*- c++ -*-  */
#ifndef EIGENVECTOTTEXTREADER_H
#define EIGENVECTOTTEXTREADER_H

#include "Reader.h"
#include "EigenvectorInfo.h"
#include "typeSelection.h"

namespace ProtoMol {

  //_________________________________________________________________PDBReader
  /**
   * Reads an EigenvectorInfo ASCII file.
   */
  class EigenvectorTextReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    EigenvectorTextReader();
    explicit EigenvectorTextReader(const std::string& filename);    
    explicit EigenvectorTextReader(const char* filename);
    EigenvectorTextReader(std::ios::openmode mode);
    virtual ~EigenvectorTextReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open(){return File::open();}
    virtual bool open(const std::string& filename){return File::open(filename);}
    virtual bool open(const char* filename){return File::open(filename);}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat() {return true;}
    virtual bool read();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Eigenvector
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(EigenvectorInfo& ei);
    //bool read(Vector3DBlock& coords, std::vector<PDB::Atom>& atoms, std::vector<PDB::Ter>& ters);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend EigenvectorTextReader& operator>>(EigenvectorTextReader& eigenvectorReader, EigenvectorInfo& info);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    EigenvectorInfo* myEigenvectorInfo;
  };

  //____________________________________________________________________________INLINES

}
#endif /* EIGENVECTOTTEXTREADER_H */
