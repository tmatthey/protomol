/*  -*- c++ -*-  */
#ifndef CRDREADER_H
#define CRDREADER_H

#include "Reader.h"
#include "Vector3DBlock.h"

namespace ProtoMol {

  //_________________________________________________________________CRDReader
  class CRDReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    CRDReader();
    explicit CRDReader(const std::string& filename);
    virtual ~CRDReader();

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
    virtual bool tryFormat();
    virtual bool read();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CRD
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //bool read(CRD& pdb);
    //bool read(Vector3DBlock& coords, std::vector<CRD::CRDAtom>& atoms);
    bool read(Vector3DBlock& coords);

    Vector3DBlock getCRD() const;
    //Vector3DBlock* orphanCoords();
    //std::vector<CRD::CRDAtom>* orphanAtoms();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //friend CRDReader& operator>>(CRDReader& pdbReader, CRD& pdb);
    //friend CRDReader& operator>>(CRDReader& pdbReader, XYZ& xyz);
    friend CRDReader& operator>>(CRDReader& pdbReader, Vector3DBlock& coords);

	void writeData();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Vector3DBlock* myCoords;
    //std::vector<CRD::CRDAtom>* myAtoms;
  };

  //____________________________________________________________________________INLINES

}
#endif /* CRDREADER_H */
