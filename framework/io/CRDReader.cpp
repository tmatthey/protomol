#include "CRDReader.h"

#include "stringutilities.h"
#include "Report.h"
#include "mathutilities.h"

using std::string;
using std::vector;
using std::endl;
using std::find;
using std::stringstream;
using namespace ProtoMol::Report;

//#define DEBUG_CRD

namespace ProtoMol {
  //_________________________________________________________________CRDReader

  CRDReader::CRDReader():
    Reader(),
    myCoords(NULL){}

  CRDReader::CRDReader(const std::string& filename):
    Reader(filename),
    myCoords(NULL){}
			//,myAtoms(NULL){}


  CRDReader::~CRDReader(){
    if(myCoords != NULL)
      delete myCoords;
    /*if(myAtoms != NULL)
      delete myAtoms;*/
  }

  bool CRDReader::tryFormat(){
    if(!open())
      return false;
	else return true;	
	/*
    do {
      string record,str;
      record = getline();
      stringstream ss(record);
      ss >> str;
      if("ATOM" == str){
	close();
	return true;
      }
    } while (!myFile.eof());
    myFile.setstate(std::ios::failbit);
    close();        
    return false;*/
  }

  bool CRDReader::read() {
    if(myCoords == NULL)
      myCoords = new Vector3DBlock();
    return read(*myCoords);
  }

  /*bool CRDReader::read(CRD& pdb){
    return read(pdb.coords,pdb.atoms);
  }*/
  
  bool CRDReader::read(Vector3DBlock& coords) {
    if(!tryFormat())
      return false;
    if(!open())
      return false;
    coords.clear();
    //atoms.clear();
	//
	
	//read the coord file.
	string line, str;
	int numRecords;
	stringstream ss;
	string xcoord, ycoord, zcoord;

	line = getline();//read the header
	while(line.empty()) line = getline(); //ignore blank lines
	
	//Ignore the lines in the header...then the first line of data should
	//have the total number of atoms.
	std::cout<<"LINE = "<<line<<endl;	
	ss << line;
	ss >> str;
	while(!isInt(str)) {
			line = getline();
		ss.clear();		
	std::cout<<"LINE = "<<line<<endl;	
			ss << line;
			ss >> str;
	}	
	numRecords = toInt(str); //total number of atoms
		for(int i=0;i<numRecords;i++) {
			myFile >> xcoord;
			myFile >> ycoord;
			myFile >> zcoord;
			coords.push_back(Vector3D(toReal(xcoord), toReal(ycoord), toReal(zcoord)));
		}	

    // Now we want to read data in until the record name is "END", then stop.
#if 0
    int big = 0;
    int toBig = 0;
    myComment = "";
    do{
      string record(getline());
      if(equalStart("END",record))
	break;
      if(equalStart("REMARK",record) || record.size() < 1){
	record = removeBeginEndBlanks(record);
	if(!record.empty())
	  myComment += (myComment.empty()?"":"\n")+record;
	continue;
      }
      if(equalStart("ATOM",record)){
	record = getRightFill(record,80);
	string str = removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_RES_SEQ,CRD::CRDAtom::L_RES_SEQ));
	int resSeq = toInt(str);
	if(!isInt(str)){
	  if(str.size() == 4 && isInt(string(&str[1],&str[4])) && str[0] >= 'A' && str[0] <= 'Z'){
	    resSeq = (str[0] - 'A') * 1000 + 10000 + toInt(string(&str[1],&str[4]));
	    ++big;
	  }
	  else {
	    resSeq = -1;
	    ++toBig;
	  }
	}
	    
	atoms.push_back(CRD::CRDAtom(removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_RECORD_NAME,CRD::CRDAtom::L_RECORD_NAME)),
				     toInt(record.substr(CRD::CRDAtom::S_SERIAL,CRD::CRDAtom::L_SERIAL)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_ATOM_NAME,CRD::CRDAtom::L_ATOM_NAME)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_ALT_LOC,CRD::CRDAtom::L_ALT_LOC)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_RES_NAME,CRD::CRDAtom::L_RES_NAME)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_CHAIN_ID,CRD::CRDAtom::L_CHAIN_ID)),
				     resSeq,
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_I_CODE,CRD::CRDAtom::L_I_CODE)),
				     toReal(record.substr(CRD::CRDAtom::S_OCCUP,CRD::CRDAtom::L_OCCUP)),
				     toReal(record.substr(CRD::CRDAtom::S_TEMP,CRD::CRDAtom::L_TEMP)),
				     toInt(record.substr(CRD::CRDAtom::S_FOOT_NOTE,CRD::CRDAtom::L_FOOT_NOTE)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_SEG_ID,CRD::CRDAtom::L_SEG_ID)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_ELEMENT_SYMBOL,CRD::CRDAtom::L_ELEMENT_SYMBOL)),
				     removeBeginEndBlanks(record.substr(CRD::CRDAtom::S_CHARGE,CRD::CRDAtom::L_CHARGE)),
				     0));
	coords.push_back(Vector3D(toReal(record.substr(CRD::CRDAtom::S_X,CRD::CRDAtom::L_X)),
				  toReal(record.substr(CRD::CRDAtom::S_Y,CRD::CRDAtom::L_Y)),
				  toReal(record.substr(CRD::CRDAtom::S_Z,CRD::CRDAtom::L_Z))));
      }
      else{
	report << recoverable << "[CRD::read] Record unknow:\'"<<record<<"\'."<<endr;
      }
    }
    while (!(myFile.eof()));

    if(big > 0)
      report << hint << "[CRD::read] Found "<<big<<" X-Plor residue number(s) starting with a character."<<endr;
    if(toBig > 0){
      report << recoverable << "[CRD::read] Found "<<toBig<<" non interger/X-Plor residue number(s)."<<endr;
    }

    close();      
    return !myFile.fail();
#endif
	return true;
  }

  /*Vector3DBlock* CRDReader::orphanCoords(){
    Vector3DBlock* tmp = myCoords;
    myCoords = NULL;
    return tmp;
  }*/

  /*std::vector<CRD::CRDAtom>* CRDReader::orphanAtoms(){
    std::vector<CRD::CRDAtom>* tmp = myAtoms;
    myAtoms = NULL;
    return tmp;
  }*/

  Vector3DBlock CRDReader::getCRD() const{
    Vector3DBlock res;
    if(myCoords != NULL)
      res = (*myCoords);
    /*if(myAtoms != NULL)
      res.atoms = (*myAtoms);*/
    return res;
  }
/*
  CRDReader& operator>>(CRDReader& pdbReader, CRD& pdb){
    pdbReader.read(pdb.coords,pdb.atoms);
    return pdbReader;
  }
  */

  CRDReader& operator>>(CRDReader& crdReader, Vector3DBlock& coords){
    /*if(pdbReader.myAtoms == NULL)
      pdbReader.myAtoms = new std::vector<CRD::CRDAtom>();*/
    crdReader.read(coords);
    return crdReader;
  }
/*
  CRDReader& operator>>(CRDReader& pdbReader, XYZ& xyz){
    if(pdbReader.myAtoms == NULL)
      pdbReader.myAtoms = new std::vector<CRD::CRDAtom>();
    if(pdbReader.read(xyz.coords,*pdbReader.myAtoms)){
      xyz.names.resize(xyz.coords.size());
      for(unsigned int i=0;i<xyz.coords.size();++i)
	xyz.names[i] = (*pdbReader.myAtoms)[i].elementName;
    }
    return pdbReader;
  }
*/
 void CRDReader::writeData()
 {
	Vector3DBlock v = *myCoords;
	for(unsigned int i=0;i<myCoords->size();i++) {
		const Vector3D& c(v[i]);
		std::cout<<c.x<<"\t"<<c.y<<"\t"<<c.z<<std::endl;
	}
	std::cout<<myCoords->size()<<std::endl;
 }
}

