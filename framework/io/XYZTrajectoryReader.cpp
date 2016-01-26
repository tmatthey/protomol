#include "XYZTrajectoryReader.h"

#include "Report.h"
#include "stringutilities.h"

using std::string;
using std::vector;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________XYZTrajectoryReader

  XYZTrajectoryReader::XYZTrajectoryReader():Reader(),myCoords(NULL),myNames(NULL),myFirst(true){}

  XYZTrajectoryReader::XYZTrajectoryReader(const std::string& filename):Reader(filename),myCoords(NULL),myNames(NULL),myFirst(true){}

  XYZTrajectoryReader::~XYZTrajectoryReader(){
    if(myCoords != NULL)
      delete myCoords;
    if(myNames != NULL)
      delete myNames;
  }

  bool XYZTrajectoryReader::tryFormat(){
    if (!open())
      return false;

    // Number of frames and atoms
    unsigned int frames,count=0;
    myFile >> frames >> count;

    // Atoms
    string str;
    Real x;
    myFile >> str >> x >> x >> x;
    close();
    return !myFile.fail();
  }

  bool XYZTrajectoryReader::read() {
    if(myCoords == NULL)
      myCoords = new Vector3DBlock();
    if(myNames == NULL)
      myNames = new vector<string>();
    return read(*myCoords,*myNames);
  }

  bool XYZTrajectoryReader::read(XYZ& xyz) {
    return read(xyz.coords,xyz.names);
  }

  bool XYZTrajectoryReader::read(Vector3DBlock& coords, std::vector<std::string>& names){
    if(myFirst){
      if(!open())
	return false;
      myFirst = false;    
      
      // Number of frames;
      unsigned int n=0;
      myFile >> n;
      if(n == 0 && myFile.good())
	return true;
    }

    // Number of atoms;
    unsigned int n=0;
    myFile >> n;
    if(myFile.fail()){
      close();
      return false;
    }

    coords.resize(n);
    names.resize(n);

    // Read atoms
    for(unsigned int i = 0; i < n && !myFile.fail(); ++i){
      myFile >> names[i] >> coords[i].x >> coords[i].y >> coords[i].z;
    }
    return !myFile.fail();
  }

  XYZ XYZTrajectoryReader::getXYZ() const{
    XYZ res;
    if(myCoords != NULL)
      res.coords = (*myCoords);
    if(myNames != NULL)
      res.names = (*myNames);
    return res;
  }

  Vector3DBlock* XYZTrajectoryReader::orphanCoords(){
    Vector3DBlock* tmp = myCoords;
    myCoords = NULL;
    return tmp;

  }
  vector<string>* XYZTrajectoryReader::orphanNames(){
    vector<string>* tmp = myNames;
    myNames = NULL;
    return tmp;
  }

  XYZTrajectoryReader& operator>>(XYZTrajectoryReader& xyzReader, XYZ& xyz){
    xyzReader.read(xyz.coords,xyz.names);
    return xyzReader;
  }

  XYZTrajectoryReader& operator>>(XYZTrajectoryReader& xyzReader, Vector3DBlock& coords){
    if(xyzReader.myNames == NULL)
      xyzReader.myNames = new vector<string>();
    xyzReader.read(coords,*xyzReader.myNames);
    return xyzReader;
  }

}
