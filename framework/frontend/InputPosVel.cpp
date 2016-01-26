#include "InputPosVel.h"
#include "Vector.h"
#include "systemutilities.h"
#include "Configuration.h"
#include "XYZBinReader.h"
#include "XYZReader.h"
#include "PDBReader.h"

using std::vector;
using std::string;

namespace ProtoMol {
  //_________________________________________________________________ File

  InputPosVel::InputPosVel():myFilename(""),myOk(true),myType(InputPosVelType::UNDEFINED){}

  InputPosVel::InputPosVel(const string& filename):myFilename(filename),myOk(isAccessible(filename)),myType(InputPosVelType::UNDEFINED){}

  
  InputPosVel::operator void*() const{
    return (!myOk ? 0 : const_cast<InputPosVel*>(this));
  }

  bool InputPosVel::operator!() const { 
    return !myOk;
  }

  void InputPosVel::setFilename(const string& filename){
    myFilename = filename;
    myOk = true;
    myType = InputPosVelType::UNDEFINED;

  }

  bool InputPosVel::open(){
    myOk = isAccessible(myFilename);
    return myOk;
  }

  bool InputPosVel::open(const string& filename){
    setFilename(filename);
    return open();
  }

  bool InputPosVel::tryFormat(InputPosVelType::Enum type){
    if(type == InputPosVelType::PDB){
      PDBReader reader(myFilename);
      return reader.tryFormat();
    }
    else if(type == InputPosVelType::XYZ){
      XYZReader reader(myFilename);
      return reader.tryFormat();
    }
    else if(type == InputPosVelType::XYZBIN){
      XYZBinReader reader(myFilename);
      return reader.tryFormat();
    }
    else{
      return isAccessible(myFilename);
    }
  }


  InputPosVelType InputPosVel::getType() const{
    return myType;
  }

  InputPosVel& operator>>(InputPosVel& posReader, PDB& pdb){
    posReader.myType = InputPosVelType::UNDEFINED;

    // PDB
    PDBReader pdbReader(posReader.myFilename);
    posReader.myOk = pdbReader.tryFormat();
    if(posReader.myOk){
		posReader.myOk = (pdbReader >> pdb?true:false);
      posReader.myType = InputPosVelType::PDB;
    }

    return posReader;
  }

  InputPosVel& operator>>(InputPosVel& posReader, XYZ& xyz){
    posReader.myType = InputPosVelType::UNDEFINED;

    // XYZ
    XYZReader xyzReader(posReader.myFilename);
    posReader.myOk = xyzReader.tryFormat();
    if(posReader.myOk){
		posReader.myOk = (xyzReader >> xyz?true:false);
      posReader.myType = InputPosVelType::XYZ;
    }
    // PDB
    if(!posReader.myOk){
      PDBReader pdbReader(posReader.myFilename);
      if(pdbReader.tryFormat()){
		  posReader.myOk = (pdbReader >> xyz?true:false);
	posReader.myType = InputPosVelType::PDB;
      }
    }   

    return posReader;
  }

  InputPosVel& operator>>(InputPosVel& posReader, Vector3DBlock& coords){
    posReader.myType = InputPosVelType::UNDEFINED;

    // XYZ
    XYZReader xyzReader(posReader.myFilename);
    posReader.myOk = xyzReader.tryFormat();
    if(posReader.myOk){
		posReader.myOk = (xyzReader >> coords?true:false);
      posReader.myType = InputPosVelType::XYZ;
    }
    // XYZ binary
    if(!posReader.myOk){
      XYZBinReader xyzBinReader(posReader.myFilename);
      if(xyzBinReader.tryFormat()){
		  posReader.myOk = (xyzBinReader >> coords?true:false);      
	posReader.myType = InputPosVelType::XYZBIN;
      }
    }

    // PDB
    if(!posReader.myOk){
      PDBReader pdbReader(posReader.myFilename);
      if(pdbReader.tryFormat()){
		  posReader.myOk = (pdbReader >> coords?true:false);
	posReader.myType = InputPosVelType::PDB;
      }
    }   

    return posReader;
  }

}
