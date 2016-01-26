#include "OutputXYZTrajectoryPos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZTrajectoryWriter.h"
#include "inputValueDefinitions.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputXYZTrajectoryPos
  const string  OutputXYZTrajectoryPos::keyword("XYZPosFile");

  OutputXYZTrajectoryPos::OutputXYZTrajectoryPos():Output(),myXYZ(NULL),myMinimalImage(false){}

  OutputXYZTrajectoryPos::OutputXYZTrajectoryPos(const string& filename, int freq, bool minimal):Output(freq),myXYZ(new XYZTrajectoryWriter(filename)),myMinimalImage(minimal){}

  OutputXYZTrajectoryPos::~OutputXYZTrajectoryPos(){
    if(myXYZ != NULL)
      delete myXYZ;
  }

  void OutputXYZTrajectoryPos::doInitialize(){
    if(myXYZ == NULL || !myXYZ->open())
      report << error <<" Can not open \'"<<(myXYZ!=NULL?myXYZ->getFilename():"")<<"\' for "<<getId()<<"."<<endr;
  }

  void OutputXYZTrajectoryPos::doRun(int){
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);

    if(!myXYZ->write(*pos,myTopology->atoms,myTopology->atomTypes))
       report << error << "Could not write "<<getId()<<" \'"<<myXYZ->getFilename()<<"\'."<<endr;    
  }

  void OutputXYZTrajectoryPos::doFinalize(int){
    myXYZ->close();
  }

  Output* OutputXYZTrajectoryPos::doMake(string&, const vector<Value>& values) const{
    return (new OutputXYZTrajectoryPos(values[0],values[1],values[2]));
  }

  void OutputXYZTrajectoryPos::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myXYZ!=NULL?myXYZ->getFilename():"",ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
  }

  bool OutputXYZTrajectoryPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    if(config->valid(InputMinimalImage::keyword) && !values[2].valid())
      values[2] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }

}
