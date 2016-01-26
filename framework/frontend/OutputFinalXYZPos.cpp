#include "OutputFinalXYZPos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZWriter.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalXYZPos
  const string  OutputFinalXYZPos::keyword("finXYZPosFile");

  OutputFinalXYZPos::OutputFinalXYZPos():Output(1),myFilename(""),myMinimalImage(false){}

  OutputFinalXYZPos::OutputFinalXYZPos(const string& filename, bool minimal):Output(1),myFilename(filename),myMinimalImage(minimal){}

  void OutputFinalXYZPos::doFinalize(int step){
    XYZWriter writer;
    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
    writer.setComment("Time : "+toString(myCache->time())+", step : "+toString(step)+(myMinimalImage?", minimal Image":"")+".");
    if(!writer.write(*pos,myTopology->atoms,myTopology->atomTypes))
      report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr;    
  }

  Output* OutputFinalXYZPos::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalXYZPos(values[0],values[1]));
  }

  void OutputFinalXYZPos::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
  }

  bool OutputFinalXYZPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputMinimalImage::keyword) && !values[1].valid())
      values[1] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }
}
