#include "OutputFinalXYZBinPos.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "XYZBinWriter.h"
#include "inputValueDefinitions.h"
#include "OutputCache.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalXYZBinPos
  const string  OutputFinalXYZBinPos::keyword("finXYZBinPosFile");

  OutputFinalXYZBinPos::OutputFinalXYZBinPos():Output(1),myFilename(""),myMinimalImage(false){}

  OutputFinalXYZBinPos::OutputFinalXYZBinPos(const string& filename, bool minimal):Output(1),myFilename(filename),myMinimalImage(minimal){}

  void OutputFinalXYZBinPos::doFinalize(int){
    XYZBinWriter writer;
    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
    if(!writer.write(*pos))
      report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr;    
  }

  Output* OutputFinalXYZBinPos::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalXYZBinPos(values[0],values[1]));
  }

  void OutputFinalXYZBinPos::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
  }

  bool OutputFinalXYZBinPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputMinimalImage::keyword) && !values[1].valid())
      values[1] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }
}
