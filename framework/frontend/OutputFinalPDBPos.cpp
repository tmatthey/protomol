#include "OutputFinalPDBPos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "PDBWriter.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalPDBPos
  const string  OutputFinalPDBPos::keyword("finPDBPosFile");

  OutputFinalPDBPos::OutputFinalPDBPos():Output(1),myFilename(""),myMinimalImage(false){}

  OutputFinalPDBPos::OutputFinalPDBPos(const string& filename, bool minimal):Output(1),myFilename(filename),myMinimalImage(minimal){}

  void OutputFinalPDBPos::doFinalize(int step){
    PDBWriter writer;
    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    writer.setComment("Time : "+toString(myCache->time())+", step : "+toString(step)+(myMinimalImage?", minimal Image":"")+".");
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
    if(!writer.write(*pos, myCache->pdb()))
      report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr;    
  }

  Output* OutputFinalPDBPos::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalPDBPos(values[0],values[1]));
  }

  void OutputFinalPDBPos::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
  }

  bool OutputFinalPDBPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputMinimalImage::keyword) && !values[1].valid())
      values[1] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }
}
