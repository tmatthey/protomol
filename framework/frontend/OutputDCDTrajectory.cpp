#include "OutputDCDTrajectory.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "DCDTrajectoryWriter.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputDCDTrajectory
  const string  OutputDCDTrajectory::keyword("DCDFile");

  OutputDCDTrajectory::OutputDCDTrajectory():Output(),myDCD(NULL),myMinimalImage(false){}

  OutputDCDTrajectory::OutputDCDTrajectory(const string& filename, int freq, bool minimal):Output(freq),myDCD(new DCDTrajectoryWriter(filename)),myMinimalImage(minimal){}

  OutputDCDTrajectory::~OutputDCDTrajectory(){
    if(myDCD != NULL)
      delete myDCD;
  }

  void OutputDCDTrajectory::doInitialize(){
    if(myDCD == NULL || !myDCD->open())
      report << error <<" Can not open \'"<<(myDCD!=NULL?myDCD->getFilename():"")<<"\' for "<<getId()<<"."<<endr;
  }

  void OutputDCDTrajectory::doRun(int){
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
    if(!myDCD->write(*pos))
       report << error << "Could not write "<<getId()<<" \'"<<myDCD->getFilename()<<"\'."<<endr;    
  }

  void OutputDCDTrajectory::doFinalize(int){
    myDCD->close();
  }

  Output* OutputDCDTrajectory::doMake(string&, const vector<Value>& values) const{
    return (new OutputDCDTrajectory(values[0],values[1],values[2]));
  }

  void OutputDCDTrajectory::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myDCD!=NULL?myDCD->getFilename():"",ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
  }

  bool OutputDCDTrajectory::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    if(config->valid(InputMinimalImage::keyword) && !values[2].valid())
      values[2] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }

}
