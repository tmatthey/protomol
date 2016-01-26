#include "OutputXYZBinTrajectoryVel.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZBinWriter.h"
#include "inputValueDefinitions.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputXYZBinTrajectoryVel
  const string  OutputXYZBinTrajectoryVel::keyword("XYZBinVelFile");

  OutputXYZBinTrajectoryVel::OutputXYZBinTrajectoryVel():Output(),myXYZ(){}

  OutputXYZBinTrajectoryVel::OutputXYZBinTrajectoryVel(const string& filename, int freq):Output(freq),myXYZ(new XYZBinWriter(filename)){}

  OutputXYZBinTrajectoryVel::~OutputXYZBinTrajectoryVel(){
    if(myXYZ != NULL)
      delete myXYZ;
  }

  void OutputXYZBinTrajectoryVel::doInitialize(){
    if(myXYZ == NULL || !myXYZ->open())
      report << error <<" Can not open \'"<<(myXYZ!=NULL?myXYZ->getFilename():"")<<"\' for "<<getId()<<"."<<endr;
  }

  void OutputXYZBinTrajectoryVel::doRun(int){
    if(!myXYZ->write(*myVelocities))
       report << error << "Could not write "<<getId()<<" \'"<<myXYZ->getFilename()<<"\'."<<endr;    
  }

  void OutputXYZBinTrajectoryVel::doFinalize(int){
    myXYZ->close();
  }

  Output* OutputXYZBinTrajectoryVel::doMake(string&, const vector<Value>& values) const{
    return (new OutputXYZBinTrajectoryVel(values[0],values[1]));
  }

  void OutputXYZBinTrajectoryVel::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myXYZ!=NULL?myXYZ->getFilename():"",ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
  }

  bool OutputXYZBinTrajectoryVel::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    return checkParameters(values);
  }

}
