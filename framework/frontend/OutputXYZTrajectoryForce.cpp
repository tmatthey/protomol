#include "OutputXYZTrajectoryForce.h"
#include "Configuration.h"
#include "Integrator.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZTrajectoryWriter.h"
#include "inputValueDefinitions.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputXYZTrajectoryForce
  const string  OutputXYZTrajectoryForce::keyword("XYZForceFile");

  OutputXYZTrajectoryForce::OutputXYZTrajectoryForce():Output(),myXYZ(){}

  OutputXYZTrajectoryForce::OutputXYZTrajectoryForce(const string& filename, int freq):Output(freq),myXYZ(new XYZTrajectoryWriter(filename)){}

  OutputXYZTrajectoryForce::~OutputXYZTrajectoryForce(){
    if(myXYZ != NULL)
      delete myXYZ;
  }

  void OutputXYZTrajectoryForce::doInitialize(){
    if(myXYZ == NULL || !myXYZ->open())
      report << error <<" Can not open \'"<<(myXYZ!=NULL?myXYZ->getFilename():"")<<"\' for "<<getId()<<"."<<endr;
  }

  void OutputXYZTrajectoryForce::doRun(int){
    if(!myXYZ->write(*(myIntegrator->getForces()),myTopology->atoms,myTopology->atomTypes))
       report << error << "Could not write "<<getId()<<" \'"<<myXYZ->getFilename()<<"\'."<<endr;    
  }

  void OutputXYZTrajectoryForce::doFinalize(int){
    myXYZ->close();
  }

  Output* OutputXYZTrajectoryForce::doMake(string&, const vector<Value>& values) const{
    return (new OutputXYZTrajectoryForce(values[0],values[1]));
  }

  void OutputXYZTrajectoryForce::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myXYZ!=NULL?myXYZ->getFilename():"",ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
  }

  bool OutputXYZTrajectoryForce::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    return checkParameters(values);
  }

}
