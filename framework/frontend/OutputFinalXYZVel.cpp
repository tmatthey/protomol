#include "OutputFinalXYZVel.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZWriter.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalXYZVel
  const string  OutputFinalXYZVel::keyword("finXYZVelFile");

  OutputFinalXYZVel::OutputFinalXYZVel():Output(1),myFilename(""){}

  OutputFinalXYZVel::OutputFinalXYZVel(const string& filename):Output(1),myFilename(filename){}

  void OutputFinalXYZVel::doFinalize(int step){
    XYZWriter writer;
    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    writer.setComment("Time : "+toString(myCache->time())+", step : "+toString(step)+".");
    if(!writer.write(*myVelocities,myTopology->atoms,myTopology->atomTypes))
      report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr;    
  }

  Output* OutputFinalXYZVel::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalXYZVel(values[0]));
  }

  void OutputFinalXYZVel::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
  }

}
