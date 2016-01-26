#include "OutputFinalXSC.h"
#include "Configuration.h"
#include "OutputCache.h"      //
#include "stringutilities.h"
#include "iSGIntegrator.h"
#include "Topology.h"          //
#include "XSCWriter.h"
#include "XSC.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalXSC
  const string  OutputFinalXSC::keyword("finXSCFile");

  OutputFinalXSC::OutputFinalXSC():Output(1),myFilename(""){}

  OutputFinalXSC::OutputFinalXSC(const string& filename):Output(1),myFilename(filename){}

  void OutputFinalXSC::doFinalize(int step){
    XSCWriter writer;

    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    writer.setComment("Time : "+toString(myCache->time())+", step : "+toString(step)+".");

    const iSGIntegrator* isg = dynamic_cast<const iSGIntegrator*>(myIntegrator->bottom());
    if(isg == NULL){
      report << error << keyword << " : Level 0 integrator is not an iSGIntegrator!"<< endr;
    }
    else {
      if(!writer.write(isg->getXSC()))
	report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr; 
    }
  }

  Output* OutputFinalXSC::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalXSC(values[0]));
  }

  void OutputFinalXSC::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
  }

}
