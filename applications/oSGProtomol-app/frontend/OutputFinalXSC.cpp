#include "OutputFinalXSC.h"
#include "Configuration.h"
#include "OutputCache.h"      //
#include "stringutilities.h"
#include "Integrator.h"
#include "oSGIntegrator.h"
#include "muVTIntegrator.h"
#include "NfPTIntegrator.h"
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

    const oSGIntegrator* osg;
    if (myIntegrator->getIdNoAlias() == "muVTVerlet" )
      osg = dynamic_cast<const muVTIntegrator*>( myIntegrator->bottom() );
    else if (myIntegrator->getIdNoAlias() == "NfPTVerlet")
      osg = dynamic_cast<const NfPTIntegrator*>( myIntegrator->bottom() );

    XSC myXSC;
    myXSC = osg->getXSC();

    if(osg == NULL){
      report << error << keyword << " : Level 0 integrator is not an oSGIntegrator!"<< endr;
    }
    else {
      if (!writer.write(myXSC)) 
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
