#include "OutputScreen.h"
#include "OutputCache.h"
#include "Integrator.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "inputValueDefinitions.h"
#include "topologyutilities.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputScreen
  const string  OutputScreen::keyword("Screen");

  OutputScreen::OutputScreen():Output(),myUnit("fs"),myFactor(1.0){}
  OutputScreen::OutputScreen(int freq):Output(freq),myUnit("fs"),myFactor(1.0){}

  void OutputScreen::doInitialize(){
    Real step = myIntegrator->getTimestep()*std::max(1,std::min(myOutputFreq,(int)(*myConfig)[InputNumsteps::keyword]));
    if(step >= 1e13){
      myUnit   = "s";
      myFactor = 1e-15;
    }
    else if(step >= 1e10){
      myUnit   = "ms";
      myFactor = 1e-12;
    }
    else if(step >= 1e7){
      myUnit   = "us";
      myFactor = 1e-9;
    }
    else if(step >= 1e4){
      myUnit   = "ns";
      myFactor = 1e-6;
    }
    else if(step >= 1e1){
      myUnit   = "ps";
      myFactor = 1e-3;
    }
  }

  void OutputScreen::doRun(int step){
    report <<plain <<"Step : ";
    report.setf(std::ios::right);
    report.width(10);
    report <<step<<", Time : ";
    report.width(18);
    report.setf(std::ios::showpoint|std::ios::fixed);
    report.precision(3);
    report << myCache->time()*myFactor<<" ["<<myUnit<<"], TE : ";
    report.precision(4);
    report.width(16);
    report <<myCache->totalEnergy()<<" [kcal/mol]";
    report << ", T : ";
    report.precision(4);
    report.width(10);
    report <<myCache->temperature()<<" [K]";
    report << ", V : ";
    report.precision(2);
    report.width(16);
    report <<myCache->volume()<<" [AA^3]"<<endr;
    report.reset();
  }

  Output* OutputScreen::doMake(string&, const  vector<Value>& values) const{
    return (new OutputScreen(values[1]));
  }

  bool OutputScreen::isIdDefined(const Configuration* config) const{
    return (config->valid("outputFreq") && 
	    !config->empty(getId()) &&
	    (!config->valid(getId()) || ((*config)[getId()] == true))
	    );
  }

  void OutputScreen::getParameters(vector<Parameter> & parameter) const{
    parameter.push_back(Parameter(getId(),Value(true),true));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
  }

  bool OutputScreen::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    return checkParameters(values);
  }

}
