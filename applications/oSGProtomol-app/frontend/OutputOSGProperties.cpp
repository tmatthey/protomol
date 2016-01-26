#include "OutputOSGProperties.h"
#include "Configuration.h"
#include "GenericTopology.h" 
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "OutputCache.h" 
#include "inputValueDefinitions.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::setw;
using std::endl;
using std::flush;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;

namespace ProtoMol {
  //________________________________________________________ Output
  const string  OutputOSGProperties::keyword("oSGPropertiesfile");

  OutputOSGProperties::OutputOSGProperties():OutputFile(),myDoMolecularTemperature(true){}
  OutputOSGProperties::OutputOSGProperties(const string& filename, 
					   int freq, 
					   int cacheFreq, 
					   int cacheSize,
					   Real closeTime, 
					   bool doMolTemp):OutputFile(filename,freq,cacheFreq,cacheSize,closeTime),
							   myDoMolecularTemperature(doMolTemp){}

  void OutputOSGProperties::doInitialize(){
    ofstream oSGPropertiesHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
    if(!oSGPropertiesHeaderFile)
      report << error <<" Can not open \'"<<myFilename<<".header\' for "<<getId()<<"."<<endr;

    oSGPropertiesHeaderFile << setw(18)
			    << "Time(fs)" << " "
			    << setw(14)
                            << "E_kinetic" << " "
                            << setw(14)
			    << "E_potential" << " "
			    << setw(14)
			    << "Temperature" << " ";
    if(myEnergies->virial())
      oSGPropertiesHeaderFile << setw(14)
			      << "Pressure(bar)" << " ";
			      
    if(myEnergies->molecularVirial())
      oSGPropertiesHeaderFile << setw(14)
			      << "Mol_Pres(bar)" << " ";
    oSGPropertiesHeaderFile << setw(12)
			    << "Volume(A^3)";
    for (unsigned int i=0; i<myTopology->iSGNumMols.size(); i++) {
      oSGPropertiesHeaderFile << setw(9)
			      << "N[" << i << "]";}
    oSGPropertiesHeaderFile << setw(12) 
			    << "deltaCQ"
			    << setw(14)
			    << "Steps"
                            << setw(14)
                            << "T_lambda(K)"
                            << endl;
    
    oSGPropertiesHeaderFile.close();
    open();
    close();

  }

  void OutputOSGProperties::doRunCached(int){

    if (myEnergies->output()) {

      myBuffer << resetiosflags(std::ios::showpoint |  std::ios::fixed | std::ios::floatfield)
	       << setw(18)
	       << setprecision(2)
	       << setiosflags(std::ios::showpoint | std::ios::fixed)
	       << myCache->time() << " "
	       << resetiosflags(std::ios::showpoint |  std::ios::fixed | std::ios::floatfield)
	       << setiosflags(std::ios::floatfield)
	       << setprecision(8)
	       << setw(14)
	       << myCache->kineticEnergy()
	       << setw(14)
	       << myCache->potentialEnergy()
	       << setw(14)
	       << myCache->temperature()
	       << setw(14)
	       << myCache->pressure()
	       << setw(14)
	       << myCache->molecularPressure()
	       << setw(15)
	       << myCache->volume() << "   ";
      for (unsigned int i=0; i<myTopology->iSGNumMols.size(); i++) {
	myBuffer << setw(10)
		 << myTopology->iSGNumMols[i];}
      myBuffer  << resetiosflags(std::ios::showpoint |  std::ios::fixed | std::ios::floatfield)
		<< setiosflags(std::ios::floatfield)
		<< setprecision(8)  
		<< setw(16) 
		<< (*myEnergies)[ScalarStructure::CQFLUCTUATION]
		<< setw(10)
		<< (*myEnergies)[ScalarStructure::DELTATIME]
                << setw(14)
                << (*myEnergies)[ScalarStructure::LAMBDA_TEMPERATURE] << endl;

    }

  }


  Output* OutputOSGProperties::doMake(string&, const vector<Value>& values) const{
    return (new OutputOSGProperties(values[0],values[1],values[2],values[3],values[4],values[5]));
  }

  void OutputOSGProperties::getParameters(vector<Parameter> &parameter) const{
    OutputFile::getParameters(parameter);
    parameter.push_back(Parameter("molecularTemperature",Value(myDoMolecularTemperature),true));
  }


}
