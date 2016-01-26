#include "OutputDiffusion.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "OutputCache.h"
#include "Integrator.h"
#include "inputValueDefinitions.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;

namespace ProtoMol {
  //________________________________________________________ Output
  const string  OutputDiffusion::keyword("diffusionFile");

  OutputDiffusion::OutputDiffusion(): OutputFile(){}
  
  OutputDiffusion::OutputDiffusion(const string& filename, 
				   int freq, 
				   int cacheFreq, 
				   int cacheSize,
				   Real closeTime): OutputFile(filename,freq,cacheFreq,cacheSize,closeTime){}

  void OutputDiffusion::doInitialize(){
    ofstream diffusionHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
    if(!diffusionHeaderFile)
      report << error <<" Can not open \'"<<myFilename<<".header\' for "<<getId()<<"."<<endr;

    diffusionHeaderFile << setw(18)
			<< "Time(fs)" << " "
			<< setw(24)
			<< "Diffusion(10^-5cm^2/s)" << " "
			<< setw(14)
			<< "Volume(AA^3)" << " "
			<< setw(14)
			<< "Temperature(K)" << " "
			<< setw(14)
			<< "Density(kgm^-3*1e3)";
    diffusionHeaderFile << endl;
    diffusionHeaderFile.close();

    open();
    close();
  }
  
  void OutputDiffusion::doRunCached(int step){
    Real diffusion = 10000.0*myCache->diffusion()*(first()?0.0:1.0/((step-getFirstStep())*myIntegrator->getTimestep()));

    myBuffer << resetiosflags(std::ios::showpoint |  std::ios::fixed | std::ios::floatfield)
	     << setw(18)
	     << setprecision(2)
	     << setiosflags(std::ios::showpoint |
			    std::ios::fixed)
	     << myCache->time() << " "
	     << resetiosflags(std::ios::showpoint |
			      std::ios::fixed)
	     << setiosflags(std::ios::floatfield)
	     << setprecision(8)
	     << setw(24)
	     << diffusion << " "
	     << setw(14)
	     << myCache->volume()<< " "
	     << setw(14)
	     << myCache->temperature() << " "
	     << setw(14)
	     << myCache->density();
    myBuffer << endl;
  }

  Output* OutputDiffusion::doMake(string&, const vector<Value>& values) const{
    return (new OutputDiffusion(values[0],values[1],values[2],values[3],values[4]));
  }


}
