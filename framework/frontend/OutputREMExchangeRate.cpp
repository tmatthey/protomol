#include "OutputREMExchangeRate.h"
#include "inputValueDefinitions.h"
#include "OutputCache.h"

using std::string;
using std::vector;
using std::endl;

namespace ProtoMol {



  const string OutputREMExchangeRate::keyword("REMExchangeRatesFile");

  OutputREMExchangeRate::OutputREMExchangeRate(): OutputFile() {}
  
  OutputREMExchangeRate::OutputREMExchangeRate(const string &filename, int freq, int cacheFreq, int cacheSize, Real closeTime): OutputFile(filename, freq, cacheFreq, cacheSize, closeTime) {}

  OutputREMExchangeRate::~OutputREMExchangeRate() {}

  void OutputREMExchangeRate::doInitialize() {
    open();
    close();
  }

  void OutputREMExchangeRate::doRunCached(int /*step*/) {
    //Empty function block -- data only makes sense at the end of the simulation!
  }

  void OutputREMExchangeRate::doFinalize(int /*step*/) {
    if (Parallel::getId() != 0) return;
    const vector<Real> exchangeRates = myCache->REMRates();
    for (unsigned int i = 0; i < exchangeRates.size(); i++)
      myBuffer << exchangeRates[i] << endl;
  }

  Output *OutputREMExchangeRate::doMake(string & /*errMesg*/, const vector<Value> &values) const {
    return (new OutputREMExchangeRate(values[0], values[1], values[2], values[3], values[4]));
  }

  void OutputREMExchangeRate::getParameters(vector<Parameter> &parameter) const {
    OutputFile::getParameters(parameter);
  }

  bool OutputREMExchangeRate::adjustWithDefaultParameters(vector<Value> &values, const Configuration *config) const {
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    return checkParameters(values);    
  }
}
