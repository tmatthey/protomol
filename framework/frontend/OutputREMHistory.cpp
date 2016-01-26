#include "OutputREMHistory.h"
#include "inputValueDefinitions.h"
#include "OutputCache.h"

using std::string;
using std::vector;
using std::endl;

namespace ProtoMol {



  const string OutputREMHistoryFile::keyword("REMHistoryFile");

  OutputREMHistoryFile::OutputREMHistoryFile(): OutputFile() {}
  
  OutputREMHistoryFile::OutputREMHistoryFile(const string &filename, int freq, int cacheFreq, int cacheSize, Real closeTime): OutputFile(filename, freq, cacheFreq, cacheSize, closeTime) {}

  OutputREMHistoryFile::~OutputREMHistoryFile() {}

  void OutputREMHistoryFile::doInitialize() {
    open();
    close();
  }

  void OutputREMHistoryFile::doRunCached(int /*step*/) {
    //Empty function block -- data only makes sense at the end of the simulation!
  }

  void OutputREMHistoryFile::doFinalize(int /*step*/) {
    const Real *history = myCache->replicaHistory();
    myBuffer << "This configuration has occupied the following temperatures (in chronological order):" << endl;
    for (unsigned int i = 1; i < history[0]; i++)
      myBuffer << history[i] << endl;
  }

  Output *OutputREMHistoryFile::doMake(string & /*errMesg*/, const vector<Value> &values) const {
    return (new OutputREMHistoryFile(values[0], values[1], values[2], values[3], values[4]));
  }

  void OutputREMHistoryFile::getParameters(vector<Parameter> &parameter) const {
    OutputFile::getParameters(parameter);
  }

  bool OutputREMHistoryFile::adjustWithDefaultParameters(vector<Value> &values, const Configuration *config) const {
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    return checkParameters(values);    
  }
}
