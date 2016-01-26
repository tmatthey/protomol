#include "ConfigurationReader.h"
#include "InputPosVel.h"
#include "TimerStatistic.h"

#include "PSFWriter.h"
#include "PARWriter.h"


#include "pmconstants.h"
#include "inputValueDefinitions.h"
#include "mathutilities.h"
#include "parseCommandLine.h"
#include "protomol.h"
#include "stringutilities.h"
#include "systemutilities.h"

using std::endl;
using std::string;
using std::vector;

namespace ProtoMol {
  declareInputValue(OutputPAR, STRING, NOCONSTRAINTS);
  declareInputValue(OutputPSF, STRING, NOCONSTRAINTS);
  
  defineInputValue(OutputPAR, "PARFile");
  defineInputValue(OutputPSF, "PSFFile");
}

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
//_____________________________________________________________________ protomol

int main(int argc, char **argv) {

  // Redirect all output to a file
  //std::ofstream out("out.txt");
  //report.setStream(&(out));

  // Redirect all output to std::cout
  //report.setStream(&(std::cout));


  TimerStatistic::timer[TimerStatistic::WALL].start();

  report << plain <<endl
	 <<"   #                                     #####  ######  #     #"<<endl
	 <<"  # #    #    #  #####   ######  #####  #     # #     # ##   ##"<<endl
	 <<" #   #   ##  ##  #    #  #       #    #       # #     # # # # #"<<endl
	 <<"#     #  # ## #  #####   #####   #    #  #####  ######  #  #  #"<<endl
	 <<"#######  #    #  #    #  #       #####  #       #       #     #"<<endl
	 <<"#     #  #    #  #    #  #       #   #  #       #       #     #"<<endl
	 <<"#     #  #    #  #####   ######  #    # ####### #       #     #"<<endl
         << PROTOMOL_HR << endl
	 << PACKAGE_CITE <<endl
         << PROTOMOL_HR << endl
         << "ProtoMol Version "<< PACKAGE_VERSION <<endl<<endl
         << PACKAGE_WHOAMI <<" "<<__DATE__<<" "<<__TIME__<<endl
         << PACKAGE_UNAME<<endl
         << PACKAGE_COMPILER <<endl
         << sizeof(void*)*8 <<"-bit"<<endl
         << PACKAGE_COMPILER_VERSION 
         << (string(PACKAGE_COMPILER_VERSION).empty() ? "":"\n")
         << PROTOMOL_HR << endl
         << "Information:"<<endl
         << endl
         << PACKAGE_BUGREPORT <<endl
         << PACKAGE_HOMEPAGE <<endl
         << PROTOMOL_HR << endr;


  //
  // Configure Configuration
  //
  Configuration config;

  InputConfig::registerConfiguration(&config);
  InputPositions::registerConfiguration(&config);
#ifdef NDEBUG
  InputDebug::registerConfiguration(&config,0); 
#else
  InputDebug::registerConfiguration(&config,1); 
#endif
  OutputPSF::registerConfiguration(&config);
  OutputPAR::registerConfiguration(&config);


  //
  // Configuration
  //
  // Parse command line
  vector<vector<string> > args = parseCommandLine(argc,argv,&config);
  if(config.set(InputConfig::keyword,args)){
    // Read config file
    ConfigurationReader configReader;
    if(!configReader.open(config[InputConfig::keyword]))
      report << error << "Can't open configuration file \'"<<config[InputConfig::keyword]<<"\'."<<endr;
    if(!(configReader >> config))
      report << error << "Could not read configuration file \'"<<config[InputConfig::keyword]<<"\'."<<endr;
    changeDirectory(config[InputConfig::keyword]);
  }
  else {
    report << error << "Undefined configuration file."<<endr;
  }
  // Overwrite configuration with command line values
  config.set(args);


  // Check if configuration is complete
  string errMsg;
  if(!config.validConfiguration(errMsg))
    report << plain << endl << errMsg <<endr; 

  report << reportlevel((int)config[InputDebug::keyword]);

  // 
  // Positions
  //
  InputPosVel reader;
  if(!reader.open(config[InputPositions::keyword]))
    report << error << "Can't open position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;
  Vector3DBlock positions;
  vector<PDB::Atom> pdbAtoms;
  if(reader.tryFormat(InputPosVelType::PDB)){
    PDB pdb;
    if(!(reader >> pdb))
      report << error << "Could not parse PDB position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;
    swap(positions,pdb.coords);
    swap(pdbAtoms,pdb.atoms);
  }
  else if(!(reader >> positions))
    report << error << "Could not parse position file \'"<<config[InputPositions::keyword]
	   <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;
  report << plain << "Using "<<reader.getType()<<" posfile \'"<<config[InputPositions::keyword]<<"\' ("<<positions.size()<<")." << endr;







  // Print
  // Configuration
  report << plain  << PROTOMOL_HR << endr;
  report << debug(10) << PROTOMOL_HR << "\nConfiguration:\n" << config.print() <<endr;



  TimerStatistic::timer[TimerStatistic::RUN].stop();

  // Clean up


  TimerStatistic::timer[TimerStatistic::WALL].stop();
  
  report << allnodesserial << plain <<"Timing" 
	 <<" : "<<TimerStatistic()<<"."<<endr;

  return 0;
}
