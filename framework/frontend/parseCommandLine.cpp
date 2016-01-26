#include "parseCommandLine.h"
#include "Report.h"
#include "pmconstants.h"
#include "inputValueDefinitions.h"
#include "Configuration.h"

#include "TopologyFactory.h"
#include "OutputFactory.h"
#include "ForceFactory.h"
#include "IntegratorFactory.h"
#include "HelpTextFactory.h"


#include "Topology.h"
#include "VacuumBoundaryConditions.h"
#include "CubicCellManager.h"
#include "PeriodicBoundaryConditions.h"


using std::vector;
using std::string;
using std::endl;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
using namespace ProtoMol::Constant::ParseCommandLine;
using namespace ProtoMol::Constant::SI;

namespace ProtoMol {

  //________________________________________________________ parseCommandLine
  vector<vector<string> > parseCommandLine (int argc, char **argv, 
					    const Configuration* config,  
					    void (*registerForceExemplarsFunction)(const GenericTopology*)){
    
    if (argc < 2 || (argc >= 2 && (string(argv[1]) =="-h" ||
				   string(argv[1]) =="--help" ))) {  
	report << plain
	       << "Usage: "<<argv[0]<<" ["<<prefix<<"config] <filename> ["
	       << prefix<<"otherconfigs args] [...] ..."<<endl
	       << PROTOMOL_HR << endl
	       << "Other options:" << endl
	       << endl
	       << "-h/--help                This menu" << endl
	       << "-m/--man <keyword>       Man page" << endl
	       << "-m/--man --keywords      Printout of all man page keywords" << endl
	       << "-v/--version             Version number" << endl
	       << "-f/--forces              Printout of all supported forces" << endl
	       << "-c/--constants           Printout of all constants" << endl
	       << "-i/--integrators         Printout of all supported integrators" << endl
	       << "-o/--outputs             Printout of all supported outputs" << endl
	       << "-t/--topologies          Printout of all supported topologies" << endl
	       << "-k/--keywords            Printout of all supported keywords, default values and aliases" << endl
	       << "-u/--units               Printout of all ProtoMol units" << endl
	  	       << PROTOMOL_HR << quit << endr;
    }
    else  if (argc >= 3 && (string(argv[1]) =="-m" ||
			    string(argv[1]) =="--man" )) {  

      HelpTextFactory::registerExemplars(config);
      TopologyFactory::registerHelpText();
      OutputFactory::registerHelpText();
      IntegratorFactory::registerHelpText();

      if(string(argv[2]) == "--keywords"){
	report << plain <<  HelpTextFactory::keywords() << quit << endr;
      }
      else {
	report << plain <<  HelpTextFactory::search(string(argv[2])) << quit << endr;
      }

    }
    else if (string(argv[1]) == "-v" ||
	     string(argv[1]) == "--version") {  
      report << plain << quit << endr;
    }
    else if (string(argv[1]) == "-c" ||
	     string(argv[1]) == "--constants") {  
      report << plain << Constant::print() << quit << endr;
    }
    else if (string(argv[1]) == "-i" ||
	     string(argv[1]) == "--integrators") {  
      report << plain
	     << "Integrators:\n"<<IntegratorFactory::print()<<endl
	     << PROTOMOL_HR << quit << endr;
    }
    else if (string(argv[1]) == "-o" ||
	     string(argv[1]) == "--outputs") {  
      report << plain
	     << "Ouputs:\n"<<OutputFactory::print()<<endl
	     << PROTOMOL_HR << quit << endr;
    }
    else if (string(argv[1]) == "-t" ||
	     string(argv[1]) == "--topologies") {  
      report << plain
	     << "Topologies:\n"<<TopologyFactory::print()<<endl
	     << PROTOMOL_HR << quit << endr;
    }
    else if (string(argv[1]) == "-k" ||
	     string(argv[1]) == "--keywords") {  
      report << plain
	     << (config?config->print():string("Undefined configuration"))<<endl
	     << PROTOMOL_HR << quit << endr;
    }
    else if (string(argv[1]) == "-u" ||
	     string(argv[1]) == "--units") {  
      report << plain
	     << "Units:"<< endl
	     << endl
	     << "Time      : [fs]"<< endl
	     << "Length    : [AA]"<< endl
	     << "Velocity  : [fs/AA]"<< endl
	     << "Energy    : [kcal/mol]"<< endl
	     << "Force     : [kcal/mol AA]"<< endl
	     << "Mass      : [AMU]"<< endl
	     << "Temperatur: [K]"<< endl
	     << "Charge    : [e]"<< endl
	     <<endl
	     << "ProtoMol Constants:"<<endl
	     << endl
	     << PROTOMOL_HR << quit << endr;
    }
    else if (string(argv[1]) == "-f" ||
	     string(argv[1]) == "--forces") {  
      if(registerForceExemplarsFunction == NULL){
	report << plain
	       << "Forces:\n"<<ForceFactory::print()<<endl
	       << PROTOMOL_HR << quit <<endr;
      }
      else {
	for(TopologyFactory::const_iterator itr=TopologyFactory::begin();itr!=TopologyFactory::end();++itr){
	  ForceFactory::unregisterAllExemplars();
	  (*registerForceExemplarsFunction)(*itr);
	  report << plain
		 << "Forces with \'"<<(*itr)->getId()<<"\' topology:\n"
		 << PROTOMOL_HR <<endl	  
		 << ForceFactory::print()<<endl
		 << PROTOMOL_HR <<endr;	  
	}
	report << quit <<endr;
      }
    }

    // Parse input ...
    vector<vector<string> > res;
    bool first = true;
    for (int i = 1; i < argc; ++i){
      string str(argv[i]); 
      if(first && !equalStart(prefix,str)){
	res.resize(res.size()+1);
	res[res.size()-1].push_back(InputConfig::keyword);
      }
      else if(equalStart(prefix,str)){
	res.resize(res.size()+1);
	str = str.substr(prefix.size());
      }
      res[res.size()-1].push_back(str);

      first = false;
    }
    return res;
  }

}

