#include "ConfigurationReader.h"
#include "InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"
#include "TimerStatistic.h"

#include "TopologyFactory.h"
#include "ForceFactory.h"
#include "IntegratorFactory.h"
#include "OutputFactory.h"

#include "registerForceExemplars.h"
#include "registerIntegratorExemplars.h"
#include "registerTopologyExemplars.h"
#include "registerOutputExemplars.h"

#include "OutputCollection.h"

#include "Parameter.h"
#include "ScalarStructure.h"
#include "Parallel.h"

#include "ModifierRemoveLinearMomentum.h"
#include "ModifierRemoveAngularMomentum.h"

#include "buildTopology.h"
#include "pmconstants.h"
#include "inputValueDefinitions.h"
#include "mathutilities.h"
#include "parseCommandLine.h"
#include "protomol.h"
#include "stringutilities.h"
#include "systemutilities.h"
#include "topologyutilities.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
using std::endl;
using std::string;
using std::vector;
//_____________________________________________________________________ protomol

int main(int argc, char **argv) {

  // Redirect all output to a file
  //std::ofstream out("out.txt");
  //report.setStream(&(out));

  // Redirect all output to std::cout
  //report.setStream(&(std::cout));

  Parallel::init(argc,argv);

  TimerStatistic::timer[TimerStatistic::WALL].start();

  report << plain <<endl
	 <<" #####"<<endl
	 <<"#     #   ####   #    #  #        ####   #    #  #####"<<endl
	 <<"#        #    #  #    #  #       #    #  ##  ##  #    #"<<endl
	 <<"#        #    #  #    #  #       #    #  # ## #  #####"<<endl
	 <<"#        #    #  #    #  #       #    #  #    #  #    #"<<endl
	 <<"#     #  #    #  #    #  #       #    #  #    #  #    #"<<endl
	 <<" #####    ####    ####   ######   ####   #    #  #####"<<endl
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

  InputTemperature::registerConfiguration(&config);
  InputConfig::registerConfiguration(&config);
  InputBoundaryConditions::registerConfiguration(&config);
  InputCellManager::registerConfiguration(&config);
  InputIntegrator::registerConfiguration(&config);
  InputSeed::registerConfiguration(&config,getTimerSeed());
  InputNumsteps::registerConfiguration(&config);
  InputFirststep::registerConfiguration(&config);
#ifdef NDEBUG
  InputDebug::registerConfiguration(&config,0); 
#else
  InputDebug::registerConfiguration(&config,1); 
#endif
  InputPositions::registerConfiguration(&config);
  InputVelocities::registerConfiguration(&config);
  InputPSF::registerConfiguration(&config);
  InputPAR::registerConfiguration(&config);
  InputRemoveLinearMomentum::registerConfiguration(&config,-1);
  InputRemoveAngularMomentum::registerConfiguration(&config,-1);
  InputUseBarrier::registerConfiguration(&config,Parallel::isBarrier());
  InputParallelPipe::registerConfiguration(&config,Parallel::getPipeSize());
  InputParallelMode::registerConfiguration(&config,Parallel::getMode().getString());
  InputMaxPackages::registerConfiguration(&config,Parallel::getMaxPackages());
  InputPDBScaling::registerConfiguration(&config,true);
  InputOutputfreq::registerConfiguration(&config,1);
  InputDihedralMultPSF::registerConfiguration(&config);
  InputVirialCalc::registerConfiguration(&config,false);
  InputMolVirialCalc::registerConfiguration(&config,false);
  InputOutput::registerConfiguration(&config,true);

  //
  // Configure OutputFactory
  //
  registerOutputExemplars();
  OutputFactory::registerAllExemplarsConfiguration(&config);
  //
  // Configure TopologyFactory
  //
  registerTopologyExemplars();
  TopologyFactory::registerAllExemplarsConfiguration(&config);

  //
  // Configure IntegratorFactory
  //
  registerIntegratorExemplars();

  //
  // Configuration
  //
  // Parse command line
  vector<vector<string> > args = parseCommandLine(argc,argv,&config,registerForceExemplars);
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

  // Fix for old topology definition
  if (!config[GenericTopology::keyword].valid()){
    config[GenericTopology::keyword] = config[InputBoundaryConditions::keyword].getString()+config[InputCellManager::keyword].getString();
  }



  // Check if configuration is complete
  string errMsg;
  if(!config.validConfiguration(errMsg))
    report << plain << endl << errMsg <<endr; 
  if(!config[InputFirststep::keyword].valid())
    report << error << "Firststep undefined."<<endr;
  if(!config[InputNumsteps::keyword].valid())
    report << error << "Numsteps undefined."<<endr;


  report << reportlevel((int)config[InputDebug::keyword]);

  //
  // Configure Parallel
  //
  if (config[InputUseBarrier::keyword].valid())
    Parallel::setBarrier(config[InputUseBarrier::keyword]);
  if (config[InputParallelPipe::keyword].valid())
    Parallel::setPipeSize(config[InputParallelPipe::keyword]);
  if (config[InputParallelMode::keyword].valid())
    Parallel::setMode(config[InputParallelMode::keyword].getString());
  if (config[InputMaxPackages::keyword].valid())
    Parallel::setMaxPackages(config[InputMaxPackages::keyword]);

  // 
  // Seed, make sure that all nodes have the same seed
  // 
  report << plain  << PROTOMOL_HR << endr;
  int seed = config[InputSeed::keyword];
  Parallel::bcast(seed);
  config[InputSeed::keyword] = seed;
  randomNumber(seed);
  report << plain << "Using seed "<<seed<<"."<<endr;
  report << plain << "Using report level "<<config[InputDebug::keyword]<<"."<<endr;

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

  // 
  // Velocities
  //
  Vector3DBlock velocities;
  if(config.valid(InputVelocities::keyword)){
    if(!reader.open(config[InputVelocities::keyword]))
      report << error << "Can't open velocity file \'"<<config[InputVelocities::keyword]<<"\'."<<endr;
    if(!(reader >> velocities))
      report << error << "Could not parse velocity file \'"<<config[InputVelocities::keyword]
	     <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;
    report << plain << "Using "<<reader.getType()<<" velfile \'"<<config[InputVelocities::keyword]<<"\' ("<<velocities.size()<<")." << endr;
    if(reader.getType() == "PDB" && (bool)config[InputPDBScaling::keyword]){
      for(unsigned int i=0;i<velocities.size();i++)
	velocities[i] /= PDBVELSCALINGFACTOR;
      report << plain <<"PDB velocities scaled."<<endr;
    }
  }
  else if(config.valid(InputTemperature::keyword)){    
    velocities.resize(positions.size());
    report << plain << "Using temperature "<<config[InputTemperature::keyword]<<"K for the velocities  ("<<velocities.size()<<")."<<endr;
    // Create velocities later, we need the topology for that ...
  }
  else {
    report << error << "Neither temperature nor velocity file specified."<<endr;
  }

  //
  // PSF
  //
  PSFReader psfReader;
  if(!psfReader.open(config[InputPSF::keyword]))
    report << error << "Can't open PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;
  PSF psf;
  if(!(psfReader >> psf))
    report << error << "Could not parse PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;
  report << plain << "Using PSF file \'"<<config[InputPSF::keyword]<<"\' ("<<psf.atoms.size()<<")." << endr;

  //
  // PAR
  //
  PARReader parReader;
  if(!parReader.open(config[InputPAR::keyword]))
      report << error << "Can't open PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;
  PAR par;
  if(!(parReader >> par))
    report << error << "Could not parse PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;
  report << plain << "Using PAR file \'"<<config[InputPAR::keyword]<<"\', "<<(parReader.getCharmmTypeDetected() != PAR::CHARMM28?"old":"new")<< " charmm force field.";
  if(!config[InputDihedralMultPSF::keyword].valid())
    config[InputDihedralMultPSF::keyword] = (parReader.getCharmmTypeDetected() != PAR::CHARMM28);
  if(config[InputDihedralMultPSF::keyword])
    report << " Dihedral multiplictity defined by PSF.";
  report << endr;
      
  //
  // Test input
  //
  if(positions.size() != velocities.size() || 
     positions.size() != psf.atoms.size())
    report << error << "Positions, velocities and PSF input have different number of atoms."<<endr;

  //
  // Create topology
  //
  GenericTopology* topo = TopologyFactory::make(errMsg,&config);
  if(topo == NULL){
    // Try to get some defaults with the postions known ...
    const GenericTopology* prototype = TopologyFactory::find(config.get(GenericTopology::getKeyword()).getString());
    if(prototype != NULL){
      vector<Parameter> parameters = prototype->getDefaults(positions);
      for(unsigned int i=0;i<parameters.size();++i){      
	if(!config.valid(parameters[i].keyword) && parameters[i].value.valid()){
	  config.set(parameters[i].keyword,parameters[i].value);
	  report << hint << parameters[i].keyword <<" undefined, using "<<parameters[i].value.getString()<<"."<<endr;
	}
      }
      topo = TopologyFactory::make(errMsg,&config);
    }
    if(topo == NULL)
      report << error << errMsg <<endr;
  }

  // Build topology
  buildTopology(topo,psf,par,config[InputDihedralMultPSF::keyword]);

  // Fix velocities
  if(!config.valid(InputVelocities::keyword)){
    randomVelocity(config[InputTemperature::keyword],topo,&velocities,config[InputSeed::keyword]);
    report << plain << "Random temperature : "<<temperature(topo,&velocities)<<"K"<<endr;
  }
  if((int)config[InputRemoveLinearMomentum::keyword] >= 0){
    report << plain << "Removed linear momentum: "
	   <<removeLinearMomentum(&velocities, topo)*Constant::INV_TIMEFACTOR<<endr;
  }
  else {
    report << plain << "Linear momentum : "
	   <<linearMomentum(&velocities, topo)*Constant::INV_TIMEFACTOR<<endr;
  }
  if((int)config[InputRemoveAngularMomentum::keyword] >= 0){
    report << plain << "Removed angular momentum : "
	   <<removeAngularMomentum(&positions,&velocities, topo)*Constant::INV_TIMEFACTOR<<endr;
  }
  else {
    report << plain << "Angular momentum : "
	   << angularMomentum(&positions, &velocities, topo)*Constant::INV_TIMEFACTOR<<endr;
  }
  report << plain << "Actual start temperature : "<<temperature(topo,&velocities)<<"K"<<endr;

 
  //
  // Create integrators and forces
  // 
  registerForceExemplars(topo);
 
  Integrator*  integrator = IntegratorFactory::make(errMsg,config[InputIntegrator::keyword]);
  if(integrator == NULL)
    report << error << errMsg <<endr;

  //
  // Add external modifiers
  //
  if((int)config[InputRemoveLinearMomentum::keyword] > 0){
    integrator->bottom()->adoptPreDriftOrNextModifier(new ModifierRemoveLinearMomentum(config[InputRemoveLinearMomentum::keyword]));
    report << plain << "Removing linear momentum with STS frequency "<<config[InputRemoveLinearMomentum::keyword]<<"."<<endr;
  }
  if((int)config[InputRemoveAngularMomentum::keyword] > 0){
    integrator->bottom()->adoptPreDriftOrNextModifier(new ModifierRemoveAngularMomentum(config[InputRemoveAngularMomentum::keyword]));
    report << plain << "Removing angular momentum with STS frequency "<<config[InputRemoveLinearMomentum::keyword]<<"."<<endr;
  }

  //
  // Create outputs
  // 
  OutputCollection* outputs = NULL;
  if(Parallel::iAmMaster() && (bool)config[InputOutput::keyword]){
    outputs = OutputFactory::makeCollection(errMsg,&config);
    if(outputs == NULL || !errMsg.empty())
      report << error << errMsg <<endr;
  }
  else {
    outputs = new OutputCollection;
  }

  // Initialize
  report << plain  << PROTOMOL_HR << endr;
  ScalarStructure scalar;
  scalar.molecularVirial(config[InputMolVirialCalc::keyword]);
  scalar.virial(config[InputVirialCalc::keyword]);
  report << plain << "Virial tensor : "<< scalar.virial()<<endr;
  report << plain << "Molecular virial tensor : "<< scalar.molecularVirial()<<endr;
  topo->time = (Real)config[InputFirststep::keyword]*integrator->getTimestep();
  integrator->initialize(topo,&positions,&velocities,&scalar);
  outputs->initialize(&config,integrator,topo,&positions,&velocities,&scalar);
  outputs->addToCache(pdbAtoms);
  outputs->addToCache(psf);
  outputs->addToCache(par);




  // Print
  // Parallel
  if(Parallel::isMPI)
    report << plain  << PROTOMOL_HR << "\n"
	   << "Using MPI with "<<Parallel::getNum()<<" node(s), "
           <<Parallel::getNum()-Parallel::getAvailableNum()<<" master(s), "
           <<Parallel::getAvailableNum()<<" slave(s).\n"
           <<"Distribution: "<<Parallel::getMode().getString()
           <<".\nPipe size: "<<Parallel::getPipeSize()
           <<".\nUsing barrier: "<<(Parallel::isBarrier()?"yes":"no")
           <<".\nMax number of packages per node per force: "<<(Parallel::getMaxPackages()<1?"inf":toString(Parallel::getMaxPackages()))
	   <<endr;
  // Output
  report << plain  << PROTOMOL_HR << endr;
  for(OutputCollection::const_iterator itr=const_cast<const OutputCollection*>(outputs)->begin();itr!=const_cast<const OutputCollection*>(outputs)->end();itr++){
    report << plain << "Output "<<(*itr)->getId();
    vector<Parameter> parameters = (*itr)->getParameters();
    for(unsigned int i=0;i<parameters.size();i++)
      report << " "<<parameters[i].value.getString();
    report <<"."<<endr;
  }
  if(!((bool)config[InputOutput::keyword]))
    report << "All output suppressed!" <<endr;

  // Integrator
  report << plain  << PROTOMOL_HR << endr;
  vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
  report << plain << InputIntegrator::keyword << " {"<<endl;
  for(int i=inter.size()-1;i>=0;i--){
    report << Constant::PRINTINDENT<<"Level "<<i<<" "<<inter[i].print() <<endl;
  }
  report << "}"<<endr;
  // Topology
  report << plain  << PROTOMOL_HR << endr;
  report << topo->print(&positions) <<endr;

  report << debug(10) << PROTOMOL_HR << "\nConfiguration:\n" << config.print() <<endr;
  report << debug(10) << PROTOMOL_HR << "\nTopology:\n" << TopologyFactory::print() <<endr;
  report << debug(10) << PROTOMOL_HR << "\nIntegrator:\n" << IntegratorFactory::print() <<endr;
  report << debug(10) << PROTOMOL_HR << "\nForce:\n" << ForceFactory::print() <<endr;
  report << debug(10) << PROTOMOL_HR << "\nOutput:\n" << OutputFactory::print() <<endr;


  //
  // Clear all factories
  //
  TopologyFactory::unregisterAllExemplars();
  IntegratorFactory::unregisterAllExemplars();
  ForceFactory::unregisterAllExemplars();
  OutputFactory::unregisterAllExemplars();

  // Run
  report << plain  << PROTOMOL_HR << endr;
  int step = config[InputFirststep::keyword];
  int last = step+(int)config[InputNumsteps::keyword];

  Parallel::sync();
  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
  TimerStatistic::timer[TimerStatistic::RUN].start();


  while(step < last) {
    outputs->run(step);
    int inc = outputs->getNext()-step;
    inc = std::min(last,step+inc)-step;
    step += inc;
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();
    integrator->run(inc);
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();
  }
  outputs->finalize(last);

  TimerStatistic::timer[TimerStatistic::RUN].stop();

  // Clean up
  delete topo;
  delete integrator;
  delete outputs;

  TimerStatistic::timer[TimerStatistic::WALL].stop();
  
  report << allnodesserial << plain <<"Timing" 
	 << (Parallel::isParallel()? string(" ("+toString(Parallel::getId())+")"): string(""))
	 <<" : "<<TimerStatistic()<<"."<<endr;

  Parallel::finalize();
  return 0;
}