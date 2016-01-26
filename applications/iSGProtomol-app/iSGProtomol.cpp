#include "ConfigurationReader.h"
#include "InputPosVel.h"
#include "iSGPARReader.h"
#include "iSGIntegrator.h"
#include "PSFReader.h"
#include "TRANSReader.h"
#include "TimerStatistic.h"

#include "TopologyFactory.h"
#include "ForceFactory.h" 
#include "IntegratorFactory.h"
#include "OutputFactory.h"  

#include "iSGregisterForceExemplars.h"
#include "iSGregisterIntegratorExemplars.h"
#include "registerTopologyExemplars.h"
#include "iSGregisterOutputExemplars.h"

#include "OutputCollection.h" 

#include "Parameter.h"
#include "ScalarStructure.h"
#include "Parallel.h"

#include "ModifierRemoveLinearMomentum.h"
#include "ModifierRemoveAngularMomentum.h"

#include "GenericTopology.h"
#include "buildISGTopology.h"
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

// adding some more keywords, which do not relate to ordinary MD
declareInputValue(InputNumComp,INT,NOTNEGATIVE);
declareInputValue(InputXSC,STRING,NOTEMPTY);
declareInputValue(InputTRANS,STRING,NOTEMPTY);
declareInputValue(InputReportFreq,BOOL,NOCONSTRAINTS);
declareInputValue(InputTICalc,BOOL,NOCONSTRAINTS);

defineInputValue(InputNumComp,"components");
defineInputValue(InputXSC,"xscfile");
defineInputValue(InputTRANS,"transfile");
defineInputValue(InputReportFreq,"reportOnFreq");
defineInputValue(InputTICalc,"TICalc");



//_____________________________________________________________________ protomol for iSGMD
//  iSGMD -- Isomolar SemiGrand ensemble Molecular Dynamics
//  see T.I. Morrow and E.J. Maginn, "Isomolar semigrand ensemble molecular dynamics:
//  Development and application to liquid-liquid equilibria" J. Chem. Phys. 122, 054504 (2005).

int main(int argc, char **argv) {

  // Redirect all output to a file
  //std::ofstream out("out.txt");
  //report.setStream(&(out));

  // Redirect all output to std::cout
  //report.setStream(&(std::cout));

  Parallel::init(argc,argv);

  TimerStatistic::timer[TimerStatistic::WALL].start();

  report << plain <<endl
	 <<"#    #####    #####    ######                                  #     #"<<endl
	 <<"    #     #  #     #   #     #  #####    ####   #####   ####   ##   ##   ####   #"<<endl
	 <<"#   #        #         #     #  #    #  #    #    #    #    #  # # # #  #    #  #"<<endl
	 <<"#    #####   #  ####   ######   #    #  #    #    #    #    #  #  #  #  #    #  #"<<endl
	 <<"#         #  #  #  #   #        #####   #    #    #    #    #  #     #  #    #  #"<<endl
	 <<"#   #     #  #     #   #        #   #   #    #    #    #    #  #     #  #    #  #"<<endl
	 <<"#    #####    #####    #        #    #   ####     #     ####   #     #   ####   ######"<<endl
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
  InputVirialCalc::registerConfiguration(&config,true);
  InputMolVirialCalc::registerConfiguration(&config,true);
  InputTRANS::registerConfiguration(&config);
  InputNumComp::registerConfiguration(&config);
  InputXSC::registerConfiguration(&config);
  InputOutput::registerConfiguration(&config,true);
  InputShake::registerConfiguration(&config,false);
  InputShakeEpsilon::registerConfiguration(&config,1e-5);
  InputShakeMaxIter::registerConfiguration(&config,30);
  InputRattle::registerConfiguration(&config,false);
  InputRattleEpsilon::registerConfiguration(&config,1e-5);
  InputRattleMaxIter::registerConfiguration(&config,30);
  InputReportFreq::registerConfiguration(&config,false);
  InputTICalc::registerConfiguration(&config,false);
  InputReducedImage::registerConfiguration(&config,false);
  InputMinimalImage::registerConfiguration(&config,true);

  //
  // Configure OutputFactory
  //
  iSGregisterOutputExemplars();
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
  vector<vector<string> > args = parseCommandLine(argc,argv,&config,iSGregisterForceExemplars);
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
  if(!config[InputNumComp::keyword].valid())
    report << error << "# of mixture components undefined." << endr;
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

  // open the positions file 
  if(!reader.open(config[InputPositions::keyword]))
    report << error << "Can't open position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;

  // create the coordinates objects
  Vector3DBlock positions;
  vector<PDB::Atom> pdbAtoms;

  // read in the positions file information
  if(reader.tryFormat(InputPosVelType::PDB)){
    PDB pdb;
    if(!(reader >> pdb))
      report << error << "Could not parse PDB position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;
    // store the xyz coordinates
    swap(positions,pdb.coords);
    // store the atom names, etc.
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

  // read in the xyz velocities
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
  // the velocities will be initialized based upon the specified temperature
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
  // create the PSF reader object
  PSFReader psfReader;
  if(!psfReader.open(config[InputPSF::keyword]))
    report << error << "Can't open PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;

  // create the storage object for the PSF information
  PSF psf;

  // read in the PSF file, which will contain the current mass and charge of each atom
  if(!(psfReader >> psf))
    report << error << "Could not parse PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;
  report << plain << "Using PSF file \'"<<config[InputPSF::keyword]<<"\' ("<<psf.atoms.size()<<")." << endr;

  //
  // TRANS
  //
  // create the TRANS reader object
  // The TRANS file contains the details of how we are to perform all of the possible
  // molecular transformations.
  TRANSReader transReader;

  if(!transReader.open(config[InputTRANS::keyword]))
    report << error << "Can't open TRANS file \'"<<config[InputTRANS::keyword]<<"\'."<<endr;

  // create the storage object for the TRANS information
  TRANS trans;

  // read in the TRANS file information
  if(!(transReader >> trans))
    report << error << "Could not parse TRANS file \'"<<config[InputTRANS::keyword]<<"\'."<<endr;
  report << plain << "Using TRANS file \'"<<config[InputTRANS::keyword]<<"\'." << endr;

  //
  // PAR
  //
  // create the PAR reader for iSGMD simulations
  // this differs from the standard PAR reader in that it must read in
  // more than one set of force constants for each bond, angle, dihedral, etc.
  // create the reader object
  iSGPARReader parReader;
  parReader.setNumComps(config[InputNumComp::keyword]);
 
  // open the PAR file
  if(!parReader.open(config[InputPAR::keyword]))
    report << error << "Can't open PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;
  
  // create the storage object for the PAR file information
  iSGPAR par;
  
  // read the PAR file information
  if(!(parReader >> par))
    report << error << "Could not parse PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;
  report << plain << "Using PAR file \'"<<config[InputPAR::keyword]<<"\', "<<(parReader.getCharmmTypeDetected() != iSGPAR::CHARMM28?"old":"new")<< " charmm force field.";
  if(!config[InputDihedralMultPSF::keyword].valid())
    config[InputDihedralMultPSF::keyword] = (parReader.getCharmmTypeDetected() != iSGPAR::CHARMM28);
  if(config[InputDihedralMultPSF::keyword])
    report << " Dihedral multiplictity defined by PSF.";
  report << endr;
      
#ifdef SHOW_PAR_INFO
  // for debugging purposes...output the PAR file information  
  for (unsigned int i=0; i<par.bonds.size(); i++) 
    report << par.bonds[i] << endr;
  for (unsigned int i=0; i<par.angles.size(); i++) 
    report << par.angles[i] << endr;
  for (unsigned int i=0; i<par.dihedrals.size(); i++) 
    report << par.dihedrals[i] << endr;
  for (unsigned int i=0; i<par.impropers.size(); i++) 
    report << par.impropers[i] << endr;
  for (unsigned int i=0; i<par.nonbondeds.size(); i++) 
    report << par.nonbondeds[i] << endr;
  for (unsigned int i=0; i<par.nbfixs.size(); i++) 
    report << par.nbfixs[i] << endr;
  for (unsigned int i=0; i<par.hbonds.size(); i++) 
    report << par.hbonds[i] << endr;  
#endif


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
  buildISGTopology(topo,psf,par);

  // Fix velocities
  if(!config.valid(InputVelocities::keyword)){
    report << debug() << "Init temp = " << config[InputTemperature::keyword] << endr;
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
  iSGregisterForceExemplars(topo);
 
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
  // Set up SHAKE if necessary
  if((bool)config[InputShake::keyword] && 
     (Real)config[InputShakeEpsilon::keyword] > 0.0 && 
     (int)config[InputShakeMaxIter::keyword] > 0) {
    // initialize the SHAKE modifier 
    integrator->bottom()->adoptPostDriftOrNextModifier(integrator->bottom()->createShakeModifier((Real)config[InputShakeEpsilon::keyword],
                                                                                                 (int)config[InputShakeMaxIter::keyword]));
    report << plain << "Shake with epsilon "<< config[InputShakeEpsilon::keyword] <<", max "<<(int)config[InputShakeMaxIter::keyword]<<" iteration(s)."<<endr;
  }   
  // Set up RATTLE if necessary
  if((bool)config[InputRattle::keyword] &&
     (Real)config[InputRattleEpsilon::keyword] > 0.0 &&
     (int)config[InputRattleMaxIter::keyword] > 0) {
    // initialize the RATTLE modifier
    integrator->bottom()->adoptPostStepModifier(integrator->bottom()->createRattleModifier((Real)config[InputRattleEpsilon::keyword],
                                                                                           (int)config[InputRattleMaxIter::keyword]));
    report << plain << "Rattle with epsilon "<< config[InputRattleEpsilon::keyword] <<", max "<<(int)config[InputRattleMaxIter::keyword]<<" iteration(s)."<<endr;
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

  // ISG integrator
  iSGIntegrator* isgIntegrator = dynamic_cast<iSGIntegrator*>( integrator->bottom() );

  // store the detailed molecule insertion/deletion information
  isgIntegrator->indexTopology(topo, trans, psf, par, config[InputDihedralMultPSF::keyword], seed);

  // XSC (eXtended System Coordinates)
  // read this only if an XSC file is specified in the config file. 
  if(config[InputXSC::keyword].valid()) isgIntegrator->readXSCs(config[InputXSC::keyword],topo);
  integrator->initialize(topo,&positions,&velocities,&scalar);
 
  // initialize the output cache
  report << plain  << PROTOMOL_HR << endr;
  topo->time = (Real)config[InputFirststep::keyword]*integrator->getTimestep();
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
  for(OutputCollection::const_iterator itr=const_cast<const OutputCollection*>(outputs)->begin();itr!=const_cast<const OutputCollection*>(outputs)->end();++itr){
    report << plain << "Output "<<(*itr)->getId();
    vector<Parameter> parameters = (*itr)->getParameters();
    for(unsigned int i=0;i<parameters.size();i++)
      report << " "<<parameters[i].value.getString();
    report <<"."<<endr;
  }
  if(!((bool)config[InputOutput::keyword]))
    report << "All output suppressed!" <<endr;

  // output the integrator details
  report << plain  << PROTOMOL_HR << endr;
  vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
  report << plain << InputIntegrator::keyword << " {"<<endl;
  for(int i=inter.size()-1;i>=0;i--){
    report << Constant::PRINTINDENT<<"Level "<<i<<" "<<inter[i].print() <<endl;
  }
  report << "}"<<endr;

  // output the topology details
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
  int myNumSteps = (int)config[InputNumsteps::keyword];
  int last = step + myNumSteps;
  
  Parallel::sync();
  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
  TimerStatistic::timer[TimerStatistic::RUN].start();

  // loop over all timesteps 
  while(step < last) {

    // Do all output at a specified regular frequency if desired,
    // otherwise only do output whenever a transformation attempt is complete
    if ( (bool)config[InputReportFreq::keyword] ) 
      scalar.output(true);

    // Do all output defined by the configuration file
    outputs->run(step);

    // Get the next step to do output and substract to get the increment
    int inc = outputs->getNext()-step;

    // Adjust the such that inc+step <= last
    inc = std::min(last,step+inc)-step;
    step += inc;
  
    // randomly pick a molecule to be transformed
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();
    integrator->run(inc);
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();
  }
  outputs->finalize(last);

  // if we are doing a TI calculation (lambda is constant at a nonzero value)
  // output the average chemical potential difference ( d(deltaMu) / d(lambda) )
  if((bool)config[InputTICalc::keyword]) {
    report.precision(8);
    report << plain << "d(E) / d(lambda) = " << isgIntegrator->getAveDeltaMu(myNumSteps) << endr;
  }
  // if we are doing a standard iSGMD simulation then report the average chemostat temperature
  else {
    report.precision(6);
    Real AverageLambdaT = isgIntegrator->getAveLambdaT(myNumSteps);
    report << plain << "Average lambda temperature (K) = " << AverageLambdaT << endr;
  }
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
