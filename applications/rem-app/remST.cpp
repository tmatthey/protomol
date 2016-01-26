#include "REMConfigurationReader.h" // Changed from #include "ConfigurationReader.h"
#include "InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"
#include "TempFileReader.h" // Added this so that we can read in temperatures.
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
#include "ModifierShake.h"
#include "ModifierRattle.h"

#include "buildTopology.h"
#include "pmconstants.h"
#include "inputValueDefinitions.h"
#include "mathutilities.h"
#include "parseCommandLine.h"
#include "protomol.h"
#include "stringutilities.h"
#include "systemutilities.h"
#include "topologyutilities.h"

#include "ForceGroup.h"


#include <unistd.h> //next four files are for file manipulation (system calls)
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>


using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
using std::endl;
using std::string;
using std::vector;

/************************ NEW REM FUNCTION DEFS ****************************/

static string changeTempDir(string temperature, Configuration &config);
static void cleanUpFiles(string origDir, Configuration &config);
static vector< vector < int > > calcPatterns(int numReplicas);
static void calcExchMat(vector< vector< Real > > &exchMatrix, Real *aTemps, Real *aPotE);
static int calcMT(vector< vector < int > > patternVec, vector< vector< Real > > exchMatrix);
static void calcSwitchPairs(int MTPat, vector< vector < int > > patternVec, vector< vector< Real > > exchMatrix, Real *sPairs);
static void doSwitch(Vector3DBlock *positions, Vector3DBlock *velocities, Real myTemp,
                     int address, Real *historyArray, int histSize, vector<Real> &countExch);
static Real calcDelta(Real myTemp, Real otherTemp, Real myEnergy, Real otherEnergy);


/************************* REM KEYWORDS *************************************/

declareInputValue(InputREMTemps, STRING, NOCONSTRAINTS);
declareInputValue(InputREMTempDir, STRING, NOCONSTRAINTS);
declareInputValue(InputNumSwitches, INT, NOCONSTRAINTS);
declareInputValue(InputEqSwitches, INT, NOCONSTRAINTS);

defineInputValue(InputREMTemps, "REMTemperatureFile");
defineInputValue(InputREMTempDir, "TempSpace");
defineInputValue(InputNumSwitches, "NumSwitches");
defineInputValue(InputEqSwitches, "EqSwitches");

//_____________________________________________________________________ rem


int main(int argc, char **argv) {


  // Redirect all output to a file
  //std::ofstream out("out.txt");
  //report.setStream(&(out));

  // Redirect all output to std::cout
  //report.setStream(&(std::cout));

  Parallel::init(argc,argv);
  Parallel::setBarrier(true);  // Added so that we can synchronize by default
  Parallel::isolateNode(); // Isolate the node until it needs to communicate
  TimerStatistic::timer[TimerStatistic::WALL].start();

  // Changed splash screen so that it isn't protomol, but is rem.
  report << plain <<endl
         <<"#####                 "<<endl
         <<"#    #                "<<endl
         <<"#    #  ######  ##   ##"<<endl
         <<"#####   #       # # # #"<<endl
         <<"#  #    ######  #  #  #"<<endl
         <<"#   #   #       #     #"<<endl
         <<"#    #  ######  #     #"<<endl
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
  InputREMTempDir::registerConfiguration(&config); // next 4 lines register REM specific keywords
  InputNumSwitches::registerConfiguration(&config);//
  InputREMTemps::registerConfiguration(&config); //
  InputEqSwitches::registerConfiguration(&config);// 
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
  InputShake::registerConfiguration(&config,false);
  InputShakeEpsilon::registerConfiguration(&config,1e-5);
  InputShakeMaxIter::registerConfiguration(&config,30);
  InputRattle::registerConfiguration(&config,false);
  InputRattleEpsilon::registerConfiguration(&config,1e-5);
  InputRattleMaxIter::registerConfiguration(&config,30);
  InputReducedImage::registerConfiguration(&config,false);
  InputMinimalImage::registerConfiguration(&config,true);

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
  
  /********************* BEGIN REM MODIFICATIONS ***************************/
  //  We require that a temperature set file be specified on the command line so that the reader can properly
  //  substitute the temperature where required.  Hence, we read in the temperature first, before we read in
  //  config file.
  config.set(args); // Need the cmd line args to use the config[] operator...
  Parallel::integrateNode(); // Turn on MPI since we need our node id
  TempFileReader TempReader(config[InputREMTemps::keyword]); // get the temp. filename & create our reader
  TempReader.setNodeNum(Parallel::getId()); // tell the reader which temp to grab (based on node number)
  Real temp = 0; // Holds the temperature that we will be simulating at
  //Real temp, tempLow, tempHigh = 0;
  string tempString;  // holds the string version of the temperature
  TempReader >> temp; // Read it in
  //tempLow = temp * 0.9;
  //tempHigh = temp * 1.1;
  int a, b;
  tempString = fcvt(temp, 0, &a, &b); // convert it to a string (used later on)
  Parallel::isolateNode(); // Done with MPI
  /************************ END REM MODIFICATIONS ******************************/
  
  if(config.set(InputConfig::keyword,args)){
    // Read config file
    REMConfigurationReader configReader(temp); // changed from ConfigurationReader configReader;
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
  
  string originalDir = changeTempDir(tempString, config); // cd() to scratch space (part of it from above) ->REM SPECIFIC<- (improves i/o time)

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
    Parallel::integrateNode(); // Begin REM specific modifications: print out the temps in order
    report << allnodesserial << "Using temperature "<<config[InputTemperature::keyword]<<"K for the velocities  ("<<velocities.size()<<")."<<endr;
    Parallel::isolateNode(); //end mods.
    // Create velocities later, we need the topology for that ...
  }
  else {
    // since we create a random temp distribution, just resize velocity array and go on.
    velocities.resize(positions.size());
    //  report << error << "Neither temperature nor velocity file specified."<<endr;
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

  topo->minimalMolecularDistances = topo->checkMoleculePairDistances(positions);
  if((bool)config[InputReducedImage::keyword] && !topo->minimalMolecularDistances){
    Vector3DBlock tmp(positions);
    topo->minimalImage(tmp);
    if(topo->checkMoleculePairDistances(tmp)){
      positions = tmp;
      report << plain <<"Fixed minimal molecule distances."<<endr;
      topo->minimalMolecularDistances = true;
    }
    else {
      report << plain <<"Could not fixed minimal molecule distances."<<endr;      
      topo->minimalMolecularDistances = false;
    }
      
  }
  
  // Fix velocities
  if(!config.valid(InputVelocities::keyword)){
    randomVelocity(temp,topo,&velocities,config[InputSeed::keyword]); // changed temp from configuration structure keyword
    Parallel::integrateNode();
    report << allnodesserial << "Random temperature : "<<temperature(topo,&velocities)<<"K"<<endr;
    Parallel::isolateNode();
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
  Parallel::integrateNode();
  report << allnodesserial << "Actual start temperature : "<<temperature(topo,&velocities)<<"K"<<endr;
  Parallel::isolateNode();

 
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
    report << plain << "Removing angular momentum with STS frequency "<<config[InputRemoveAngularMomentum::keyword]<<"."<<endr;
  }
  if((bool)config[InputShake::keyword] && 
     (Real)config[InputShakeEpsilon::keyword] > 0.0 && 
     (int)config[InputShakeMaxIter::keyword] > 0) {
    // initialize the SHAKE modifier 
    integrator->bottom()->adoptPostDriftOrNextModifier(integrator->bottom()->createShakeModifier((Real)config[InputShakeEpsilon::keyword],
                                                                                                 (int)config[InputShakeMaxIter::keyword]));
    report << plain << "Shake with epsilon "<< config[InputShakeEpsilon::keyword] <<", max "<<(int)config[InputShakeMaxIter::keyword]<<" iteration(s)."<<endr;
  }   

  // find out if the user also wants to use RATTLE
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
  topo->time = (Real)config[InputFirststep::keyword]*integrator->getTimestep();
  integrator->initialize(topo,&positions,&velocities,&scalar);
  outputs->initialize(&config,integrator,topo,&positions,&velocities,&scalar);
  outputs->addToCache(pdbAtoms);
  outputs->addToCache(psf);
  outputs->addToCache(par);




  // Print
  // Parallel
  //Parallel::setMode(ParallelType::STATIC);
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
  
  /***************************** BEGIN REM ENGINE ******************************/
  Parallel::integrateNode();
  int NumReplicas = Parallel::getNum();
  Parallel::isolateNode();
  //report << "Number of replicas: " << NumReplicas << endr;
  int numTotalSwitches = (int)config[InputNumSwitches::keyword]; // number of REM switches to attempt
  int remStepLength = (int)config[InputNumsteps::keyword]; // number of steps per switch
  int EqSwitch = (int)config[InputEqSwitches::keyword]; //number of switches to skip to allow equilibration
  vector < Real > countExchanges;
  countExchanges.resize(NumReplicas);
  for (int i = 0; i < NumReplicas; i++)
  {
      countExchanges[i] = 0.0;
  }
  Real *switchPairs = new Real[NumReplicas];
  Real myTemp = 0;
  Real myEnergy = 0;
  int stepLength = 0;
  Real numSwitches = 0; //keeps track of the number of switches performed

  // Variables Specifically for the MasterNode
  int MTPattern = 0;
  vector< vector < int > > patVec;
  Parallel::integrateNode(); // to getId info
  if (Parallel::getId() == Parallel::getMasterId()) // Master Node calculates all patterns for MTREM
  {
    patVec = calcPatterns(NumReplicas);
  }
  Parallel::isolateNode();
  vector< vector < Real > > exchMat;
  exchMat.resize(NumReplicas);
  for (int i = 0; i < NumReplicas; i++)
    exchMat[i].resize(NumReplicas);
  //report << "exchMatrix of size:" << exchMat.size() << "x" << exchMat[0].size() << endr;
  Real *allTemps = new Real[NumReplicas];
  for (int k = 0; k < NumReplicas; k++) // initialize array
    allTemps[k] = 0;
  Real *allPotE = new Real[NumReplicas];
  for (int k = 0; k < NumReplicas; k++) // initialize array
    allPotE[k] = 0;
  // End Master Node Variables 

  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
  TimerStatistic::timer[TimerStatistic::RUN].start();

  Real *tempHistory = 0; // holds the historical data for this particular configuration.
  tempHistory = new Real[numTotalSwitches + 1]; // the first data member is the "originating address", the rest, historical data.
  if (tempHistory == 0) // check array
    report << error << "Can't create history array.  Exiting." << endr;
  for (int i = 0; i < numTotalSwitches; i++) // initialize array
    tempHistory[i] = 0;

  Parallel::integrateNode(); // turn on MPI
  tempHistory[0] = Parallel::getId();

  // output the intial conditions
  outputs->run(0);
  // begin switching/simulating
  for (int i = 0; i < numTotalSwitches; i++) {

    Parallel::sync();
    tempHistory[i + 1] = temp; // record the fact that this configuration was at this temperature
    //std::cout << "Replica " << Parallel::getId() << " is simulating with temperature " 
    //       << temperature(topo, &velocities) << " round " << i << endl;
    stepLength += remStepLength;

    // Do an MD Run.
    Parallel::isolateNode(); // turn off MPI for the integrator, else we get interesting results...
    while(step < stepLength) { // changed from step < last, end addition
      // moved "outputs->run(step);" as to make sure simulation segment gets logged prior to switch
      int inc = outputs->getNext()-step;
      inc = std::min(stepLength,step+inc)-step;
      step += inc;
      TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();
      integrator->run(inc);
      TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();
      if ( (step) != (numTotalSwitches * remStepLength))
	outputs->run(step);
    }

    if (i >= EqSwitch){
	// Turn on MPI for the switch
	Parallel::integrateNode(); // begin additions
	TimerStatistic::timer[TimerStatistic::IDLE].start();
	Parallel::sync(); // sync up after integration
	TimerStatistic::timer[TimerStatistic::IDLE].stop();

	myTemp = temp; // temperature(topo,&velocities);
	myEnergy = scalar.potentialEnergy();

	//std::cout << "Replica " << Parallel::getId() << " Is done simulating.  Ending temperature: " 
	//	      << myTemp << "  Ending PotE: " << myEnergy << endl;

	// get arrays of temp and PotE info - NOTE GATHER NUM IS ONE FOR SEND AND RECIVE -- NOT TOTAL REC SIZE!
	Parallel::gather(&myTemp,1,allTemps,Parallel::getMasterId());
	Parallel::sync();
	Parallel::gather(&myEnergy,1,allPotE,Parallel::getMasterId());
	Parallel::sync();

	// Have master node report Temps and PotE
	//for (int k = 0; k < NumReplicas; k++)
	//  report << "Temperature " << k << ": " << allTemps[k] << "   PotE: " << allPotE[k] << endl;

	if  (Parallel::getId() == Parallel::getMasterId()) // MasterNode
	{
	  calcExchMat(exchMat, allTemps, allPotE);

	  if ((i % 2) == 0)
	  {
	    //MTPattern = calcMT(patVec,exchMat); 
	    MTPattern = i % patVec.size();
            std::cout << "MTPattern: " << MTPattern << endl;
	  }
	  else
	  {
	    MTPattern = MTPattern + 1;
	    std::cout << "MTPattern: " << MTPattern << endl;
	  }
	  calcSwitchPairs(MTPattern, patVec, exchMat, switchPairs);
	}

	Parallel::sync();
	Parallel::bcastSlaves(switchPairs, switchPairs + NumReplicas);

	//report << allnodes << plain << "The switch pairs: ";
	//for (int k=0; k < NumReplicas; k++)
	//  report << allnodes << plain << switchPairs[k] << " ";
	//report << allnodes << plain << endr;

	Parallel::sync();

	// doSwitch if Replica is part of switch set
	for (int i=0; i < NumReplicas; i++)
	{
	  if (Parallel::getId() == switchPairs[i])
	  {
	    if ((i % 2) == 0)
	    {
	      //report << allnodes << plain << "Switch Do - Replicas " << switchPairs[i] << " and " << switchPairs[i+1] << endl;
	      doSwitch(&positions, &velocities, myTemp, int(switchPairs[i+1]), tempHistory, numTotalSwitches + 1,countExchanges);
	      //report << allnodes << plain << "Switch Complete - Replicas " << switchPairs[i] << " and " << switchPairs[i+1] << endl;
	      Parallel::isolateNode(); // hide MPI to reinitialize the stats (e.g. potential energy, etc.)
	      integrator->initialize(topo, &positions, &velocities, &scalar);
	      Parallel::integrateNode();
	    }
	    else
	    {
	      //report << allnodes << plain << "Switch Do - Replicas " << switchPairs[i] << " and " << switchPairs[i-1] << endl;
	      doSwitch(&positions, &velocities, myTemp, int(switchPairs[i-1]), tempHistory, numTotalSwitches + 1,countExchanges);
	      //report << allnodes << plain << "Switch Complete - Replicas " << switchPairs[i] << " and " << switchPairs[i-1] << endl;
	      Parallel::isolateNode(); // hide MPI to reinitialize the stats (e.g. potential energy, etc.)
	      integrator->initialize(topo, &positions, &velocities, &scalar);
	      Parallel::integrateNode();
	      numSwitches++; // used to monitor exchange rate.  Updates counter for higher temp
	    }
	  }
	}
    }
  }

  // Report Acceptance Ratio Data
  Real *numAccepts = 0;
  if (Parallel::getId() == 0) { // data structure to hold acceptance ratio data
    numAccepts = new Real[NumReplicas];
    if (numAccepts == 0)
      report << error << "Can't create numAccepts!" << endr;    
  }

  // get the addressing information for replica histories
  Real *history = new Real[NumReplicas];
  if (history == 0)
    report << error << "Cant create address array.  Quitting!" << endr;
  Parallel::allgather(tempHistory, 1, history);  // get address data from all nodes (first element that we saved in the array.
  // the processor we're receiving from is the location of our ID in the array
  // e.g. if rank = 2, then with an array of [5 1 4 2 0 3], we receive from rank 3, and send to rank 5.
  int recv_address = 0;
  for (int i = 0; i < NumReplicas; i++)
    if (history[i] == Parallel::getId())
      recv_address = i;
  // actual send/recv call
  Parallel::sendrecv_replace(tempHistory, numTotalSwitches + 1, (int) history[Parallel::getId()], recv_address);
  // replace address space with the size of the array
  tempHistory[0] = numTotalSwitches + 1;
  outputs->addToCache<Real *>(tempHistory);

  TimerStatistic::timer[TimerStatistic::COMMUNICATION].start();
  Parallel::gather(&numSwitches, 1, numAccepts, 0); // gather the number of switches to the root node (node 0)
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].stop();
  if (Parallel::getId() == 0) {
    vector<Real> Accepts;
    for (int i = 0; i < NumReplicas; i++)
      Accepts.push_back(numAccepts[i] / (numTotalSwitches - EqSwitch)); // make into ratios
    outputs->addToCache<vector<Real> >(Accepts); // add this to the cache so the output writer can get to it.
  }
  
  std::cout << "Exchange Counts for Rep " << Parallel::getId() <<": ";
  for (int i = 0; i < Parallel::getNum(); i++)
    std::cout << countExchanges[i] << " ";
  std::cout << endl;

  outputs->finalize(numTotalSwitches * remStepLength); // changed from finalize(last)

  TimerStatistic::timer[TimerStatistic::RUN].stop();

  // Clean up
  delete topo;
  delete integrator;
  delete outputs;

  TimerStatistic::timer[TimerStatistic::WALL].stop();
  
  report << allnodesserial << plain <<"Timing" 
         << (Parallel::isParallel()? string(" ("+toString(Parallel::getId())+")"): string(""))
         <<" : "<<TimerStatistic()<<"."<<endr;

  Parallel::sync();
  cleanUpFiles(originalDir, config); // copy everything back!
  Parallel::finalize();

  return 0; // yippee! we're done! we can go home!
}

// _______________________________________________ BEGIN NEW REM FUNCTIONS

// _______________________________________________changeTempDir
string changeTempDir(string temperature, Configuration &config) {
  // Change to a temp dir.  (performance)
  char initDir[200];
  char *result = getcwd(initDir, 200);
  if (result == 0)
    report << error << "Can't save current dir.  Quitting!" << endr;
  
  string new_dir = (std::string)config[InputREMTempDir::keyword];

  // added this to remove files from previosly failed jobs
  //string sys_string = "/bin/rm -rf ";
  //sys_string += new_dir;
  //system(sys_string.c_str());

  int retval = mkdir(new_dir.c_str(), 0700);
  new_dir += temperature;
  new_dir += "k"; // adds k label to all temperature files for post processing
  retval = mkdir(new_dir.c_str(), 0700);
  if (retval == -1)
    report << error << "Problem making temp dir.  Maybe the disk is full?" << endr;
  string sys_command = "/bin/cp * ";
  sys_command += new_dir;
  system(sys_command.c_str());
  retval = chdir (new_dir.c_str());
  if (retval == -1)
    report << error << "Can't chdir() to the sim dir.  Please check your system." << endr;
  return string(result);
}

// _______________________________________________cleanUpFiles
void cleanUpFiles(string origDir, Configuration &config) {
  char curDir[200];
  if (getcwd(curDir, 200) == 0)
    report << error << "Can't get CWD.  Quitting!" << endr;

  chdir("/");
  string sys_string = "/bin/cp -r ";
  sys_string += curDir;
  sys_string += " ";
  sys_string += origDir;
  system(sys_string.c_str());
  sys_string = "/bin/rm -rf ";
  sys_string += curDir;
  sys_string += "; /bin/rmdir ";
  sys_string += string(config[InputREMTempDir::keyword]);
  system(sys_string.c_str());
}

//________________________________________________calcPatterns

vector< vector < int > > calcPatterns(int numReplicas)
{   
    int first = 0;
    int second = 0;
    int temp = 0;
    int M = numReplicas;
    vector< vector < int > > patVec;
    patVec.resize(2*(M-1));
    //cout << patVec.size() << endl;
    for (unsigned int w=0; w < patVec.size(); w++)
      patVec[w].resize(0);
    
    for (int i=0; i < (M-1); i++)
    {
      first = 2*i;
      second = 2*i+1;
      
      patVec[(first)].push_back(0);    
      patVec[(first)].push_back((M-1)-i);
      
      for (int j=1; j<M; j+=2)
      {
        if (j+1 != M-i)
        {
          patVec[second].push_back(j-1);
          patVec[second].push_back(j);
        }else
        {
          temp = first;
          first = second;
          second = temp;
        }
        if ((j+2 != M-i) && (j+2 <= M))
        {
          patVec[first].push_back(j);
          patVec[first].push_back(j+1);
        }
      }
    }
    
    for (unsigned int k=0; k < patVec.size(); k++)
    {
      report << allnodes << plain << "Pattern " << k << endr;
      for (unsigned int l=0; l < patVec[k].size(); l++)
        report << allnodes << plain << patVec[k][l] << " ";
      report << allnodes << plain << endr;
    }
    
    return patVec;
}
 
// _______________________________________________calcExchMat
void calcExchMat(vector< vector< Real > > &exchMatrix, Real *aTemps, Real *aPotE)
  {
    for (int i=0; i < Parallel::getNum(); i++)
    {
      //report << "Entering calcExchMat Function at 3" << "-- i value is: " << i << endr;
      for (int j=0; j < Parallel::getNum(); j++)
      {
        //report << "About to intitialize Matrix of size:" << exchMatrix.size() << "x" << exchMatrix[i].size() 
	//       <<" at i=" << i << "  j=" << j << endr;
        exchMatrix[i][j]=0;
	// report << "Just intialized Matrix at i=" << i << "  j=" << j << endr;
        if ( ((i==0) && (j!=0)) || (i+1 == j) )
        {
          // print out who you are attempting to switch with
          // std::cout << "Replica " << i << ": checking switch accept with replica "
          //          << j << "." << endr;
          Real delta = calcDelta(aTemps[i], aTemps[j], aPotE[i], aPotE[j]); // calculate delta
          report << allnodes << plain << "Replicas " << i << " and " << j << ": delta = " << delta << endr;
          Real w;
          // Metropolis criterion
          if (delta <= 0)
            w = 1.0;
          else
            w = exp(-delta);
          double num = drand48();
          if (num <= w)
          {
            exchMatrix[i][j] = 1;
            //report << "Replicas " << i << " and "
            //       << j << ": Switch accept succeeded, returning true" << endr;
          }else
          {
            //report << "Replicas " << i << " and "
            //       << j << " Switch accept failed, returning false" << endr;
          }
        }
      }
    }
    //report << "Leaving calcExchMat Function" << endr;
    //report << "Matrix of size:" << exchMatrix.size() << "x" << exchMatrix[1].size() << endr;
    for (int i=0; i < Parallel::getNum(); i++)
      {
	for (int j=0; j < Parallel::getNum(); j++)
	  report << allnodes << plain << exchMatrix[i][j] << " ";
	report << allnodes << plain << endr;
      }
  }

// _______________________________________________calcMT
int calcMT(vector< vector < int > > patternVec, vector< vector< Real > > exchMatrix)
  {
    Real MaxTraversal = 0;
    unsigned int MTPattern = 0;
    for (unsigned int i=0; i < patternVec.size(); i+=2)
    {
      Real CurTraversal = 0;
      for (unsigned int j=0; j < patternVec[i].size(); j+=2)
      {
        if (exchMatrix[patternVec[i][j]][patternVec[i][j+1]])
        {
          CurTraversal = CurTraversal + (patternVec[i][j+1] - patternVec[i][j]);
        }
      }
      if (CurTraversal >= MaxTraversal)
      {
        MaxTraversal = CurTraversal;
        MTPattern = i;
      }
    }
    report << allnodes << plain << "MaxTraversal: " << MaxTraversal << endr;
    return MTPattern;
  }
 
//_________________________________________________ calcSwitchPairs
void calcSwitchPairs(int MTPat, vector< vector < int > > patternVec, vector< vector< Real > > exchMatrix, Real *sPairs) {
  int i = 0;

  // re-initialize array
  for (int k=0; k < Parallel::getNum(); k++)
    sPairs[k]= -2;

  for (unsigned int j=0; j < patternVec[MTPat].size(); j+=2)
  {
    if (exchMatrix[patternVec[MTPat][j]][patternVec[MTPat][j+1]])
    {
      sPairs[i] = patternVec[MTPat][j];
      sPairs[i+1] = patternVec[MTPat][j+1];      
      i = i + 2;
    }
  }

}
 
//_________________________________________________ doSwitch
void doSwitch(Vector3DBlock *positions, Vector3DBlock *velocities, Real myTemp, int address, Real *historyArray, int histSize, vector<Real> &countExch ) {
  Real otherTemp = 0;
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].start();
  // swap temps
  //report << allnodes << plain << "Replica " << Parallel::getId() << " swap temps" << endl;
  Parallel::sendrecv(&myTemp, 1, address, &otherTemp, 1, address);
  // swap positions
  //report << allnodes << plain << "Replica " << Parallel::getId() << "swap positions" << endl;
  Parallel::sendrecv_replace(positions, address, address);
  //swap velocities
  //report << allnodes << plain << "Replica " << Parallel::getId() << "swap velocities" << endl;
  Parallel::sendrecv_replace(velocities, address, address);
  // send temp history
  //report << allnodes << plain << "Replica " << Parallel::getId() << "swap temp history" << endl;
  Parallel::sendrecv_replace(historyArray, histSize, address, address);
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].stop();
  // rescale factor
  //report << allnodes << plain << "Replica " << Parallel::getId() << "rescale temp" << endl;
  Real factor = sqrt(myTemp / otherTemp);
  for (unsigned int j = 0; j < velocities->size(); j++)
    (*velocities)[j] = (*velocities)[j] * factor; // rescale velocity
  countExch[address] = countExch[address] + 1;  
}

//_________________________________________________ calcDelta
Real calcDelta(Real myTemp, Real otherTemp, Real myEnergy, Real otherEnergy) {
  Real result = 0;
  Real Kb = Constant::BOLTZMANN;
  result = (1 / Kb / myTemp - 1 / Kb / otherTemp) * (otherEnergy - myEnergy);
  return result;
}
