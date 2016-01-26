//XXX #include "BVFReader.h"
#include "registerIntegratorExemplars.h"
#include "LeapfrogIntegrator.h"
#include "ScalarStructure.h"
#include "Parallel.h"
#include "parseCommandLine.h"
#include "TimerStatistic.h"

#include "registerTopologyExemplars.h"
#include "registerForceExemplars.h"

#include "TopologyFactory.h"
#include "IntegratorFactory.h"
#include "ForceFactory.h"
#include "ForceGroup.h"

#include "buildTopology.h"
#include "topologyutilities.h"

#include "DCDTrajectoryReader.h"
#include "ConfigurationReader.h"

#include "inputValueDefinitions.h"

#include "InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"
#include "PDBReader.h"
#include "systemutilities.h"

#include "protomol.h"

#include <vector>
#include <string>

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;

using std::vector;
using std::string;
using std::endl;

//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {

  Parallel::init(argc,argv);

  TimerStatistic::timer[TimerStatistic::WALL].start();

  report << plain <<endl
	 <<"              #"<<endl
	 <<"              #                                                "<<endl
	 <<"              #   ####           ####             ####   #    #"<<endl
	 <<" ## ##    #####  #    #   ###   #    #   ####    #    #  #    #"<<endl
	 <<"#  #  #  #    #  ######  #   #  ######  #    #    #####   #####"<<endl
	 <<"#  #  #  #    #  #       #   #  #       #             #       #"<<endl
	 <<"#  #  #   #####   #####  #   #   #####  #        #####   ##### "<<endl
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
  // Configuration
  //
  Configuration config;

  // register config file keywords
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
  //XXX  InputBVF::registerConfiguration(&config);

  // create the topology factory
  registerTopologyExemplars();
  TopologyFactory::registerAllExemplarsConfiguration(&config);

  //
  // Configure IntegratorFactory
  //
  registerIntegratorExemplars();

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

  // Check if configuration is complete
  //  FIXME later.
  string errMsg;
  //if(!config.validConfiguration(errMsg))
  //report << plain << endl << errMsg <<endr; 
  //if(!config[InputFirststep::keyword].valid())
  //report << error << "Firststep undefined."<<endr;
  //if(!config[InputNumsteps::keyword].valid())
  //report << error << "Numsteps undefined."<<endr;


  //
  // PSF
  //
  PSFReader psfReader;
  if(!psfReader.open(config[InputPSF::keyword]))
    report << error << "Can't open PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;

  // create a PSF object
  PSF psf;

  // store the PSF file information in psf
  if(!(psfReader >> psf))
    report << error << "Could not parse PSF file \'"<<config[InputPSF::keyword]<<"\'."<<endr;
  report << plain << "Using PSF file \'"<<config[InputPSF::keyword]<<"\' ("<<psf.atoms.size()<<")." << endr;

  //
  // PAR
  //
  PARReader parReader;
  if(!parReader.open(config[InputPAR::keyword]))
    report << error << "Can't open PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;

  // create a PAR object
  PAR par;

  // store the PAR file force parameters in par
  if(!(parReader >> par))
    report << error << "Could not parse PAR file \'"<<config[InputPAR::keyword]<<"\'."<<endr;
  report << plain << "Using PAR file \'"<<config[InputPAR::keyword]<<"\', "<<(parReader.getCharmmTypeDetected() != PAR::CHARMM28?"old":"new")<< " charmm force field.";
  if(!config[InputDihedralMultPSF::keyword].valid())
    config[InputDihedralMultPSF::keyword] = (parReader.getCharmmTypeDetected() != PAR::CHARMM28);
  if(config[InputDihedralMultPSF::keyword])
    report << " Dihedral multiplictity defined by PSF.";
  report << endr;


  // 
  // Positions.  First try to open as a dcd, then pdb, then xyz.  Exit if all
  // 3 fail.
  //

  DCDTrajectoryReader inputDCD( config[InputPositions::keyword] );

  // PDB coordinates will be stored here
  Vector3DBlock positions;

  // DCD trajectories will be stored here
  vector<Vector3DBlock> trajectory;

  // find out if this is a DCD file
  if(!inputDCD.tryFormat() ) {

    // it is not a DCD file, so create a PDB/XYZ reader object
    InputPosVel reader;

    // if this is not a PDB or XYZ file then report an error and quit
    if(!reader.open(config[InputPositions::keyword]))
      report << error << "Can't open position file \'"
	     << config[InputPositions::keyword]<<"\'."<<endr;

    if(!(reader >> positions))
      report << error << "Could not parse position file \'"
	     <<config[InputPositions::keyword]
	     <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;

    report << plain << "Using "<<reader.getType()<<" posfile \'"
	   <<config[InputPositions::keyword]<<"\' ("<<positions.size()<<")." << endr;

    trajectory.push_back( positions );

  }

  else{

    // Read dcd file and store each snapshot in trajectory.
    while((inputDCD >> positions ))
      trajectory.push_back( positions );

    // Check for no frames.
    if(trajectory.size() < 1)
      report << error << "Did not find any frames." << endr;

  }

  // Fix for old topology definition
  if (!config[GenericTopology::keyword].valid()){
    config[GenericTopology::keyword] = config[InputBoundaryConditions::keyword].getString()+config[InputCellManager::keyword].getString();
  }

  // create the appropriate topology (PBC or VBC, exclusion type, etc.)
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


  // Build the topology (uses the information in psf and par to create the molecule list, bonding lists,
  // Lennard-Jones parameter table and exclusion list)
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
    
  // get the # of atoms in the system
  unsigned int numAtoms = topo->atoms.size();

  // Fix velocities
  Vector3DBlock velocities(numAtoms);
  if(!config.valid(InputVelocities::keyword)){
    randomVelocity(config[InputTemperature::keyword],topo,&velocities,config[InputSeed::keyword]);
    report << plain << "Random temperature : "<<temperature(topo,&velocities)<<"K"<<endr;
  }

  // create a vector block of forces.  This is just a placeholder for the atomic forces that will be computed by
  // the nonbonded force objects.
  Vector3DBlock forces(numAtoms);

  // create the energy storage object (scalarstructure)
  ScalarStructure energies;


  // Create the forces and integrators  
  registerForceExemplars(topo);
  Integrator*  integrator = IntegratorFactory::make(errMsg,config[InputIntegrator::keyword]);
  if(integrator == NULL)
    report << error << errMsg <<endr;
  Real timeStep = integrator->top()->getTimestep()* (int)config[InputOutputfreq::keyword];
  report << plain << "Time step : "<< timeStep << endl
	 << PROTOMOL_HR << endr;
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
  // Topology
  report << topo->print(&positions) <<endr;
  report << plain  << PROTOMOL_HR << endr;

  // Collect all forces ...
  ForceGroup fg;    
  vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
  for(unsigned int i=0;i<inter.size();i++){
    const vector<MakeableDefinition>& forces(inter[i].forces);
    for(unsigned j=0;j<forces.size();j++){
      const MakeableDefinition& force(forces[j]);
      report << plain << Force::scope<<" "<<force.id<<endr;;
      string str(force.id);
      for(unsigned int k=0;k<force.parameters.size();k++){
	if(force.parameters[k].value.valid()){
	  str += " "+force.parameters[k].keyword+" "+force.parameters[k].value.getString();
	}

	// print
	if(!force.parameters[k].keyword.empty()){
	  report << plain <<"      "<<force.parameters[k].keyword;
	}
	report << plain <<" "<<force.parameters[k].value.getString();
	if(!force.parameters[k].text.empty())
	  report << plain <<"\t # "<<force.parameters[k].text;
	report << endr;

      }
      Force* f = ForceFactory::make(errMsg,str);
      if(f){
	fg.addForce(f);
      }
      else {
	report << errMsg <<endr;
      }
    }
  }
  delete integrator;
  report << PROTOMOL_HR << endr;


  // Clear all factories
  TopologyFactory::unregisterAllExemplars();
  IntegratorFactory::unregisterAllExemplars();
  ForceFactory::unregisterAllExemplars();

  Parallel::sync();
  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
  TimerStatistic::timer[TimerStatistic::RUN].start();

  // loop over all trajactories and process
  for(unsigned int i=0; i<trajectory.size(); ++i) {
    if(i > 0){
      for(unsigned int j=0; j<numAtoms;j++)
	velocities[j] = (trajectory[i][j]-trajectory[i-1][j])* Constant::TIMEFACTOR / timeStep;
    }
    energies.clear();
    forces.zero();
    Parallel::distribute(&energies,&forces);
    fg.evaluateSystemForces(topo, &(trajectory[i]), &forces, &energies);
    fg.evaluateExtendedForces(topo, &(trajectory[i]), &velocities, &forces, &energies);
    Parallel::reduce(&energies,&forces);
    report <<plain <<"Step : ";
    report.setf(std::ios::right);
    report.width(10);
    report <<i<<", Time : ";
    report.width(18);
    report.setf(std::ios::showpoint|std::ios::fixed);
    report.precision(3);
    report << i*timeStep<<" [fs], PE : ";
    report.precision(4);
    report.width(16);
    report <<energies.potentialEnergy()<<" [kcal/mol]";
    report << ", Tavg : ";
    report.precision(4);
    report.width(10);
    report <<temperature(topo,&velocities)<<" [K]";
    report << ", V : ";
    report.precision(2);
    report.width(16);
    report <<topo->getVolume(trajectory[i])<<" [AA^3]"<<endr;
    report.reset();
  }
  TimerStatistic::timer[TimerStatistic::RUN].stop();
  TimerStatistic::timer[TimerStatistic::WALL].stop();


  // Clean up
  delete topo;

  report << allnodesserial << plain <<"Timing" 
	 << (Parallel::isParallel()? string(" ("+toString(Parallel::getId())+")"): string(""))
	 <<" : "<<TimerStatistic()<<"."<<endr;
  Parallel::finalize();

  return 0;

}
