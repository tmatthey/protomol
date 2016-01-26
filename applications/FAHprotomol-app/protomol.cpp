#include "ConfigurationReader.h"
#include "InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"
#include "EigenvectorReader.h"
#include "StateRestore.h"
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

#include "SemiGenericTopology.h"
#include "PeriodicBoundaryConditions.h"

#include "protomol_to_FAHgui.h"

using namespace ProtoMol;

#ifdef USE_CHECKPOINT
#include <cstdlib>
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ using ::rand; }
#endif

// the following is required to be sure the "EXPORT" works if it is used
#define CUSTOM_ARCHIVE_TYPES portable_binary_oarchive,portable_binary_iarchive
#include <boost/archive/basic_binary_oarchive.hpp>
#include <boost/archive/basic_binary_oprimitive.hpp>
#include <boost/archive/impl/basic_binary_oarchive.ipp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/detail/common_oarchive.hpp>
#include <boost/archive/detail/oserializer.hpp>
#include <boost/archive/binary_oarchive.hpp>
class portable_binary_oarchive :
    // don't derive from binary_oarchive !!!
    public boost::archive::binary_oarchive_impl<portable_binary_oarchive>
{
  typedef portable_binary_oarchive derived_t;
  friend class boost::archive::detail::common_oarchive<derived_t>;
  friend class boost::archive::basic_binary_oarchive<derived_t>;
  friend class boost::archive::basic_binary_oprimitive<derived_t, std::ostream>;
  friend class boost::archive::save_access;
  using boost::archive::binary_oarchive_impl<derived_t>::save;
 public:
  portable_binary_oarchive(std::ostream & os, unsigned int flags) :
    boost::archive::binary_oarchive_impl<portable_binary_oarchive>(os, flags) {}
};

#include <boost/archive/basic_binary_iarchive.hpp>
#include <boost/archive/basic_binary_iprimitive.hpp>
#include <boost/archive/impl/basic_binary_iarchive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>
#include <boost/archive/detail/common_iarchive.hpp>
#include <boost/archive/detail/iserializer.hpp>
#include <boost/archive/binary_iarchive.hpp>
class portable_binary_iarchive :
    // don't derive from binary_iarchive !!!
    public boost::archive::binary_iarchive_impl<portable_binary_iarchive>
{
public:
  ~portable_binary_iarchive() {}
private:
  typedef portable_binary_iarchive derived_t;
  friend class boost::archive::detail::common_iarchive<derived_t>;
  friend class boost::archive::basic_binary_iarchive<derived_t>;
  friend class boost::archive::basic_binary_iprimitive<derived_t, std::istream>;
  friend class boost::archive::load_access;
  using boost::archive::binary_iarchive_impl<derived_t>::load;
 public:
  portable_binary_iarchive(std::istream & is, unsigned int flags) :
    boost::archive::binary_iarchive_impl<portable_binary_iarchive>(is, flags) {}
};



void readArchive(portable_binary_iarchive& myInArchive,
		 StateRestore& rst)
{
    myInArchive >> rst.mySystemsize;
    // Integrator levels
    myInArchive >> rst.myLevels;
    // Then the timestep
    myInArchive >> rst.myTimestep;
    // Simulation box (All are 0 for VBC)
    myInArchive >> rst.myE1.x;
    myInArchive >> rst.myE1.y;
    myInArchive >> rst.myE1.z;
    myInArchive >> rst.myE2.x;
    myInArchive >> rst.myE2.y;
    myInArchive >> rst.myE2.z;
    myInArchive >> rst.myE3.x;
    myInArchive >> rst.myE3.y;
    myInArchive >> rst.myE3.z;
    myInArchive >> rst.myOrigin.x;
    myInArchive >> rst.myOrigin.y;
    myInArchive >> rst.myOrigin.z;
//     // Atom information
     rst.myPositions->resize(rst.mySystemsize);
     rst.myVelocities->resize(rst.mySystemsize);
     rst.myForces.resize(rst.myLevels+1);
     for (unsigned int i = 0; i <= rst.myLevels; i++) {
       rst.myForces[i] = new Vector3DBlock();
       rst.myForces[i]->resize(rst.mySystemsize);
     }
     for (unsigned int i = 0; i < rst.mySystemsize; i++) {
       // Position
       myInArchive >> (*(rst.myPositions))[i].x;
       myInArchive >> (*(rst.myPositions))[i].y;
       myInArchive >> (*(rst.myPositions))[i].z;
       // Velocity
       myInArchive >> (*(rst.myVelocities))[i].x;
       myInArchive >> (*(rst.myVelocities))[i].y;
       myInArchive >> (*(rst.myVelocities))[i].z;
       // Force
       for (unsigned int j = 0; j <= rst.myLevels; j++) {
	 myInArchive >> (*((rst.myForces)[j]))[i].x;
	 myInArchive >> (*((rst.myForces)[j]))[i].y;
	 myInArchive >> (*((rst.myForces)[j]))[i].z;
       }
     }
//     // Energies
     unsigned int first = (rst.myEnergies)->FIRST;
     unsigned int llast = (rst.myEnergies)->LAST;
     for (unsigned int i = first; i < llast; i++) {
       Real r;
       myInArchive >> r;
       (*(rst.myEnergies))[(ScalarStructure::Index)i] = r;
     }
}


void writeArchive(portable_binary_oarchive& myArchive,
		  GenericTopology* topo,
		  Vector3DBlock positions,
		  Vector3DBlock velocities,
		  Integrator* integrator,
		  ScalarStructure scalar)
{
      unsigned int syssize = positions.size();
      myArchive << syssize;
      unsigned int levels = integrator->level(); // THIS IS THE HIGHEST INTEGRATOR
      myArchive << levels;
      Real ts = topo->time;
      //Real ts = integrator->getTimestep();
      myArchive << ts;
      if (((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.getKeyword() == "Periodic") {
	myArchive << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e1().z;
	myArchive << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e2().z;
	myArchive << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.e3().z;
	myArchive << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().x << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().y << ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.origin().z;
      }
      else {
	Real zero = 0.0;
	myArchive << zero << zero << zero;
	myArchive << zero << zero << zero;
	myArchive << zero << zero << zero;
	myArchive << zero << zero << zero;
      }
      for (unsigned int i = 0; i < positions.size(); i++) {
	myArchive << positions[i].x << positions[i].y << positions[i].z;
	myArchive << velocities[i].x << velocities[i].y << velocities[i].z;
	Integrator* tmp = integrator;
	for (unsigned int j = 0; j <= levels; j++) {
	  myArchive << (*(tmp->getForces()))[i].x << (*(tmp->getForces()))[i].y << (*(tmp->getForces()))[i].z;
	  if (j != levels)
	    tmp = tmp->next();
	}
      }
      unsigned int first = scalar.FIRST;
      unsigned int last = scalar.LAST;
      for (unsigned int i = first; i < last; i++) {
	Real val = scalar[(ScalarStructure::Index)i];
	myArchive << val;
      }

}
		 
#endif

using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
using std::endl;
using std::string;
using std::vector;
//_____________________________________________________________________ protomol
#ifdef FAHCORE
extern "C" int fah_main( int argc, char **argv ){
#else
int main(int argc, char **argv) {
#endif

  // Redirect all output to a file
  //std::ofstream out("out.txt");
  //report.setStream(&(out));

  // Redirect all output to std::cout
  //report.setStream(&(std::cout));

  Parallel::init(argc,argv);

  TimerStatistic::timer[TimerStatistic::WALL].start();

  report << plain <<endl
	 <<"######                                  #     #"<<endl
	 <<"#     #  #####    ####   #####   ####   ##   ##   ####   #"<<endl
	 <<"#     #  #    #  #    #    #    #    #  # # # #  #    #  #"<<endl
	 <<"######   #    #  #    #    #    #    #  #  #  #  #    #  #"<<endl
	 <<"#        #####   #    #    #    #    #  #     #  #    #  #"<<endl
	 <<"#        #   #   #    #    #    #    #  #     #  #    #  #"<<endl
	 <<"#        #    #   ####     #     ####   #     #   ####   ######"<<endl
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
  InputShake::registerConfiguration(&config,false);
  InputShakeEpsilon::registerConfiguration(&config,1e-5);
  InputShakeMaxIter::registerConfiguration(&config,30);
  InputRattle::registerConfiguration(&config,false);
  InputRattleEpsilon::registerConfiguration(&config,1e-5);
  InputRattleMaxIter::registerConfiguration(&config,30);
  InputReducedImage::registerConfiguration(&config,false);
  InputMinimalImage::registerConfiguration(&config,true);
  InputShadow::registerConfiguration( &config, false );
  InputShadowOrder::registerConfiguration( &config, 4 );
  InputShadowFreq::registerConfiguration( &config, 1 );
  InputEigenVectors::registerConfiguration(&config);
  InputEigenValues::registerConfiguration(&config);
  InputCheckpointFreq::registerConfiguration(&config);
  InputCheckpointFile::registerConfiguration(&config);
  InputRestoreState::registerConfiguration(&config);

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

  InputPosVel reader;
  Vector3DBlock positions;
  Vector3DBlock velocities;
#ifdef USE_CHECKPOINT
  // ARCHIVE READING
  bool doRestore = false;
  StateRestore rst;  
  std::string fn;
  if (config.valid(InputRestoreState::keyword)) {
    config[InputRestoreState::keyword].get(fn);
    std::ifstream inFile(fn.c_str(), ios::in);
    portable_binary_iarchive myInArchive(inFile, ios::in);
    doRestore = true;
    readArchive(myInArchive, rst);
    positions.resize(rst.mySystemsize);
    velocities.resize(rst.mySystemsize);
    positions.intoAssign(*(rst.myPositions));
    velocities.intoAssign(*(rst.myVelocities));
  }
#endif
  // 
  // Positions
  //
  if(!reader.open(config[InputPositions::keyword]))
    report << error << "Can't open position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;
  vector<PDB::Atom> pdbAtoms;
  if(reader.tryFormat(InputPosVelType::PDB)){
    PDB pdb;
    if(!(reader >> pdb))
      report << error << "Could not parse PDB position file \'"<<config[InputPositions::keyword]<<"\'."<<endr;
#ifdef USE_CHECKPOINT
    if (!doRestore)
      swap(positions,pdb.coords);
#endif
    swap(pdbAtoms,pdb.atoms);
  }
  else if(!(reader >> positions))
    report << error << "Could not parse position file \'"<<config[InputPositions::keyword]
	   <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;
  report << plain << "Using "<<reader.getType()<<" posfile \'"<<config[InputPositions::keyword]<<"\' ("<<positions.size()<<")." << endr;

  // 
  // Velocities
  //
#ifdef USE_CHECKPOINT
  if (!doRestore) {
#endif
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
#ifdef USE_CHECKPOINT
  }
#endif
  //
  // Eigenvectors/values
  //
  EigenvectorReader evReader;
  EigenvectorInfo ei;
  if(config.valid(InputEigenVectors::keyword)){
    if(!evReader.open(config[InputEigenVectors::keyword]))
      report << error << "Can't open eigenvector file \'"<<config[InputEigenVectors::keyword]<<"\'."<<endr;
    if(!(evReader >> ei)) {
      //report << error << "Could not parse eigenvector file \'"<<config[InputEigenVectors::keyword]<<endr;
      if (ei.myEigenvectorLength != (double)positions.size())
	report << error << "Eigenvector length is wrong, should be " << positions.size() << " got "<<ei.myEigenvectorLength<<"."<<endr;
      if (ei.myNumEigenvectors < 1 || ei.myNumEigenvectors > (double)positions.size())
	report << error << "Wrong number of eigenvectors (" << ei.myNumEigenvectors << ")." << endr;
    }
    report << plain << "Using "<<reader.getType()<<" eigfile \'"<<config[InputEigenVectors::keyword]<<"\' ("<<ei.myEigenvectorLength<<")." << endr;
  }

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
#ifdef USE_CHECKPOINT
    if (!doRestore)
#endif
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

  //Normal mode?
  string::size_type loc = ((string)config[InputIntegrator::keyword]).find( "NormMode", 0 );
  //string::size_type loc = integrator->keyword.find( "NormMode", 0 );
  if( loc != string::npos ){
	  if(!config.valid(InputEigenVectors::keyword)){ 
		  //now allow no file for re-diaganalization code, trap error in NormModeInt if still no Eigenvectors
          //report << error << "No Eigenvector file for NormMode integrator"<<endr;
		  integrator->mhQ = NULL;
		  integrator->maxEigval = 0;
		  integrator->numEigvects = 0;
	  }else{
		  integrator->mhQ = ei.myEigenvectors;
		  cout << "NULL: " << (integrator->mhQ == NULL) << endl;
		  cout << "VALUE: " << integrator->mhQ[0] << endl;
		  integrator->maxEigval = ei.myMaxEigenvalue;
		  integrator->numEigvects = ei.myNumEigenvectors;
	  }
      report << plain  << "Using Normal Mode integrator. " << endr;
  }else{
      report << plain  << "Using MD integrator. " << endr;
  }

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

  //  ----------------------------------------------------------------------  //
  //  Create a shadow Hamiltonian object.                                     //
  //  ----------------------------------------------------------------------  //
  if( ( bool )config[ InputShadow::keyword ] &&
      ( int )config[ InputShadowOrder::keyword ] > 0 &&
      ( int )config[ InputShadowFreq::keyword ] > 0 ) {

      integrator->adoptPostStepModifier(
              integrator->createShadowModifier(
                      (int)config[ InputShadowOrder::keyword ],
                      (int)config[ InputShadowFreq::keyword ] ) );

      report << plain << "Calculating shadow energy with order 2k = "
             << (int)config[ InputShadowOrder::keyword ] << ", and frequency = "
             << (int)config[ InputShadowFreq::keyword ] << "."
             << endr;

  }

  if( (bool)config[InputShake::keyword] &&
      (Real)config[InputShakeEpsilon::keyword] > 0.0 &&
      (int)config[InputShakeMaxIter::keyword] > 0) {

      // initialize the SHAKE modifier 
      integrator->bottom()->adoptPostDriftOrNextModifier(
                            integrator->bottom()->createShakeModifier(
                                (Real)config[InputShakeEpsilon::keyword],
                                (int)config[InputShakeMaxIter::keyword] ) );

      report << plain << "Shake with epsilon "
             << config[InputShakeEpsilon::keyword] << ", max "
             << (int)config[InputShakeMaxIter::keyword] << " iteration(s)."
             << endr;

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
#ifdef USE_CHECKPOINT
  if (doRestore) {
    scalar.intoAssign(*(rst.myEnergies));
    outputs->setRestore();
  }
#endif
  scalar.molecularVirial(config[InputMolVirialCalc::keyword]);
  scalar.virial(config[InputVirialCalc::keyword]);
  report << plain << "Virial tensor : "<< scalar.virial()<<endr;
  report << plain << "Molecular virial tensor : "<< scalar.molecularVirial()<<endr;
  
#ifdef USE_CHECKPOINT
  int myLevels = integrator->level();
  Integrator* tmp = integrator->bottom();
  if (doRestore) {
    topo->time = rst.myTimestep;
    config[InputFirststep::keyword] = rst.myTimestep / integrator->getTimestep();
  }
  else
#endif
    topo->time = (Real)config[InputFirststep::keyword]*integrator->getTimestep();
  integrator->initialize(topo,&positions,&velocities,&scalar);
#ifdef USE_CHECKPOINT
  // RESET THE STATE OF THE INTEGRATOR AFTER INITIALIZE IS CALLED.
  if (doRestore) {
    for (unsigned int i = 0; i <= myLevels; i++)
      {
 	// MAKE SURE THIS IS ONLY DONE ONE TIME!
	if (i == 0) {
	  swap(positions, *(rst.myPositions));
	  swap(velocities, *(rst.myVelocities));
	}
	swap(*(tmp->getForces()), *(rst.myForces[myLevels-i]));
	scalar.intoAssign(*(rst.myEnergies));
	buildMolecularCenterOfMass(&positions,topo);
	buildMolecularMomentum(&velocities,topo);
 	tmp->deleteInternalModifiers();
 	tmp->initializeModifiers();
	if (i != myLevels)
	  tmp = tmp->previous();
      }
  }
#endif
  outputs->initialize(&config,integrator,topo,&positions,&velocities,&scalar);
  outputs->addToCache(pdbAtoms);
  outputs->addToCache(psf);
  outputs->addToCache(par);

  // Print
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
#ifdef USE_CHECKPOINT
  if (doRestore) {
    string bc = config[InputBoundaryConditions::keyword].getString();
    if (bc == "Periodic") {
      ((SemiGenericTopology<PeriodicBoundaryConditions>*)topo)->boundaryConditions.set(rst.myE1, rst.myE2, rst.myE3, rst.myOrigin);
    }
  }
#endif
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
  if (doRestore)
    step = rst.myTimestep / integrator->getTimestep();

  int last = step+(int)config[InputNumsteps::keyword];

  Parallel::sync();
  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
  TimerStatistic::timer[TimerStatistic::RUN].start();

#ifdef FAHCORE
  char *tempName="protomolFAH_BETA";
  fcSetSystemName( tempName );
  gsInitData( (topo->atoms).size, (topo->bonds).size, tempName );
#endif
  

#ifdef USE_CHECKPOINT
  int ckptFreq, nextCkpt;
  if (config.valid(InputCheckpointFreq::keyword)) {
    ckptFreq = config[InputCheckpointFreq::keyword];
    if (config.valid(InputOutputfreq::keyword))
      {
	int outputFreq = config[InputOutputfreq::keyword];
	if (ckptFreq % outputFreq != 0)
	  report << error << "Output frequency of " << outputFreq << " is not an even factor of checkpoint frequency of " << ckptFreq << ".  This will uneven time jumps in output files." << endr;
      }
    nextCkpt = ckptFreq;
  }
  else {
    ckptFreq = -1;
    nextCkpt = last+1;
  }
#endif

  while(step < last) {
    outputs->run(step);
	#ifdef USE_CHECKPOINT
    Real myTime = topo->time;
    int inc = min(outputs->getNext(), nextCkpt)-step;
	#else
    int inc = outputs->getNext()-step;
	#endif
    inc = std::min(last,step+inc)-step;
    step += inc;
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();
    integrator->run(inc);
    TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();
	//PRB FAH Code to checkpoint and update GUI
	#ifdef FAHCORE
    //Handle Checkpointing

	//Update the GUI if needed
	//Currently energy and temperature are not computed correctly.
	gsPollRequest( gui_callback, topo, &positions, step, last, 0.0f, 0.0f );
    fcReportProgress( last, step );
	#endif
    
	#ifdef USE_CHECKPOINT
    if (doRestore && (inc > 0 && topo->time == myTime))
       topo->time += inc*integrator->getTimestep();
    if (step % ckptFreq == 0 && ckptFreq != -1) {
      // DO CHECKPOINTING
      string filename;
      if (config[InputCheckpointFile::keyword].valid()) {
	filename = config[InputCheckpointFile::keyword].getString();
      }
      else {
	time_t myTime;
	time(&myTime);
	filename = ctime(&myTime);
	for (unsigned int i = 0; i < filename.size(); i++)
	  if (filename[i] == ' ') filename[i] = '_';
	filename += "STEP";
	filename += toString(step);
	filename[filename.size()-1] = '.'; // Eliminate ?
	filename += "ar";
      }
      report << plain << "WRITING ARCHIVE..." << endr;
      std::ofstream myFile(filename.c_str(), ios::out);
      portable_binary_oarchive myArchive(myFile, ios::out);      
      writeArchive(myArchive, topo, positions, velocities, integrator, scalar);
      report << plain << "DONE." << endr;
      nextCkpt += ckptFreq;
    }
	#endif
	//PRB end this FAH code seg
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
