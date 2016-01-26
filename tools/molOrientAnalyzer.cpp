#include "Topology.h"
#include "oneAtomContraints.h"
#include "OneAtomPair.h"
#include "CoulombForce.h"
#include "NonbondedSimpleFullSystemForce.h"
#include "UniversalSwitchingFunction.h"
#include "CubicCellManager.h"
#include "VacuumBoundaryConditions.h"
#include "PaulTrapExtendedForce.h"
#include "buildTopology.h"

#include "InputPosVel.h"
#include "PARReader.h"
#include "PSFReader.h"
#include "TimerStatistic.h"
#include "XYZWriter.h"

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


  // 
  // Positions
  //
  InputPosVel reader;
  if(!reader.open(argv[1]))
    report << error << "Can't open position file \'"<<argv[1]<<"\'."<<endr;
  Vector3DBlock positions;
  vector<PDB::Atom> pdbAtoms;
  if(reader.tryFormat(InputPosVelType::PDB)){
    PDB pdb;
    if(!(reader >> pdb))
      report << error << "Could not parse PDB position file \'"<<argv[1]<<"\'."<<endr;
    swap(positions,pdb.coords);
    swap(pdbAtoms,pdb.atoms);
  }
  else if(!(reader >> positions))
    report << error << "Could not parse position file \'"<<argv[1]
	   <<"\'. Supported formats are : "<<InputPosVelType::getPossibleValues(", ")<<"."<<endr;
  report << plain << "Using "<<reader.getType()<<" posfile \'"<<argv[1]<<"\' ("<<positions.size()<<")." << endr;

  // 
  // Velocities
  //
  Vector3DBlock velocities;
  velocities.zero(positions.size());

  //
  // PSF
  //
  PSFReader psfReader;
  if(!psfReader.open(argv[2]))
    report << error << "Can't open PSF file \'"<<argv[2]<<"\'."<<endr;
  PSF psf;
  if(!(psfReader >> psf))
    report << error << "Could not parse PSF file \'"<<argv[2]<<"\'."<<endr;
  report << plain << "Using PSF file \'"<<argv[2]<<"\' ("<<psf.atoms.size()<<")." << endr;

  //
  // PAR
  //
  PARReader parReader;
  if(!parReader.open(argv[3]))
      report << error << "Can't open PAR file \'"<<argv[3]<<"\'."<<endr;
  PAR par;
  if(!(parReader >> par))
    report << error << "Could not parse PAR file \'"<<argv[3]<<"\'."<<endr;
  report << plain << "Using PAR file \'"<<argv[3]<<"\', "<<(parReader.getCharmmTypeDetected() != PAR::CHARMM28?"old":"new")<< " charmm force field.";
  if(parReader.getCharmmTypeDetected() != PAR::CHARMM28)
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
  GenericTopology* topo = new Topology<VacuumBoundaryConditions,CubicCellManager>(1.0,ExclusionType::ONE3,VacuumBoundaryConditions(),CubicCellManager(5.0));

  // Build topology
  buildTopology(topo,psf,par,parReader.getCharmmTypeDetected() != PAR::CHARMM28);

  if(!topo->checkMoleculePairDistances(positions)){
    Vector3DBlock tmp(positions);
    topo->minimalImage(tmp);
    if(topo->checkMoleculePairDistances(tmp)){
      positions = tmp;
      report << plain <<"Fixed minimal molecule distances."<<endr;
    }
    else {
      report << plain <<"Could not fixed minimal molecule distances."<<endr;      
    }    
  }

  report << "p = [ ";

  
  XYZ pos;
  for(int i=0;i<topo->molecules.size();++i){
    Vector3DBlock mol(topo->molecules[i].size());
    for(int j=0;j<topo->molecules[i].size();++j){
      mol[j]=positions[topo->molecules[i][j]];
    }
    Vector3D n;
    Real e,dc;
    if(!mol.fitplane(n,dc,e))
      report << "[0 0 0] ";
    else
      report << "["<<n.x <<" "<<n.y <<" "<<n.z<<"] ";          
    pos.coords.push_back(n);
    pos.names.push_back("M"+toString(topo->molecules[i].size()));
  }
  report << "]"<<endr;
  XYZWriter xyz("n.xyz");
  xyz << pos;

  // Clean up
  delete topo;
  return 0;
}
