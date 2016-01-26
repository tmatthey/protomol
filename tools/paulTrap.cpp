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

#include "DCDTrajectoryReader.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::string;
//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {


  // Parse
  if(argc < 2)
    report << quit << "usage: " << argv[0] << " <dcd input file>"<< endr;

  // Open
  string filename(argv[1]);
  if(!DCDTrajectoryReader(filename).tryFormat())
    report << error << "'" << filename << "' not DCD." << endr;
  DCDTrajectoryReader in(filename);
  if(!in)
    report << error << "Could not open '" << filename << "'." << endr;


  // Read
  std::vector<Vector3DBlock> trajectory;
  Vector3DBlock xyz;
  int n = 0;
  while((in >> xyz)){
    if(argc < 2 || atoi(argv[2]) == n)
      trajectory.push_back(xyz);
    n++;
  }
  if(trajectory.size() < 1)
    report << error << "Empty DCD file"<< endr;
  if(trajectory[0].size() < 1)
    report << error << "Empty DCD frame file"<< endr;

  unsigned int atomCount = trajectory[0].size();


  SystemForce* f1 = new NonbondedSimpleFullSystemForce<OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForce,false,NoConstraint> >(OneAtomPair<VacuumBoundaryConditions,UniversalSwitchingFunction,CoulombForce,false,NoConstraint>(CoulombForce(),UniversalSwitchingFunction()),64);

  PaulTrapExtendedForce<VacuumBoundaryConditions>* f2 = new PaulTrapExtendedForce<VacuumBoundaryConditions>(2.6e-09,2.6e-09,0.0,0.0,0.0,Vector3D(0.0,0.0,0.0)); 

  GenericTopology* topo = new Topology<VacuumBoundaryConditions,CubicCellManager>(1.0,ExclusionType::NONE,VacuumBoundaryConditions(),CubicCellManager(1.0));

  AtomType atomType;
  atomType.name   = "CA";
  atomType.mass   = 40.8;
  atomType.charge = 1;


  topo->atomTypes.push_back(atomType);

  Atom atom;
  atom.type = 0;
  atom.scaledCharge = topo->atomTypes[atom.type].charge*Constant::SQRTCOULOMBCONSTANT;
  atom.scaledMass    = topo->atomTypes[atom.type].mass;

  topo->atoms.resize(atomCount,atom);

  buildExclusionTable(topo,topo->exclude);

  Vector3DBlock velocities(atomCount);
  Vector3DBlock forces(atomCount);
  ScalarStructure energies;

  //Process
  report.precision(15);
  for(unsigned int i=0;i<trajectory.size();++i){
    energies.clear();
    f1->evaluate(topo,&(trajectory[i]),&forces,&energies);
    f2->evaluate(topo,&(trajectory[i]),&velocities,&forces,&energies);
    Real pe = energies[ScalarStructure::OTHER]+energies[ScalarStructure::COULOMB];
    report << i << " \t " << pe <<" \t " <<  f2->cohesive(pe) <<" \t " <<  f2->epsilonEnergy(pe) <<" \t " << endr;
  }

  delete topo;
  delete f1;
  delete f2;
  
  return 0;
}
