#include "scpenergy.h"
#include "SystemForce.h"
#include "buildExclusionTable.h"
#include "CellListEnumerator_standard.h"
#include "CubicCellManager.h"
#include "C1SwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "CnSwitchingFunction.h"
#include "CmpCnCnSwitchingFunction.h"
#include "CoulombBornRadiiForce.h"
#include "CoulombSCPISMForce.h"
#include "CoulombSCPISMParameterTable.h"
#include "CutoffSwitchingFunction.h"
#include "ExclusionType.h"
#include "ShiftSwitchingFunction.h"
#include "NonbondedCutoffBornForce.h"
#include "NonbondedCutoffSystemForce.h"
#include "OneAtomPair.h"
#include "VacuumBoundaryConditions.h"
#include "Vector.h"
#include "Torsion.h"
#include "Report.h"
using namespace ProtoMol::Report;

SCPEnergy scpenergy(int numatoms, 
		    Real positions[], 
		    Real charges[], 
		    SCPBond bonds[],
		    int numbonds,
		    SCPAngle angles[],
		    int numangles,
		    SCPDihedral dihedrals[],
		    int numdihedrals,
		    string* atomTypes, 
		    int* residue_seq,
		    string exclusion,
		    Real cutoff,
		    Real D,
		    int s,
		    string switching,
		    Real switchon,
		    Real switchoff,
		    Real order) {
  SCPEnergy retval;
  retval.force = new Real[numatoms*3];
  
  // Initialize the SCPISM forces based on the switching function passed
  SystemForce* myForce; // CoulombSCPISM
  SystemForce* myForce2; // CoulombBornRadii

  if (switching == "C1") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),C1SwitchingFunction(1000)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,C1SwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),C1SwitchingFunction(cutoff)));
  }
  else if (switching == "C2") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),C2SwitchingFunction(switchon,1000)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,C2SwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),C2SwitchingFunction(switchon,cutoff)));
  }
  else if (switching == "Cn") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),CnSwitchingFunction(switchon,1000,order,switchoff)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,CnSwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),CnSwitchingFunction(switchon,cutoff,order,switchoff)));
  }
  else if (switching == "CmpCnCn") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),CmpCnCnSwitchingFunction(switchon,1000,order,switchoff)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,CmpCnCnSwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),CmpCnCnSwitchingFunction(switchon,cutoff,order,switchoff)));
  }
  else if (switching == "Cutoff") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),CutoffSwitchingFunction(1000)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,CutoffSwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),CutoffSwitchingFunction(cutoff)));
  }
  else if (switching == "Shift") {
    myForce = new NonbondedCutoffSystemForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombSCPISMForce> >(1000, OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombSCPISMForce>(CoulombSCPISMForce(D),ShiftSwitchingFunction(1000)));
    myForce2 = new NonbondedCutoffBornForce<CubicCellManager,OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombBornRadiiForce> >(cutoff, OneAtomPair<VacuumBoundaryConditions,ShiftSwitchingFunction,CoulombBornRadiiForce>(CoulombBornRadiiForce(s),ShiftSwitchingFunction(cutoff)));
  }

  // Populate positions
  Vector3DBlock posvec;
  for (unsigned int i = 0; i < numatoms; i++)
    posvec.push_back(Vector3D(positions[i*3],positions[i*3+1],positions[i*3+2]));
  // Forces and energies will be populated when we evaluate.
  // Just make them empty for now
  Vector3DBlock* forces = new Vector3DBlock(numatoms);
  ScalarStructure ss;

  // Set up topology.  Just need the atoms array.
  // JI: and the exclusion table!!!
  //Topology<VacuumBoundaryConditions, CubicCellManager>* topo = new Topology<VacuumBoundaryConditions,CubicCellManager>(),Vector<std::string>("NormalCubic");
  Topology<VacuumBoundaryConditions, CubicCellManager>* topo = new Topology<VacuumBoundaryConditions,CubicCellManager>(1.0, ExclusionType(exclusion), VacuumBoundaryConditions(), CubicCellManager(5));


  topo->atoms.resize(numatoms);

  // I realize this is redundant
  // But for this simple application it doesn't matter much.
  topo->atomTypes.resize(numatoms);
  topo->doSCPISM = true;

  // Populate SCPISM table
  CoulombSCPISMParameterTable mySCPISM;
  mySCPISM.populateTable();

  
  for (unsigned int i = 0; i < numatoms; i++) {
    topo->atomTypes[i].mySCPISM = new SCPISMAtomTypeParameters();
    topo->atoms[i].mySCPISM = new SCPISMAtomParameters();
    topo->atomTypes[i].mySCPISM->sqrt_alpha = mySCPISM.myData[atomTypes[i]].sqrt_alpha_i;
    topo->atomTypes[i].mySCPISM->g_i = mySCPISM.myData[atomTypes[i]].hbond_factor;
    topo->atomTypes[i].mySCPISM->isHbonded = mySCPISM.myData[atomTypes[i]].isHbonded;
    topo->atomTypes[i].mySCPISM->A_i = mySCPISM.myData[atomTypes[i]].A_i;
    topo->atomTypes[i].mySCPISM->B_i = mySCPISM.myData[atomTypes[i]].B_i;
    topo->atomTypes[i].mySCPISM->C_i = mySCPISM.myData[atomTypes[i]].C_i;
    Real R_vdw = mySCPISM.myData[atomTypes[i]].R_vdw + 1.40;
    topo->atoms[i].mySCPISM->dR_vdw2 = 1.0 / (4.0*M_PI*R_vdw*R_vdw);
    topo->atoms[i].mySCPISM->dR_vdw2 = 1.0 / (4.0*M_PI*R_vdw*R_vdw);
    topo->atoms[i].mySCPISM->r_cov = mySCPISM.myData[atomTypes[i]].r_cov;
    topo->atoms[i].type = i;
    // R_i,w is r_cov + 0.35 if charge negative, 0.85 for positive charge
    // Check this charge here....
    topo->atoms[i].mySCPISM->R_w = mySCPISM.myData[atomTypes[i]].r_cov + (charges[i] > 0 ? 0.85 : 0.35);
    topo->atoms[i].mySCPISM->R_p = topo->atoms[i].mySCPISM->R_w + mySCPISM.myData[atomTypes[i]].R_iw;
    topo->atoms[i].mySCPISM->R_iw = mySCPISM.myData[atomTypes[i]].R_iw;


    topo->atoms[i].mySCPISM->R_w = topo->atoms[i].mySCPISM->r_cov + (charges[i] > 0 ? 0.85 : 0.35);
    topo->atoms[i].mySCPISM->R_p = topo->atoms[i].mySCPISM->R_iw + topo->atoms[i].mySCPISM->R_w;
    //tempatom.residue_name = atom->residue_name;
    topo->atoms[i].residue_seq = residue_seq[i];
    // Now, the scaled charge.  This is straightforward.      
    topo->atoms[i].scaledCharge = (charges[i])*Constant::SQRTCOULOMBCONSTANT;
    //tempatom.scaledMass = atom->mass;
    topo->atoms[i].mySCPISM->sqrtalphaSCPISM = topo->atomTypes[i].mySCPISM->sqrt_alpha;
    topo->atoms[i].mySCPISM->sasaFrac = 0.0;
    topo->atoms[i].mySCPISM->polarFrac = 0.0;
    topo->atoms[i].mySCPISM->bornRadius = 0.0;
    topo->atoms[i].molecule = 0;
    // Now we need the size of the group for heavy atom ordering
    // We need to parse the name for any H's then any numbers following    
    // First, if the atom is an H then this is 0
    if (atomTypes[i] == "H"){
      topo->atoms[i].hvyAtom = 0;
    }
    else{
      // Otherwise, we need to parse..
      // Initialize to 1
      topo->atoms[i].hvyAtom = 1;
      for (unsigned int pos = 0; pos < atomTypes[i].size(); ++pos){
	if (atomTypes[i][pos] == 'H'){
	  string number = "";
	  while (isdigit(atomTypes[i][++pos]))   {
	    number += atomTypes[i][pos];
	  }
	  if (number == "") // never entered loop, default is 1
	    number = "1";
	  topo->atoms[i].hvyAtom += atoi(number.c_str());
	}
      }
    }    
    // C/C++ starts at 0, where PSF/PDB at 1
    topo->atoms[i].atomNum = i;
  }

  for (int i = 0; i < numbonds; i++)
    {
      topo->bonds.push_back(Bond());
      topo->bonds[i].atom1 = bonds[i].atom1-1;
      topo->bonds[i].atom2 = bonds[i].atom2-1;
    }

  for (int i = 0; i < numangles; i++)
    {
      topo->angles.push_back(Angle());
      topo->angles[i].atom1 = angles[i].atom1-1;
      topo->angles[i].atom2 = angles[i].atom2-1;
      topo->angles[i].atom3 = angles[i].atom3-1;
    }

  for (int i = 0; i < numdihedrals; i++)
    {
      topo->dihedrals.push_back(Torsion());
      topo->dihedrals[i].atom1 = dihedrals[i].atom1-1;
      topo->dihedrals[i].atom2 = dihedrals[i].atom2-1;
      topo->dihedrals[i].atom3 = dihedrals[i].atom3-1;
      topo->dihedrals[i].atom4 = dihedrals[i].atom4-1;
    }

  buildExclusionTable(topo,topo->exclude);

  //JI: debug
  report <<plain << topo->print()<< endl;
  //JI: end debug
  for (unsigned int i = 0; i < numatoms; i++)
    report << plain <<" "<< posvec[i] << endl;
  myForce2->evaluate(topo, &posvec, forces, &ss);
  myForce->evaluate(topo, &posvec, forces, &ss);

  retval.energy = ss[ScalarStructure::COULOMB];
  for (unsigned int i = 0; i < numatoms; i++) {
    retval.force[i*3] = (*forces)[i].x;
    retval.force[i*3+1] = (*forces)[i].y;
    retval.force[i*3+2] = (*forces)[i].z;
  }

  delete myForce;
  delete myForce2;
  return retval;
}

