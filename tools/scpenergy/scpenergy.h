#ifndef SCPENERGY_H
#define SCPENERGY_H

#include "Real.h"
#include <string>
using namespace std;
using namespace ProtoMol;

// SCPEnergy is the structure returned;
// Contains a floating point value for the SCPISM energy
// And a floating point array of size 3*N for the forces
struct SCPEnergy {
  Real energy;
  Real* force;
};

struct SCPBond {
  SCPBond(int a1, int a2) : atom1(a1), atom2(a2) {}
  int atom1;
  int atom2;
};

struct SCPAngle {
  SCPAngle(int a1, int a2, int a3) : atom1(a1), atom2(a2), atom3(a3) {}
  int atom1;
  int atom2;
  int atom3;
};

struct SCPDihedral {
  SCPDihedral(int a1, int a2, int a3, int a4) : atom1(a1), atom2(a2), atom3(a3), atom4(a4) {}
  int atom1;
  int atom2;
  int atom3;
  int atom4;
};


// numatoms = N
// positions = atomic positions (size 3N)
// charges = atomic charges (size N)
// atomTypes = atom type names (size N)
// residue_seq = residue # 
// D = Dielectric constant
// s = Born switch
// switching = Switching function
// switchon = Distance to turn on switching
// switchoff = Distance to turn off switching
// order = For Cn switching functions, the order
// Possible switching: Universal, C1, C2, Cn, CmpCnCn, Shift
// choices for exclusion are: "none", "1-2", "1-3", "1-4", "scaled1-4"
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
		    Real cutoff=5.0,
		    Real D=80,
		    int s=1,
		    string switching="Cutoff",
		    Real switchon=0.0,
		    Real switchoff=0.0,
		    Real order=2);

#endif
