// A sample driver file for Alanine Dipeptide

#include "scpenergy.h"
#include "Report.h"
using namespace ProtoMol::Report;


int main(int argc, char **argv) {
  const int numatoms = 22;
  const int numbonds = 21;
  const int numangles = 36;
  const int numdihedrals = 42;
  Real positions[numatoms*3] = {-2.4791, 10.3596, -4.2879,
				-2.2157, 9.59246, -3.5288,
				-1.8183, 10.2215, -5.1702,
				-3.5299, 10.1911, -4.6069,
				-2.1984, 11.7472, -3.8232,
				-1.6780, 12.5811, -4.5625,
				-2.4869, 12.0468, -2.5441,
				-2.8563, 11.3945, -1.8868,
				-2.1397, 13.3006, -1.9507,
				-1.7579, 13.9638, -2.7129,
				-3.3526, 13.9409, -1.2548,
				-3.0982, 14.9534, -0.8747,
				-4.2467, 14.0262, -1.9087,
				-3.6121, 13.2959, -0.3883,
				-1.0123, 13.0612, -1.0064,
				-0.9779, 12.0400, -0.3216,
				-0.0407, 13.9911, -0.9842,
				-0.0611, 14.8201, -1.5377,
				1.08845, 13.8258, -0.1226,
				1.35114, 14.7907, 0.36140,
				0.88264, 13.0787, 0.67341,
				1.96441, 13.4788 ,-0.7113};

  
  Real charges[numatoms] = {-0.270000,0.090000,0.090000,0.090000,
		      0.510000,-0.510000,-0.470000,0.310000,
		      0.070000,0.090000,
		      -0.270000,0.090000,0.090000,0.090000,
		      0.510000,-0.510000,-0.470000,0.310000,
		      -0.110000,0.090000,0.090000,0.090000};

  string atomTypes[numatoms] = {"CT3","HA","HA","HA",
			  "C","O","NH1","H",
			  "CT1","HB",
			  "CT3","HA","HA","HA",
			  "C","O","NH1","H",
			  "CT3","HA","HA","HA"};
  int rs[numatoms] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

  SCPBond bonds[numbonds] = {SCPBond(5, 1), SCPBond(5, 7), 
			     SCPBond(1, 2), SCPBond(1, 3),
			     SCPBond(1, 4), SCPBond(6, 5), 
			     SCPBond(11, 9), SCPBond(7, 8),
			     SCPBond(7, 9), SCPBond(15, 9), 
			     SCPBond(9, 10), SCPBond(11, 12),
			     SCPBond(11, 13), SCPBond(11, 14), 
			     SCPBond(16, 15), SCPBond(15, 17),
			     SCPBond(17, 18), SCPBond(17, 19), 
			     SCPBond(19, 20), SCPBond(19, 21),
			     SCPBond(19, 22)};  

  SCPAngle angles[numangles] = {SCPAngle(2,1,3),SCPAngle(2,1,4),
				SCPAngle(2,1,5),SCPAngle(3,1,4),
				SCPAngle(3,1,5),SCPAngle(4,1,5),
				SCPAngle(1,5,6),SCPAngle(1,5,7),
				SCPAngle(6,5,7),SCPAngle(5,7,8),
				SCPAngle(5,7,9),SCPAngle(8,7,9),
				SCPAngle(7,9,10),SCPAngle(7,9,11),
				SCPAngle(7,9,15),SCPAngle(10,9,11),
				SCPAngle(10,9,15),SCPAngle(11,9,15),
				SCPAngle(9,11,12),SCPAngle(9,11,13),
				SCPAngle(9,11,14),SCPAngle(12,11,13),
				SCPAngle(12,11,14),SCPAngle(13,11,14),
				SCPAngle(9,15,16),SCPAngle(9,15,17),
				SCPAngle(16,15,17),SCPAngle(15,17,18),
				SCPAngle(15,17,19),SCPAngle(18,17,19),
				SCPAngle(17,19,20),SCPAngle(17,19,21),
				SCPAngle(17,19,22),SCPAngle(20,19,21),
				SCPAngle(20,19,22),SCPAngle(21,19,22)};

  SCPDihedral dihedrals[numdihedrals] = {SCPDihedral(1,5,7,8),
					SCPDihedral(1,5,7,9),
					SCPDihedral(2,1,5,6),
					SCPDihedral(2,1,5,7),
					SCPDihedral(3,1,5,6),
					SCPDihedral(3,1,5,7),
					SCPDihedral(4,1,5,6),
					SCPDihedral(4,1,5,7),
					SCPDihedral(5,7,9,10),
					SCPDihedral(5,7,9,11),
					SCPDihedral(5,7,9,15),
					SCPDihedral(6,5,7,8),
					SCPDihedral(6,5,7,9),
					SCPDihedral(7,9,11,12),
					SCPDihedral(7,9,11,13),
					SCPDihedral(7,9,11,14),
					SCPDihedral(7,9,15,16),
					SCPDihedral(7,9,15,17),
					SCPDihedral(8,7,9,10),
					SCPDihedral(8,7,9,11),
					SCPDihedral(8,7,9,15),
					SCPDihedral(9,15,17,18),
					SCPDihedral(9,15,17,19),
					SCPDihedral(10,9,11,12),
					SCPDihedral(10,9,11,13),
					SCPDihedral(10,9,11,14),
					SCPDihedral(10,9,15,16),
					SCPDihedral(10,9,15,17),
					SCPDihedral(11,9,15,16),
					SCPDihedral(11,9,15,17),
					SCPDihedral(12,11,9,15),
					SCPDihedral(13,11,9,15),
					SCPDihedral(14,11,9,15),
					SCPDihedral(15,17,19,20),
					SCPDihedral(15,17,19,21),
					SCPDihedral(15,17,19,22),
					SCPDihedral(16,15,17,18),
					SCPDihedral(16,15,17,19),
					SCPDihedral(18,17,19,20),
					SCPDihedral(18,17,19,21),
					SCPDihedral(18,17,19,22),
					 SCPDihedral(19,17,15,9)};
  // choices are: "none", "1-2", "1-3", "1-4", "scaled1-4"
  string exclusion = "1-3";

  SCPEnergy myValues = scpenergy(numatoms, positions, charges, bonds, numbonds, angles, numangles, dihedrals, numdihedrals, atomTypes, rs,exclusion);
  report << plain << "ENERGY: " << myValues.energy << endr;
  report << plain << "FORCES: " << endr;
  for (unsigned int i = 0; i < numatoms; i++)
    report << plain <<atomTypes[i]<<" "<< myValues.force[i*3] << " " << myValues.force[i*3+1] << " " << myValues.force[i*3+2] << endr;
  delete myValues.force;
  return 0;
}
