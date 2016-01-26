#include "iSGPAR.h"

#include "mathutilities.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ iSGPAR
  const Real iSGPAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28(pow(2.0,-5.0/6.0));
  const Real iSGPAR::Nonbonded::SIGMA_CHARMM28_TO_CHARMM19(pow(2.0, 5.0/6.0));

  void iSGPAR::clear() {
    bonds.clear();
    angles.clear();
    dihedrals.clear();
    impropers.clear();
    nonbondeds.clear();
    nbfixs.clear();
    hbonds.clear();
  }


  //_________________________________________________________________globals
  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Bond & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2;
 
    for (unsigned int i=0; i<p.forceConstant.size(); i++) 
      OS <<","<<p.forceConstant[i]<<","<<p.distance[i];
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Angle & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3;

    for (unsigned int i=0; i<p.forceConstant.size(); i++)    
      OS<<","<<p.forceConstant[i]<<","<<p.angleval[i]<<","<<p.ub_flag<<","<<p.k_ub[i]<<","<<p.r_ub[i];
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Dihedral & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3<<","<<p.atom4;

    for (unsigned int i=0; i<p.forceConstant.size(); i++) {
      OS<<","<<p.multiplicity[i];
      for(unsigned int j=0;j<p.forceConstant[i].size();++j)
	OS << ","<<p.forceConstant[i][j];
      for(unsigned int j=0;j<p.periodicity[i].size();++j)
	OS << ","<<p.periodicity[i][j];
      for(unsigned int j=0;j<p.phaseShift[i].size();++j)
	OS << ","<<p.phaseShift[i][j];
    }
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Improper & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3<<","<<p.atom4;
    for (unsigned int i=0; i<p.forceConstant.size(); i++) 
      OS<<","<<p.forceConstant[i]<< ","<<p.periodicity[i]<<","<<p.phaseShift[i];
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Nonbonded & p) {
    OS <<p.number<<","<<p.atom;
    for (unsigned int i=0; i<p.epsilon.size(); i++) 
      OS<<","<<p.polarizability[i]<<","<<p.epsilon[i]<<","<<p.sigma[i]<<","<<p.negative[i]<< ","
	<<p.vdw[i]<<","<<p.polarizability2[i]<<","<<p.epsilon14[i]<<","<<p.sigma14[i]<<","<<p.negative2[i];
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Nbfix & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.a<<","<<p.b<<","<<p.a14<< ","<<p.b14;
    return OS;
  }

  MyStreamer& operator<< (MyStreamer& OS, const iSGPAR::Hbond & p) {
    OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.emin<<","<<p.rmin;
    return OS;
  }
}
