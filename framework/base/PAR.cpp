#include "PAR.h"

#include "mathutilities.h"

//using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ PAR
  const Real PAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28(pow(2.0,-5.0/6.0));
  const Real PAR::Nonbonded::SIGMA_CHARMM28_TO_CHARMM19(pow(2.0, 5.0/6.0));

  void PAR::clear(){
    bonds.clear();
    angles.clear();
    dihedrals.clear();
    impropers.clear();
    nonbondeds.clear();
    nbfixs.clear();
    hbonds.clear();
  }


  //_________________________________________________________________globals
//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Bond & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.forceConstant<<","<<p.distance;
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Angle & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3<<","<<p.forceConstant<<","<<p.angleval<<","<<p.ub_flag<<","<<p.k_ub<<","<<p.r_ub;
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Dihedral & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3<<","<<p.atom4<<","<<p.multiplicity;
//     for(unsigned int i=0;i<p.forceConstant.size();++i)
//       OS << ","<<p.forceConstant[i];
//     for(unsigned int i=0;i<p.periodicity.size();++i)
//       OS << ","<<p.forceConstant[i];
//     for(unsigned int i=0;i<p.phaseShift.size();++i)
//       OS << ","<<p.forceConstant[i];
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Improper & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.atom3<<","<<p.atom4<<","<<p.forceConstant<< ","<<p.periodicity<<","<<p.phaseShift;
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Nonbonded & p) {
//     OS <<p.number<<","<<p.atom<<","<<p.polarizability<<","<<p.epsilon<<","<<p.sigma<<","<<p.negative<< ","
//        <<p.vdw<<","<<p.polarizability2<<","<<p.epsilon14<<","<<p.sigma14<<","<<p.negative2;
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Nbfix & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.a<<","<<p.b<<","<<p.a14<< ","<<p.b14;
//     return OS;
//   }

//   MyStreamer& operator<< (MyStreamer& OS, const PAR::Hbond & p) {
//     OS <<p.number<<","<<p.atom1<<","<<p.atom2<<","<<p.emin<<","<<p.rmin;
//     return OS;
//   }
}
