#include "XSC.h"

namespace ProtoMol {
  //_____________________________________________________________________ XSC
  void XSC::clear(){
    Lambda = 0.;
    Lambda_vel = 0.;
    myMolecule = 0;
    old_type = 0;
    new_type = 0;
    Eta = EtaVol = Eta_vel = EtaVol_vel = EtaLambda = EtaLambda_vel = 0.;
    Vol = Epsilon_vel = 0.;
  }

}

