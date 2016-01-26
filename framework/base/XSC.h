/*  -*- c++ -*-  */
#ifndef XSC_H
#define XSC_H

#include "Real.h"
#include <string>

namespace ProtoMol {
  //_________________________________________________________________XSC
  struct XSC{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XSC
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // XSC container
    std::string simType;

    // for chemostat
    Real Lambda, Lambda_vel;
    unsigned int myMolecule, old_type, new_type;
    
    // for thermostat
    Real Eta, EtaVol, Eta_vel, EtaVol_vel, EtaLambda, EtaLambda_vel;

    // for barostat
    Real Vol, Epsilon_vel;
  };

  //____________________________________________________________________________INLINES

}
#endif /* XSC_H */

