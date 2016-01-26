#include "XSCWriter.h"

#include <iomanip>

#include "Report.h"
#include "stringutilities.h"
#include "systemutilities.h"

using std::string;
using std::endl;
using std::setprecision;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________XSCWriter

  //~~~ empty constructor ~~~
  XSCWriter::XSCWriter():Writer(){}

  //~~~ constructor ~~~
  XSCWriter::XSCWriter(const std::string& filename):Writer(filename){}

  //~~~ write function that writes the XSC values to a file ~~~
  bool XSCWriter::write(const XSC& xsc){
    if (!open())
      return false;

    // integrator type
    myFile << "# Type of integrator used" << endl;
    myFile << "integrator   " << xsc.simType << endl;

    // Comment
    myFile << "# ProtoMol (built on "<< __DATE__ << " at " << __TIME__<< ") generated this XSC file by "<<getUserName()<<". " <<myComment<< endl;

    // Write xsc info
    myFile << setprecision(15); // This should be some FLT_DIG or DBL_DIG ...
 
    myFile << "# Molecule being transformed, it's old and new identities:" << endl;
    myFile << "molecule     " << xsc.myMolecule << endl;
    myFile << "oldtype      " << xsc.old_type << endl;
    myFile << "newtype      " << xsc.new_type << endl;

    myFile << "# Final Lambda and Lambda velocity:" << endl;
    myFile << "Lambda       " << xsc.Lambda << endl;
    myFile << "LambdaVel    " << xsc.Lambda_vel << endl;

    myFile << "# Final atom-thermostat value and velocity:" << endl;
    myFile << "Eta          " << xsc.Eta << endl;
    myFile << "EtaVel       " << xsc.Eta_vel << endl;

    myFile << "# Final Lambda-thermostat value and velocity:" << endl;
    myFile << "EtaLambda    " << xsc.EtaLambda << endl;
    myFile << "EtaLambdaVel " << xsc.EtaLambda_vel << endl;

    if (xsc.simType == "NfPTVerlet") {
      myFile << "# Final volume-thermostat value and velocity:" << endl;
      myFile << "EtaVol       " << xsc.EtaVol << endl;
      myFile << "EtaVolVel    " << xsc.EtaVol_vel << endl;

      myFile << "# Final volume and barostat velocity:" << endl;
      myFile << "Volume       " << xsc.Vol << endl;
      myFile << "EpsilonVel   " << xsc.Epsilon_vel << endl;
    }
    else {
      myFile << "# Final volume:" << endl;
      myFile << "Volume       " << xsc.Vol << endl;
    }
    close();
    return !myFile.fail();
  }

  XSCWriter& operator<<(XSCWriter& xscWriter, const XSC& xsc){
    xscWriter.write(xsc);
    return xscWriter;
  }
}
