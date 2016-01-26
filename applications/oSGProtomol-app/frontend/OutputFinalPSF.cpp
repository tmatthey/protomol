#include "OutputFinalPSF.h"
#include "Configuration.h"
#include "OutputCache.h"        //change
#include "stringutilities.h"
#include "GenericTopology.h"    
#include "PSFWriter.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputFinalPSF
  const string  OutputFinalPSF::keyword("finPSFFile");

  OutputFinalPSF::OutputFinalPSF():Output(1),myFilename(""){}

  OutputFinalPSF::OutputFinalPSF(const string& filename):Output(1),myFilename(filename){}

  void OutputFinalPSF::doFinalize(int step){
    PSFWriter writer;
    if(!writer.open(myFilename))
      report << error << "Can't open "<<getId()<<" \'"<<myFilename<<"\'."<<endr;
    writer.setComment("Time : "+toString(myCache->time())+", step : "+toString(step)+".");

    // create a temporary PSF object
    PSF tempPSF = myCache->psf();
    
    // set the array sizes for tempPSF
    tempPSF.atoms.resize(myTopology->atoms.size());
    tempPSF.bonds.resize(myTopology->bonds.size());
    tempPSF.angles.resize(myTopology->angles.size());
    tempPSF.dihedrals.resize(myTopology->dihedrals.size());
    tempPSF.impropers.resize(myTopology->impropers.size());
    
    // loop over all atoms and copy the topology info into tempPSF
    for (unsigned int i=0; i<tempPSF.atoms.size(); i++) {
      tempPSF.atoms[i].identity         = myTopology->molecules[ myTopology->atoms[i].molecule ].type;
      tempPSF.atoms[i].mass             = myTopology->atoms[i].scaledMass;
      tempPSF.atoms[i].charge           = myTopology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
      tempPSF.atoms[i].number           = i + 1;
      tempPSF.atoms[i].seg_id           = myTopology->molecules[ myTopology->atoms[i].molecule ].name;
      tempPSF.atoms[i].atom_name        = myTopology->atoms[i].name;
      tempPSF.atoms[i].residue_sequence = myTopology->atoms[i].molecule + 1;
      tempPSF.atoms[i].residue_name     = myTopology->atomTypes[ myTopology->atoms[i].type ].name;
      tempPSF.atoms[i].atom_type        = myTopology->atomTypes[ myTopology->atoms[i].type ].name;
    }      

    ///for molecular systems I must add loops over bonds, etc.
      
    const PSF MyiSGPSF = tempPSF;
 
    if(!writer.write(MyiSGPSF))
      report << error << "Could not write "<<getId()<<" \'"<<myFilename<<"\'."<<endr;  
  }

  Output* OutputFinalPSF::doMake(string&, const vector<Value>& values) const{
    return (new OutputFinalPSF(values[0]));
  }

  void OutputFinalPSF::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
  }

}
