#include "OutputOSGXYZTrajectoryPos.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "XYZTrajectoryWriter.h"
#include "inputValueDefinitions.h"
#include "ScalarStructure.h"

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ OutputOSGXYZTrajectoryPos
  const string  OutputOSGXYZTrajectoryPos::keyword("oSGXYZPosFile");

  OutputOSGXYZTrajectoryPos::OutputOSGXYZTrajectoryPos():Output(),myXYZ(NULL),myMinimalImage(false),myComponent("all"){}

  OutputOSGXYZTrajectoryPos::OutputOSGXYZTrajectoryPos(const string& filename,
                                                       int freq,
                                                       bool minimal,
                                                       const string& compToWrite):Output(freq),
                                                                                  myXYZ(new XYZTrajectoryWriter(filename)),
                                                                                  myMinimalImage(minimal),
                                                                                  myComponent(compToWrite) {}

  OutputOSGXYZTrajectoryPos::~OutputOSGXYZTrajectoryPos(){
    if(myXYZ != NULL)
      delete myXYZ;
  }

  void OutputOSGXYZTrajectoryPos::doInitialize(){
  
    if (myComponent != "all") {
      // loop thru all the molecules in the topology and find the name of the component(molecule) whose
      // coordinates we are to output.  This is a check to make sure the user gave us a valid component
      // name for the parameter oSGXYZposfileComponent
      bool validCompName = false;
    
      for (unsigned int m = 0; m < myTopology->molecules.size(); m++) {
        if (myComponent == myTopology->molecules[m].name) {
          validCompName = true;
          break;}
      } // end loop over molecules
    
      if (validCompName)
        report << plain << "The coordinates of all " << myComponent << " molecules will be written to an XYZ trajectory file." << endr;
      else
        report << error << "Could not find any " << myComponent << " molecules to output!" << endr;

      numAtoms = 0;
      // loop thru each molecule to be outputted and count the # of atoms we need to write output
      for (unsigned int m = 0; m < myTopology->molecules.size(); m++) {
        if (myComponent == myTopology->molecules[m].name) {
          // we have found a molecule to be outputted, loop thru the list of atoms on this molecule
          for (unsigned int a=0; a<myTopology->molecules[m].atoms.size(); a++) numAtoms++;
        } // if statement
      } // end loop over molecules
      report << plain << "The XYZ coordinates of " << numAtoms << " atoms will be written to the trajectory file." << endr;
    } // end if myComponent statement
    else report << plain << "The XYZ coordinates of every atom will be written to the trajectory file." << endr;
  
    if(myXYZ == NULL || !myXYZ->open())
      report << error <<" Can not open \'"<<(myXYZ!=NULL?myXYZ->getFilename():"")<<"\' for "<<getId()<<"."<<endr;
  }

  // here is where we write the coordinate set
  void OutputOSGXYZTrajectoryPos::doRun(int mystep) {

    if (myEnergies->trajectory()) {

      /// here we copy over only the necessary coordinates from myPositions
      // if myComponent = "all" then we are to output the coordinates of all the atoms in the system
      if (myComponent == "all") {
        const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
        if(!myXYZ->write(*pos,myTopology->atoms,myTopology->atomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myXYZ->getFilename()<<"\'."<<endr;
      }

      // otherwise we output only the coordinates of the specified component (molecule type)
      else {
        const Vector3DBlock*       TempPos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
        Vector3DBlock*             ComponentPos;
        std::vector<Atom>          ComponentAtoms;
        std::vector<AtomType>      ComponentAtomTypes;
        int                        count = 0;

        ComponentPos = new Vector3DBlock();
        ComponentPos->resize(numAtoms); 

        // first loop over all molecules and find the ones we need to output
        for (unsigned int m=0; m<myTopology->molecules.size(); m++) {

          if (myTopology->molecules[m].name == myComponent) {
            // we have found a molecule to be outputted, loop thru the list of atoms on this molecule
            for (unsigned int a=0; a<myTopology->molecules[m].atoms.size(); a++) {

              // get the ID# of this atom
              unsigned int myID = myTopology->molecules[m].atoms[a];

              // get the coordinates, atom object and atomType object for this atom
              (*ComponentPos)[count] = (*TempPos)[myID];
              ComponentAtoms.push_back(myTopology->atoms[myID]);
              ComponentAtomTypes.push_back(myTopology->atomTypes[ myTopology->atoms[myID].type ]);
              count++;
            } // end for loop over atoms
          } // end if statement
        } // end loop over molecules

        // now output the data we just compiled
        if(!myXYZ->write(*ComponentPos, ComponentAtoms, ComponentAtomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myXYZ->getFilename()<<"\'."<<endr;

      } // end else statement
    } // end if (myEnergies->trajectory()) statement
  } // end function doRun(int)

  void OutputOSGXYZTrajectoryPos::doFinalize(int){
    myXYZ->close();
  }

  Output* OutputOSGXYZTrajectoryPos::doMake(string&, const vector<Value>& values) const{
    return (new OutputOSGXYZTrajectoryPos(values[0],values[1],values[2],values[3]));
  }

  void OutputOSGXYZTrajectoryPos::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myXYZ!=NULL?myXYZ->getFilename():"",ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),Text("whether the coordinates should be transformed to minimal image or not")));
    
    // the following string parameter will determine if we are to output the trajectory of a single
    // component or if we are to output the trajectory of all the atoms in the system (the default)
    parameter.push_back(Parameter(keyword+"Component",Value(myComponent,ConstraintValueType::NotEmpty()),"all"));
  }

  bool OutputOSGXYZTrajectoryPos::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    if(config->valid(InputMinimalImage::keyword) && !values[2].valid())
      values[2] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }

}
