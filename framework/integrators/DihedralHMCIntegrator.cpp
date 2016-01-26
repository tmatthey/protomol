#include "DihedralHMCIntegrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "Topology.h"

#include <fstream>

using namespace ProtoMol::Report;
using std::vector;
using std::string;
using std::endl;


namespace ProtoMol {
  //____________________________________________________________ DihedralHMCIntegrator

  const string DihedralHMCIntegrator::keyword("DihedralHMC");
  DihedralHMCIntegrator::DihedralHMCIntegrator()    
    : MTSIntegrator(),
      myOldPositions(NULL),
      myOldVelocities(NULL),
      myOldEnergies(NULL),
      myInitialTemperature(-1),
      myDihedralIndex(-1){}

  DihedralHMCIntegrator::DihedralHMCIntegrator(int cycles,
					       Real initialTemperature,
					       bool /*randomCycLen*/,
					       bool dihset,
					       std::string dsetfile,
					       bool angset,
					       std::string asetfile,                                               
					       ForceGroup *overloadedForces,
					       StandardIntegrator *nextIntegrator)
    
    : MTSIntegrator(cycles,overloadedForces,nextIntegrator),
      myOldPositions(new Vector3DBlock),
      myOldVelocities(new Vector3DBlock),
      myOldEnergies(new ScalarStructure),
      myInitialTemperature(initialTemperature),
      myDihedralsSet(dihset),
      myDhmcDiSetFile(dsetfile),
      myAnglesSet(angset),
      myDhmcAnSetFile(asetfile),
      myDihedralIndex(-1),      
      myDihedrals(new vector< int >()),
      myAngles(new vector< Real >()) {}
  
  DihedralHMCIntegrator::~DihedralHMCIntegrator(){
    if(myOldPositions != NULL)
      delete myOldPositions;
    if(myOldVelocities != NULL)
      delete myOldVelocities;
    if(myOldEnergies != NULL)
      delete myOldEnergies;
  }

  void DihedralHMCIntegrator::initialize(GenericTopology *topo,
					 Vector3DBlock *positions,
					 Vector3DBlock *velocities,
					 ScalarStructure *energies){

    MTSIntegrator::initialize(topo,positions,velocities,energies);

    if(myDihedralsSet){
      std::ifstream dihedralsSetinput(string(myDhmcDiSetFile).c_str(), std::ios::in);
      if(!dihedralsSetinput)
        report << error <<" Could not open \'"<<myDhmcDiSetFile<<"\' for "<<getId()<<"."<<endr;
      int tempdihedral = 1;
      while (dihedralsSetinput >> tempdihedral)
      {
        if(tempdihedral < 0 || tempdihedral >= (int)myTopo->dihedrals.size())
        {
          report << error << "[DihedralHMCIntegrator::initialize] Dihedral index "<< tempdihedral
                 <<" out of range [0,"<<myTopo->dihedrals.size()<<"[."<<endr;
        }
        myDihedrals->push_back(tempdihedral);
      }
    }

    // Function to read in myAngleSet
    if(myAnglesSet){
      std::ifstream anglesSetinput(string(myDhmcAnSetFile).c_str(), std::ios::in);
      if(!anglesSetinput)
        report << error <<" Could not open \'"<<myDhmcAnSetFile<<"\' for "<<getId()<<"."<<endr;
      Real tempangle = 1.0;
      while (anglesSetinput >> tempangle)
      {
        if(tempangle < 0 || tempangle >= 6.284)
        {
          report << error << "[DihedralHMCIntegrator::initialize] Dihedral Angle "<< tempangle
                 <<" out of range [0,6.284]."<<endr;
        }
        myAngles->push_back(tempangle);
      }
    }
    
  }

  void DihedralHMCIntegrator::run(int numTimesteps){
    report.precision(3);
    Real myOldPhi = 0;
    Real myOldKineticEnergy = 0;
    Real trialPhi = 0;
    Real angle = 0;

    //  -----------------------------------------------------------------------  //
    //  This run method calls the integrator below it to do the MD runs.  The    //
    //  length of the MD runs is passed to DHMC integrator thru its constructor.  //
    //  Do fileIO after every numTimesteps DHMC cycles.                           //
    //  -----------------------------------------------------------------------  //
    
    for(int i = 0; i < numTimesteps; i++) 
    {
        preStepModify();

        // --- This section of code should become the perturbation function
        
        // determine the dihedral to rotate
        if(myDihedralsSet){
          int randIndexMyDihedrals = int( randomNumber() * (myDihedrals->size()));
          if (randIndexMyDihedrals == static_cast<int>(myDihedrals->size())) // just incase randomNumber could be 1
            randIndexMyDihedrals = myDihedrals->size() - 1;
          myDihedralIndex = (*myDihedrals)[randIndexMyDihedrals];  
        }
        else // no dihedral set given... pick an arbitrary dihedral
        {
          int randIndexMyDihedrals = int( randomNumber() * (myTopo->dihedrals.size()));
          if (randIndexMyDihedrals == static_cast<int>(myTopo->dihedrals.size())) //just incase randomNumber could be 1
            randIndexMyDihedrals = myTopo->dihedrals.size() - 1;
          myDihedralIndex = randIndexMyDihedrals;   
        }

        saveValues();
        myOldPhi = computePhiDihedral(myTopo,myPositions,myDihedralIndex);
        myOldKineticEnergy = kineticEnergy(myTopo,myVelocities);
        // determine the number of radians to rotate the dihedral
        if(myAnglesSet){
          unsigned int angleindex = int ( randomNumber() * (myAngles->size()));
          if (angleindex == myAngles->size())
            angleindex = myAngles->size() - 1;
          angle = (*myAngles)[angleindex];
          report << debug(1) << "Trial angle: " << angle << endr;
          if ( (angle - myOldPhi) >= 0 )
            angle = angle - myOldPhi;
          else
            angle = 2*M_PI + (angle - myOldPhi);
        }
        else // no angle set given... pick an arbitrary angle value
        {
          angle = randomNumber()*2*M_PI;
        }
        rotateDihedral(myTopo, myPositions, myVelocities, myDihedralIndex, angle);
        buildMolecularMomentum(myVelocities,myTopo);

        // --- end the perturbation function section of code

        // Intializes forces (energies) for LF although no LF is done
        // Calculates forces for use in pre MD metropolis test
        next()->initialize(myTopo,myPositions,myVelocities,myEnergies);

        trialPhi = computePhiDihedral(myTopo,myPositions,myDihedralIndex);
        report << debug(1) << "Dihedral Index = " << myDihedralIndex << " OldPhi = " << myOldPhi
               << " Rotation Angle = " << angle << " TrialPhi = " << trialPhi << endr;

        //  ---------------------------------------------------------------  //
        //  metropolis "pre"Test - Filters "collision energies"              //
        //  ---------------------------------------------------------------  //

        Real probAccept = 0.0;
        Real oldPot = myOldEnergies->potentialEnergy();
        Real pot = (myEnergies->potentialEnergy());
        Real potPostHMC = 0;

        report << debug(1) << "DMCoriginalPot = " << oldPot << " DMCtrialPot = " << pot  << endr;

        if( (!isnan(pot)) && (metropolisTest(pot,oldPot,myInitialTemperature,probAccept)) )
        {
          report << debug(1) << "DMC: Accepted" << endr;
          report << debug(1) << "DMCoriginalPot = " << oldPot << " DMCtrialPot = " << pot  << endr;
        }
        else if( (!isnan(pot)) && (metropolisTest(pot/10,oldPot,myInitialTemperature,probAccept)) )
        {
          // restore the energies to their old value so that HMC makes the proper comparison
          // the starting potential energy value for HMC should be oldPot
          //(*myEnergies)   = (*myOldEnergies);

          randomVelocity(myInitialTemperature, myTopo,myVelocities);
          //buildMolecularMomentum(myVelocities,myTopo);

          next()->run(myCycleLength);

          probAccept = 0.0;
          potPostHMC = (myEnergies->potentialEnergy());

          if( !isnan(potPostHMC) && (metropolisTest( myEnergies->potentialEnergy() +
                                                     kineticEnergy(myTopo,myVelocities),
                                                     myOldEnergies->potentialEnergy() +
                                                     myOldKineticEnergy,
                                                     myInitialTemperature,
                                                     probAccept) ))
          {
            report << debug(1) << "Dw/HMC: Accepted" << endr;
            report << debug(1) << "DMCoriginalPot = " << oldPot << " Dw/HMCtrialPot = " << potPostHMC  << endr;
          }
          else
          {
            report << debug(1) << "Dw/HMC: Rejected" << endr;
            restoreValues();
            report << debug(1) << "DMCoriginalPot = " << oldPot << " Dw/HMCtrialPot " << potPostHMC  << endr;
          }
        }
        else
        {
          if (isnan(pot)) {
            report << warning << "Not a number reported for Energy. Not an error..." << endr;
          }
          report << debug(1) << "DMC: Rejected" << endr;
          restoreValues();
          report << debug(1) << "DMCoriginalPot = " << oldPot << " DMCtrialPot = " << pot  << endr;
        }

        postStepModify();
    }
  }


  MTSIntegrator* DihedralHMCIntegrator::doMake(string& , const vector<Value>& values, 
						ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new DihedralHMCIntegrator(values[0],values[1],values[2],values[3],values[4],values[5],values[6],fg,nextIntegrator);
  }

  void DihedralHMCIntegrator::getParameters(vector< Parameter> &parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("temperature", Value(myInitialTemperature,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("dihedralsSet",Value(myDihedralsSet),false,Text("If false a dihedral is selected randomly")));
    parameters.push_back(Parameter("dhmcDiSetFile",Value(myDhmcDiSetFile,ConstraintValueType::NotEmpty())));
    parameters.push_back(Parameter("anglesSet",Value(myAnglesSet),false,Text("If false the dihedral angle is selected randomly")));
    parameters.push_back(Parameter("dhmcAnSetFile",Value(myDhmcAnSetFile,ConstraintValueType::NotEmpty())));
  }

  bool DihedralHMCIntegrator::metropolisTest(Real newEnergy, Real oldEnergy, Real theTemperature, Real &acceptProb) {

    const Real deltaE = newEnergy - oldEnergy;

    if(deltaE < 0) {
      acceptProb = 1; 
      return true;
    } 

    acceptProb = exp(-deltaE / (Constant::BOLTZMANN * theTemperature));
    if(randomNumber() < acceptProb)
      return  true;

    return false;
  }

  void DihedralHMCIntegrator::saveValues() {
    // myOldPositions to have the val of the new positions
    (*myOldPositions)  = (*myPositions);
    (*myOldVelocities) = (*myVelocities);
    (*myOldEnergies)   = (*myEnergies);
  }


  void DihedralHMCIntegrator::restoreValues(){
    // *myPositions has to be reset to myOldPositions for the next HMC cycle
    (*myPositions)  = (*myOldPositions);
    (*myVelocities) = (*myOldVelocities);
    (*myEnergies)   = (*myOldEnergies);
    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);
  }

}

