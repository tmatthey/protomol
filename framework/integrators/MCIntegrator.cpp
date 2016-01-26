#include "MCIntegrator.h"
#include "Report.h"
#include "GenericTopology.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "ForceGroup.h"

using namespace ProtoMol::Report;
using std::vector;

namespace ProtoMol {
  //____________________________________________________________ MCIntegrator

  MCIntegrator::MCIntegrator()    
    : MTSIntegrator(),
      myInitialTemperature(-1),
      myOldKineticEnergy(-1),
      myOldPositions(NULL),
      myOldVelocities(NULL),
      myOldEnergies(NULL){
  }

  //  ---------------------------------------------------------------------  //

  MCIntegrator::MCIntegrator(int cycles,
			     Real initialTemperature,
			     ForceGroup *overloadedForces,
			     StandardIntegrator *nextIntegrator)
    
    : MTSIntegrator(cycles,overloadedForces,nextIntegrator),
      myInitialTemperature(initialTemperature),
      myOldKineticEnergy(-1),
      myOldPositions(new Vector3DBlock),
      myOldVelocities(new Vector3DBlock),
      myOldEnergies(new ScalarStructure) {
    if(overloadedForces->anyForces()){
      report << error << "Monte Carlo integrators do not expect any forces!"<<endr;
    }

  }
  
  //  ---------------------------------------------------------------------  //

  MCIntegrator::~MCIntegrator(){
    if(myOldPositions != NULL)
      delete myOldPositions;
    if(myOldVelocities != NULL)
      delete myOldVelocities;
    if(myOldEnergies != NULL)
      delete myOldEnergies;
  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::initialize(GenericTopology *topo,
				Vector3DBlock *positions,
				Vector3DBlock *velocities,
				ScalarStructure *energies){

    MTSIntegrator::initialize(topo,positions,velocities,energies);

    //Report::report <<"[MCIntegrator::initialize]"<<Report::endr;

    myOldPositions->resize(positions->size());
    myOldVelocities->resize(velocities->size());


  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::run(int numTimesteps) {

    saveValues();

    for(int i = 0; i < numTimesteps; i++) {

      preStepModify();

      perturbSystem();

      //  The metropolis test is based in part on the change in KE from the new
      //  velocities.
      myOldKineticEnergy = kineticEnergy(myTopo,myVelocities);

      walk(myCycleLength); 

      if(metropolisTest()) {
	saveValues();

        report << debug(1) << "Move accepted" << endr;

      }
      else {
	restoreValues();

        report << debug(1) << "Move rejected" << endr;

      } 
      postStepModify();
    }  

  }

  //  ---------------------------------------------------------------------  //

  bool MCIntegrator::metropolisTest() {

    return metropolisTest( myEnergies->potentialEnergy() +
                           kineticEnergy(myTopo,myVelocities),
			   myOldEnergies->potentialEnergy() +
                           myOldKineticEnergy );
  }

  //  ---------------------------------------------------------------------  //

  bool MCIntegrator::metropolisTest( Real newValue, Real oldValue ) {

    Real diff = newValue - oldValue;

    if( diff < 0 )
      return true;

    Real acceptProb = exp( -diff /
                           ( Constant::BOLTZMANN * getInitialTemperature() ) );

    report << debug(2) << "Metropolis prob: = " << acceptProb << endr;

    return( randomNumber() < acceptProb );

  }

  //  ---------------------------------------------------------------------  //

  bool MCIntegrator::metropolisTest( Real newValue, Real oldValue,
                                     Real & acceptProb ) {

      Real diff = newValue - oldValue;

      if( diff < 0 ) {

          acceptProb = 1.0;

          return true;

      }

      acceptProb = exp( -diff / ( Constant::BOLTZMANN * getInitialTemperature() ) );

      return( randomNumber() < acceptProb );

  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::getParameters(vector< Parameter> &parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter( "temperature", Value( myInitialTemperature,
                    ConstraintValueType::NotNegative() ) ) );
  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::saveValues() {
    (*myOldPositions)  = (*myPositions);
    (*myOldVelocities) = (*myVelocities);
    (*myOldEnergies)   = (*myEnergies);
    next()->saveForces();
  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::restoreValues() {
    (*myPositions)  = (*myOldPositions);
    (*myEnergies)   = (*myOldEnergies);
    (*myVelocities) = (*myOldVelocities);
    next()->restoreForces();
    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);
  }

  //  ---------------------------------------------------------------------  //

  void MCIntegrator::walk(int steps){
    myNextIntegrator->run(steps);
  }

}

