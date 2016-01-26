#include "StandardIntegrator.h"
#include "Report.h"
#include "Parallel.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"

namespace ProtoMol {
  //_________________________________________________________________ StandardIntegrator
  StandardIntegrator::StandardIntegrator() :
      Integrator(), myPreviousIntegrator(NULL) {}

  StandardIntegrator::StandardIntegrator(ForceGroup* forceGroup)
          : Integrator(forceGroup),myPreviousIntegrator(NULL) {}

  void StandardIntegrator::run(int numTimesteps){
    for(int i = 0; i < numTimesteps; i++){
      preStepModify();
      doHalfKick();
      doDriftOrNextIntegrator();
      calculateForces();
      doHalfKick();
      postStepModify();
    }
  }

  void StandardIntegrator::initialize(GenericTopology *topo,
				      Vector3DBlock   *positions, 
				      Vector3DBlock   *velocities, 
				      ScalarStructure *energies){
    Integrator::initialize(topo, positions, velocities,energies);
    //Report::report <<"[StandardIntegrator::initialize]"<<Report::endr;
  }

  void StandardIntegrator::initializeForces() {
    //Report::report <<"[StandardIntegrator::initializeForces]"<<Report::endr;
    addModifierBeforeInitialize();
    calculateForces();
    addModifierAfterInitialize();
  } 

  void StandardIntegrator::calculateForces(){
    //Report::report <<"[StandardIntegrator::calculateForces]"<<Report::endr;

    //  Save current value of potentialEnergy().
    myPotEnergy = myEnergies->potentialEnergy();

    myForces->zero();
    preForceModify();
    if(!anyMediForceModify()){
      Parallel::distribute(myEnergies,myForces);
      //Report::report <<"calculateForces dist"<<Report::endr;
    }
    myForcesToEvaluate->evaluateSystemForces(myTopo, myPositions, myForces, myEnergies);
    mediForceModify();
    myForcesToEvaluate->evaluateExtendedForces(myTopo, myPositions, myVelocities, myForces, myEnergies);
    if(!anyMediForceModify()){
      Parallel::reduce(myEnergies,myForces);
      //Report::report <<"calculateForces red"<<Report::endr;
    }
    postForceModify();

    //  Compute my potentialEnergy as the difference before/after the call to
    //  calculateForces().
    myPotEnergy = myEnergies->potentialEnergy() - myPotEnergy;

  }

  void StandardIntegrator::doHalfKick(){

    Real h = 0.5 * getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = myPositions->size();

    updateBeta(h);

    for( unsigned int i = 0; i < count; ++i) {
      (*myVelocities)[i] += (*myForces)[i] * h / myTopo->atoms[i].scaledMass;
    }  

    buildMolecularMomentum(myVelocities,myTopo);

  }

  void StandardIntegrator::doKick(){

    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = myPositions->size();

    updateBeta(h);

    for( unsigned int i = 0; i < count; ++i) {
      (*myVelocities)[i] += (*myForces)[i] * h / myTopo->atoms[i].scaledMass;
    }  

    buildMolecularMomentum(myVelocities,myTopo);

  }

  Integrator* StandardIntegrator::previous(){
    return myPreviousIntegrator;
  }

  const Integrator* StandardIntegrator::previous() const{
    return myPreviousIntegrator;
  }


}  
