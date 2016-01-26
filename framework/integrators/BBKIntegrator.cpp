#include "BBKIntegrator.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"


using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
  //________________________________________________________ BBKIntegrator

  const string BBKIntegrator::keyword("BBK");


  BBKIntegrator::BBKIntegrator():
    STSIntegrator(),
    myLangevinTemperature(-1.0),
    myGamma(-1.0), 
    mySeed(-1){}

  BBKIntegrator::BBKIntegrator(Real timestep,
                   Real langevinTemperature,
                   Real gamma,
                   int seed,
                   ForceGroup *overloadedForces):
    STSIntegrator(timestep,overloadedForces),
    myLangevinTemperature(langevinTemperature),
    myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)), 
    mySeed(seed){
  }

  void BBKIntegrator::initialize(GenericTopology *topo,
                 Vector3DBlock *positions,
                 Vector3DBlock *velocities,
                 ScalarStructure *energies){
  
    STSIntegrator::initialize(topo,positions,velocities,energies);
    initializeForces();
  }

  void BBKIntegrator::run (int numTimesteps) {

    for (int i = 0; i < numTimesteps; i++) {
      preStepModify();
      doFirstHalfKick();
      doDriftOrNextIntegrator();
      calculateForces();
      doSecondHalfKick();
      postStepModify();
    }  
  }

  void BBKIntegrator::doFirstHalfKick(){
    // Our range. This also assumed by the methods in Parallel!!!
    Real dt = getTimestep()/Constant::TIMEFACTOR;
    //Thermostat factor
    Real forceConstant = 2*Constant::BOLTZMANN*myLangevinTemperature*myGamma;
    // We only compute the doKick & doDrift on our range
    for(unsigned int i=0;i<myPositions->size();i++){
      (*myVelocities)[i] *= 1.0-0.5*dt*myGamma;
      (*myVelocities)[i] += (*myForces)[i]*dt*0.5/myTopo->atoms[i].scaledMass;
      //Add random pertubation      
      Vector3D randomForce(randomGaussianNumber(mySeed),randomGaussianNumber(mySeed),randomGaussianNumber(mySeed));
      (*myVelocities)[i] += randomForce * sqrt(dt * 0.5 * forceConstant / myTopo->atoms[i].scaledMass);
    }
    buildMolecularMomentum(myVelocities,myTopo);
  }

  void BBKIntegrator::doSecondHalfKick(){
    // Our range. This also assumed by the methods in Parallel!!!
    Real dt = getTimestep()/Constant::TIMEFACTOR;
    //Thermostat factor
    Real forceConstant = 2*Constant::BOLTZMANN*myLangevinTemperature*myGamma;
    // We only compute the doKick & doDrift on our range
    for(unsigned int i=0;i<myPositions->size();i++){
      (*myVelocities)[i] += (*myForces)[i]*dt*0.5/myTopo->atoms[i].scaledMass;
      //Add random pertubation      
      Vector3D randomForce(randomGaussianNumber(mySeed),randomGaussianNumber(mySeed),randomGaussianNumber(mySeed));
      (*myVelocities)[i] += randomForce * sqrt(dt * 0.5 * forceConstant / myTopo->atoms[i].scaledMass);
      (*myVelocities)[i] /= (1.0+0.5*dt*myGamma);
    }
    buildMolecularMomentum(myVelocities,myTopo);
  }

  void BBKIntegrator::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("temperature",Value(myLangevinTemperature,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234));
  }

  STSIntegrator* BBKIntegrator::doMake(string& , const vector<Value>& values,ForceGroup* fg)const{
    return new BBKIntegrator(values[0],values[1],values[2],values[3],fg);

  }

}
