#include "LangevinLeapfrogIntegrator.h"
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
  //_________________________________________________________________ LangevinLeapfrogIntegrator

  const string LangevinLeapfrogIntegrator::keyword("LangevinLeapfrog");

  LangevinLeapfrogIntegrator::LangevinLeapfrogIntegrator():
    STSIntegrator(),
    myLangevinTemperature(-1.0),
    myGamma(-1.0), // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
    mySeed(-1){
  }

  LangevinLeapfrogIntegrator::LangevinLeapfrogIntegrator(Real timestep,
                               Real LangevinTemperature,
                               Real gamma,
                               int seed,
                               ForceGroup *overloadedForces):
    STSIntegrator(timestep,overloadedForces),
    myLangevinTemperature(LangevinTemperature),
    myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)), // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
    mySeed(seed){
  }

  void LangevinLeapfrogIntegrator::initialize(GenericTopology *topo,
                         Vector3DBlock *positions,
                         Vector3DBlock *velocities,
                         ScalarStructure *energies){
  
    STSIntegrator::initialize(topo,positions,velocities,energies);
    initializeForces();
  }

  void LangevinLeapfrogIntegrator::doDrift() {
      const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
      myPositions->intoWeightedAdd(h, *myVelocities);
      buildMolecularCenterOfMass(myPositions,myTopo);
      buildMolecularMomentum(myVelocities,myTopo);
  }

  void LangevinLeapfrogIntegrator::doHalfKick() {
    const unsigned int count = myPositions->size();
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
    const Real vdt = exp(-0.5*myGamma*dt);
    const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) );
    const Real forceConstant = 2 * Constant::BOLTZMANN * myLangevinTemperature * myGamma;

    for( int i = 0; i < count; i++ ) {
        //  Generate gaussian random numbers for each spatial direction
        Vector3D gaussRandCoord1(randomGaussianNumber(mySeed), randomGaussianNumber(mySeed), randomGaussianNumber(mySeed));
        Real mass = myTopo->atoms[i].scaledMass;
        Real sqrtFCoverM = sqrt(forceConstant / mass);
        // semi-update velocities
        (*myVelocities)[i] = (*myVelocities)[i]*vdt
                                +(*myForces)[i] * fdt / mass
                                    +gaussRandCoord1*sqrtFCoverM*ndt;
    }
    buildMolecularMomentum(myVelocities,myTopo);
  }

  void LangevinLeapfrogIntegrator::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("temperature",Value(myLangevinTemperature,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234));
  }

  STSIntegrator* LangevinLeapfrogIntegrator::doMake(string& , const vector<Value>& values,ForceGroup* fg)const{
    return new LangevinLeapfrogIntegrator(values[0],values[1],values[2],values[3],fg);

  }

}
