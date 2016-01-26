#include "DMDLeapfrogIntegrator.h"
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

  //_________________________________________________________________ DMDLeapfrogIntegrator

  const string DMDLeapfrogIntegrator::keyword("DMDLeapfrog");

  DMDLeapfrogIntegrator::DMDLeapfrogIntegrator()
    :STSIntegrator(),
     myDissipativeForces(NULL),
     myRandomForces(NULL),
     myVhat(NULL),
     myGamma(-1.0),
     myTemperature(-1.0),
     myNumIter(-1),
     mySigma(0.0),
     mySeed(-1)
  {}

  DMDLeapfrogIntegrator::DMDLeapfrogIntegrator(Real timestep,
					       int numIter,
					       Real gamma,
					       Real temperature,
					       int seed,
					       ForceGroup *overloadedForces)
    :STSIntegrator(timestep,overloadedForces),
     myDissipativeForces(new Vector3DBlock),
     myRandomForces(new Vector3DBlock),
     myVhat(new Vector3DBlock),
     myGamma(gamma*0.001),
     myTemperature(temperature),
     myNumIter(numIter),
     mySigma(sqrt(2*myGamma*myTemperature*Constant::BOLTZMANN)),
     mySeed(seed)
  {}


  DMDLeapfrogIntegrator::~DMDLeapfrogIntegrator()
  {
    if(myDissipativeForces != NULL)
      delete myDissipativeForces;
    if(myRandomForces != NULL)
      delete myRandomForces;
    if(myVhat != NULL)
      delete myVhat;
  }

  void DMDLeapfrogIntegrator::initialize(GenericTopology *topo,
					 Vector3DBlock *positions,
					 Vector3DBlock *velocities,
					 ScalarStructure *energies) 
  {  
    STSIntegrator::initialize(topo,positions,velocities,energies);
    initializeForces();

    myDissipativeForces->zero(positions->size());
    myRandomForces->zero(positions->size());
    myVhat->zero(positions->size());
    calculateDissipativeAndRandomForces();
  }


  void DMDLeapfrogIntegrator::doHalfKick() 
  {

    Real f = 0.5*getTimestep()/Constant::TIMEFACTOR;
    Real g = 0.5*sqrt(2*f);
    Vector3DBlock myTempCoordBlock(myPositions->size());
    myTempCoordBlock.zero();
    myTempCoordBlock.intoWeightedAdd(f, *myForces); 
    myTempCoordBlock.intoWeightedAdd(f, *myDissipativeForces); 
    myTempCoordBlock.intoWeightedAdd(g, *myRandomForces); 
    // We only compute the doHalfKick and doDrift in our range
    for(unsigned int i=0;i<myPositions->size();i++)
      (*myVelocities)[i] += myTempCoordBlock[i]/myTopo->atoms[i].scaledMass;

    buildMolecularMomentum(myVelocities,myTopo);
  }


  void DMDLeapfrogIntegrator::doHalfKickVhat() 
  {

    Real f = 0.5*getTimestep()/Constant::TIMEFACTOR;
    Real g = 0.5*sqrt(2*f);
    Vector3DBlock myTempCoordBlock(myPositions->size());
    myTempCoordBlock.zero();
    myTempCoordBlock.intoWeightedAdd(f, *myForces);
    myTempCoordBlock.intoWeightedAdd(g, *myRandomForces); 
    // We only compute the doHalfKick and doDrift in our range
    for(unsigned int i=0;i<myPositions->size();i++)
      (*myVhat)[i] = (*myVelocities)[i] + 
	myTempCoordBlock[i]/myTopo->atoms[i].scaledMass;

  }

  void DMDLeapfrogIntegrator::doHalfKickIterate() {

    Real f = 0.5*getTimestep()/Constant::TIMEFACTOR;
    Vector3DBlock myTempCoordBlock(myPositions->size());
    Vector3DBlock vHat(myPositions->size());
    vHat.zero();
    vHat.intoAdd(*myVelocities);
    for(int j=0;j<myNumIter;j++){
      myTempCoordBlock.zero();
      myTempCoordBlock.intoWeightedAdd(f, *myDissipativeForces); 
      // We only compute the doHalfKick and doDrift in our range
      for(unsigned int i=0;i<myPositions->size();i++)
	(*myVelocities)[i] = (*myVhat)[i] + 
	  myTempCoordBlock[i]/myTopo->atoms[i].scaledMass;
      calculateDissipativeForces();
    }
    buildMolecularMomentum(myVelocities,myTopo);
  
  }


  //Here I did a little different than the litterature: 
  //"Towards better integrators for dissipative particle dynamics simulations"
  //by Gerhard Besold, etc. In their paper, they defined rij=ri-rj and vij=vi-vj,
  //which does not conform with mathematical convention. I changed them into
  //rij=rj-ri, and vij=vj-vi. Thus the dissipitive pairwise force F_ij,
  //F_ij = - gamma*weight^2*(v_ij \dot e_ij) e_ij, where e_ij is the unit vector.
  //This is the force exerted on particle "j" by particle "i".
  void DMDLeapfrogIntegrator::calculateDissipativeForces()
  {
    myDissipativeForces->zero();  
    for (unsigned int i = 0; i < myTopo->angles.size(); i++){
      int a1 = myTopo->angles[i].atom1;
      int a2 = myTopo->angles[i].atom2;
      int a3 = myTopo->angles[i].atom3;

      const Vector3D& posi1 = (*myPositions)[a1];
      const Vector3D& posi2 = (*myPositions)[a2];
      const Vector3D& posi3 = (*myPositions)[a3];
      const Vector3D& vel1 = (*myVelocities)[a1];
      const Vector3D& vel2 = (*myVelocities)[a2];
      const Vector3D& vel3 = (*myVelocities)[a3];

      //now handle atoms 1 and 2
      Vector3D unitVec = posi2 - posi1; // now it is only a vector
      unitVec.normalize(); 
      // now the vector is unit vector, the dist is the length of the original vector
      // We assume the weight, w1, and w2 to be 1 always, thus the following still holds true
      // w1 = w2^2;
      Real coeff = -myGamma*((vel2-vel1).dot(unitVec));
      (*myDissipativeForces)[a2] += unitVec * coeff;
      (*myDissipativeForces)[a1] -= unitVec * coeff;

      //now handle atoms 2 and 3
      unitVec = posi3 - posi2; // now it is only a vector
      unitVec.normalize(); 

      coeff = -myGamma*((vel3-vel2).dot(unitVec));
      (*myDissipativeForces)[a3] += unitVec * coeff;
      (*myDissipativeForces)[a2] -= unitVec * coeff;

      //now handle atoms 1 and 3
      unitVec = posi3 - posi1; // now it is only a vector
      unitVec.normalize(); 
      coeff = -myGamma*((vel3-vel1).dot(unitVec));
      (*myDissipativeForces)[a3] += unitVec * coeff;
      (*myDissipativeForces)[a1] -= unitVec * coeff;
    }
  }


  void DMDLeapfrogIntegrator::calculateDissipativeAndRandomForces()
  {
    Real sdv = 1;
    myRandomForces->zero();  
    myDissipativeForces->zero();  
    for (unsigned int i = 0; i < myTopo->angles.size(); i++){
      int a1 = myTopo->angles[i].atom1;
      int a2 = myTopo->angles[i].atom2;
      int a3 = myTopo->angles[i].atom3;

      const Vector3D& posi1 = (*myPositions)[a1];
      const Vector3D& posi2 = (*myPositions)[a2];
      const Vector3D& posi3 = (*myPositions)[a3];
      const Vector3D& vel1 = (*myVelocities)[a1];
      const Vector3D& vel2 = (*myVelocities)[a2];
      const Vector3D& vel3 = (*myVelocities)[a3];

      //now handle atoms 1 and 2
      Vector3D unitVec = posi2 - posi1; // now it is only a vector
      unitVec.normalize(); 
      // now the vector is unit vector, the dist is the length of the original vector
      // We assume the weight, w1, and w2 to be 1 always, thus the following still holds true
      // w1 = w2^2;

      Real coeff = -myGamma*((vel2-vel1).dot(unitVec));
      (*myDissipativeForces)[a2] += unitVec * coeff;
      (*myDissipativeForces)[a1] -= unitVec * coeff;
      Real randNum = randomGaussian(sdv, mySeed);
      coeff = mySigma*randNum;
      (*myRandomForces)[a2] += unitVec * coeff;
      (*myRandomForces)[a1] -= unitVec * coeff;

      //now handle atoms 2 and 3
      unitVec = posi3 - posi2; // now it is only a vector
      unitVec.normalize(); 

      coeff = -myGamma*((vel3-vel2).dot(unitVec));
      (*myDissipativeForces)[a3] += unitVec * coeff;
      (*myDissipativeForces)[a2] -= unitVec * coeff;
      randNum = randomGaussian(sdv, mySeed);
      coeff = mySigma*randNum;
      (*myRandomForces)[a3] += unitVec * coeff;
      (*myRandomForces)[a2] -= unitVec * coeff;

      //now handle atoms 1 and 3
      unitVec = posi3 - posi1; // now it is only a vector
      unitVec.normalize(); 

      coeff = -myGamma*((vel3-vel1).dot(unitVec));
      (*myDissipativeForces)[a3] += unitVec * coeff;
      (*myDissipativeForces)[a1] -= unitVec * coeff;
      randNum = randomGaussian(sdv, mySeed);
      coeff = mySigma*randNum;
      (*myRandomForces)[a3] += unitVec * coeff;
      (*myRandomForces)[a1] -= unitVec * coeff;
    }
  }


  void DMDLeapfrogIntegrator::run (int numTimesteps) {
    for (int i = 0; i < numTimesteps; i++) {
      preStepModify();
      doHalfKick();
      doDriftOrNextIntegrator();
      calculateForces(); // only calculate the conservative forces.
      calculateDissipativeAndRandomForces();
      doHalfKickVhat();
      doHalfKickIterate();//iteratively calculate the dissipative force and velocities.
      postStepModify();
    }

  }

  STSIntegrator* DMDLeapfrogIntegrator::doMake(string& , const vector<Value>& values,ForceGroup* fg)const{
    return new DMDLeapfrogIntegrator(values[0],values[1],values[2],values[3],values[4],fg);

  }

  void DMDLeapfrogIntegrator::getParameters(std::vector<Parameter>& parameters) const{
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("iterations",Value(myNumIter,ConstraintValueType::NotNegative()),2));
    parameters.push_back(Parameter("gamma",Value(myGamma*1000,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("temperature",Value(myTemperature,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::Positive()),1234));
  }
}


