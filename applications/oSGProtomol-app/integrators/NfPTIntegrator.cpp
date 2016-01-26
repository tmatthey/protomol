//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for NfPT dynamics                   //
//                                                                           //
//  Unless modified, this integrator uses the molecular virial to control    //
//  the pressure (barostat modifier is in oSGIntegrator.cpp). -- TIM         //
//  -----------------------------------------------------------------------  //

#include "NfPTIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierPostForceBarostat.h"
#include "oSGModifierPostForceChemostat.h"
#include "ModifierPostForceThermostat.h"
#include "ModifierPreForceBarostat.h"
#include "oSGModifierPreForceChemostat.h"
#include "ModifierPreForceThermostat.h"
#include "ModifierNPTShake.h"
#include "ModifierNPTRattle.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Keyword.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const string NfPTIntegrator::keyword("NfPTVerlet");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Default or empty constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  NfPTIntegrator::NfPTIntegrator(): oSGIntegrator() {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  NfPTIntegrator::NfPTIntegrator(Real timestep,
                                 unsigned int numComp,
			         Real temperature,
                                 Real pressure,
                                 Real tauT,
			         Real tauV,
                                 Real tauP,
			         Real tauD,
                                 Real tauL,
                                 Real MuTemp,
                                 ForceGroup *overloadedForces)
    : oSGIntegrator(timestep, numComp, temperature, pressure, tauT, tauV,
                    tauP, tauD, tauL, MuTemp,overloadedForces)  {}

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Chemostat -- prior to force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::PreForceChemostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    //  Calculate kT ln (V)
    Real kTlnV;
    if (Insert) kTlnV = kbT * log(myVolume);
    else kTlnV = -kbT * log(myVolume);

    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;
    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;

    //  Advance the chemostat velocity
    myLambdaVel += (myTargetMu + kTlnV - (*myEnergies).deltaMu()) * halfDeltaT / Qd;

    //  Advance the chemostat velocity
    myLambdaVel *= exp(-myEtaLambdaVel * halfDeltaT);

    //  Advance the chemostat
    myLambda += myLambdaVel * 2. * halfDeltaT;

    //  If Lambda should become negative, reset it to be just above zero,
    //  and if Lambda should become larger than myNumStages, reset it to be just under myNumStages.
    if (myLambda < 0) myLambda = 0.000001;
    else if (myLambda > myNumStages) myLambda = static_cast<Real>(myNumStages) - 0.000001;
  
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Chemostat -- after force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::PostForceChemostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    //  Calculate kT ln (V)
    Real kTlnV;
    if (Insert) kTlnV = kbT * log(myVolume);
    else kTlnV = -kbT * log(myVolume);

    //  Advance the chemostat velocity
    myLambdaVel *= exp(-myEtaLambdaVel * halfDeltaT);
    
    //  Advance the chemostat velocity
    myLambdaVel += (myTargetMu + kTlnV - (*myEnergies).deltaMu()) * halfDeltaT / Qd;

    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;
    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;

#ifdef DEBUG_CHEMOSTAT
    report.precision(6);
    report << hint << "*******************************************************" << endr;
    report << hint << "Molecule being transformed: " << T+1 << ", current step: " << NumTransSteps << endr;
    report << hint << "Current identity is: " << OldType << ", being inserted? " << Insert << endr;
    report << hint << "TargetMu = " << myTargetMu + kTlnV << endr;
    report << hint << "Lambda = " << myLambda << endr;
    report << hint << "Drift Current DeltaMu = " << (*myEnergies).deltaMu() << endr;
    report << hint << "LennardJones DeltaMu = " << (*myEnergies)[ScalarStructure::LENNARDJONES_DELTAMU] << endr;
    report << hint << "Coulomb DeltaMu = " << (*myEnergies)[ScalarStructure::COULOMB_DELTAMU] << endr;
    report << hint << "LambdaVel = " << myLambdaVel << endr;
    report << hint << "*******************************************************" << endr;
#endif

    //  New chemostat kinetic energy
    Real pLambda = (Qd * 0.5) * (myLambdaVel * myLambdaVel);
    //  New chemostat potential energy
    Real VLambda = -myLambda * myTargetMu;
    //  New chemostat thermostat kinetic energy 
    Real pEtaL = (Ql * 0.5) * (myEtaLambdaVel * myEtaLambdaVel);
    //  New chemostat thermostat potential energy
    Real VEtaL = myEtaLambda * Constant::BOLTZMANN * myMuTemp;
        
    //  Add the energy from the extended system chemostat.
    (*myEnergies)[ScalarStructure::INTEGRATOR] += pLambda + VLambda + pEtaL + VEtaL;
    //  add to the average chemical potential difference
    AveDeltaMu += (*myEnergies).deltaMu();

    //  add to the number of timesteps used to transform this molecule
    NumTransSteps++;

    //  compute the chemostat temperature and add to the average
    LambdaT += (Qd * myLambdaVel * myLambdaVel) / Constant::BOLTZMANN;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  doHalfKick().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::doHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // calculate the COM and momentum of each molecule
    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);

    // ---------------------------------------------------------------------
    //  Do the first update of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {
      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to the thermostat and barostat forces.
        (*myVelocities)[atom] *= exp(-(mass / myTopo->molecules[i].mass * (1.0 + 3.0/myNumFree)
				       * myEpsilonVel + myEtaVel ) * halfDeltaT);

        //  Add to the new momentum of this molecule
        Momentum += (*myVelocities)[atom] * mass;

      } //  end loop over atoms

      //  Store the updated molecular momentum
      myTopo->molecules[i].momentum = Momentum;
    } //  end loop over molecules

    // -------------------------------------------------------------------------
    //  Do the second and third updates of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to barostat force.
        (*myVelocities)[atom] -= (myTopo->molecules[i].momentum - (*myVelocities)[atom] * mass)
	  * (1.0 + 3.0/myNumFree) * (1.0 / myTopo->molecules[i].mass)
	  * myEpsilonVel * halfDeltaT;

        //  Advance the velocities due to atomic forces.
        (*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;
      }  // end loop over atoms
    } //  end loop over molecules
  }  //  End doHalfKick().

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  doDrift().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::doDrift() {

    //  Timestep.  Units: (fs)
    const Real deltaT = getTimestep();

    // -------------------------------------------------------------------------
    //  Do the first update of the atom positions
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the positions due to change in box length.
        (*myPositions)[atom] *= exp((mass / myTopo->molecules[i].mass)
				    * myEpsilonVel * deltaT);
      }  // end loop over atoms
    } //  end loop over molecules

    // -------------------------------------------------------------------------
    //  Do the second and final update of the positions
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Temporary storage element for the updated molecular COM
      Vector3D COM(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the positions due to box volume or barostat (similar to COM shifting).
        (*myPositions)[atom] += (myTopo->molecules[i].position - (*myPositions)[atom]
				 * (mass / myTopo->molecules[i].mass)) * deltaT * myEpsilonVel;

        //  Advance the positions due to velocity.
        (*myPositions)[atom] += (*myVelocities)[atom] * deltaT * Constant::INV_TIMEFACTOR;

        //  Add to the new COM of this molecule
        COM += (*myPositions)[atom] * mass;
      }  //  end loop over atoms

      //  Store the updated molecular COM
      myTopo->molecules[i].position = COM / (myTopo->molecules[i].mass);
    } //  end loop over molecules

    //  update the COM of each molecule
    buildMolecularCenterOfMass(myPositions,myTopo);
  }  // end doDrift()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  do2ndHalfkick()
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::do2ndHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // -----------------------------------------------------------------------------
    //  Do the first update of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to atomic forces.
        (*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;

        //  Add to the new momentum of this molecule
        Momentum += (*myVelocities)[atom] * mass;
      }  //  end loop over atoms

      //  Store the updated molecular momentum
      myTopo->molecules[i].momentum = Momentum;
    } //  end loop over molecules

    // ------------------------------------------------------------------------------
    //  Do the second and third updates of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to barostat force.
        (*myVelocities)[atom] -= (myTopo->molecules[i].momentum - (*myVelocities)[atom] * mass)
	  * (1.0 + 3.0/myNumFree) * (1.0 / myTopo->molecules[i].mass)
	  * myEpsilonVel * halfDeltaT;

        //  Advance the velocities due to the thermostat and barostat forces.
        (*myVelocities)[atom] *= exp(-(mass / myTopo->molecules[i].mass * (1.0 + 3.0/myNumFree)
          * myEpsilonVel + myEtaVel ) * halfDeltaT);

        //  Add to the new momentum of this molecule
        Momentum += (*myVelocities)[atom] * mass;
      }  // end loop over atoms

      //  Store the updated molecular momentum
      myTopo->molecules[i].momentum = Momentum;
    } //  end loop over molecules
  }  // end do2ndHalfKick
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  initialize().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::initialize(GenericTopology *topo,
                                 Vector3DBlock *positions,
                                 Vector3DBlock *velocities,
                                 ScalarStructure *energies) {

    oSGIntegrator::initialize(topo, positions, velocities, energies);

    //  store the total number of atoms and molecules
    NumAtoms = myTopo->atoms.size();
    NumMols = myTopo->molecules.size();
    myNumFree = myTopo->degreesOfFreedom;
    
    //  determine the # of molecules of each component in the mixture
    // initialize the # of component to zero
    N.resize(myNumComp);
    myTopo->iSGNumMols.resize(myNumComp);
    for (unsigned int C=0; C < myNumComp; C++) {
      N[C] = 0;
      myTopo->iSGNumMols[C] = 0;
    }

    // using the type ID of each molecule, compute the # of molecules
    // of each component
    for (unsigned int M=0; M < NumMols; M++) {
      N[myTopo->molecules[M].type]++;
      myTopo->iSGNumMols[myTopo->molecules[M].type]++;
    }

    //  Initialize integrator variables to something sane.
    if (T == -1) {
      myEpsilonVel = 0.;
      myEta        = 0.;
      myEtaV       = 0.;
      myEtaLambda  = 0.;
      myEtaVel     = 0.;
      myEtaVolVel  = 0.;
      myLambdaVel  = 0.;
      myEtaLambdaVel = 0.;

      //  Initialize myVolume
      myVolume = myTopo->getVolume(*positions);
    }
    else {
      // Initialize the volume
      Real fac = myVolume / myTopo->getVolume();
      myTopo->rescaleVolume(fac);
      uncache();

      // Initialize the transforming molecule identities
      myTopo->molecules[T].lambda = myLambda;
      myTopo->molecules[T].type = OldType;
      if (Insert) {
        myTopo->molecules[T].newtype = OldType;
        
        // since this molecule has not been completely inserted, we must subract one from the # of molecules
        N[myTopo->molecules[T].type]--;
        myTopo->iSGNumMols[myTopo->molecules[T].type]--;
      }
      else myTopo->molecules[T].newtype = OldType + 1;
      
      // Compute the residual chemical potential
      // Mu = kT*ln(fi/kT * V/Ni), the fugacity part was already added by ModifierOSG::ReadXSCs,
      // so here we just add the N or N+1 term since the volume term is added by the chemostat modifier
      if (Insert) myTargetMu -= kbT * log (static_cast<Real>(N[OldType]+1));
      else myTargetMu += kbT * log (static_cast<Real>(N[OldType]));
      myTargetMu /= static_cast<Real>(myNumStages);      
    }

    //  Compute the fixed barostat mass and thermostat masses.
    Qo = myNumFree * kbT * (myTauT * myTauT);
    Qv = kbT * (myTauV * myTauV);
    W = (myNumFree + 3) * kbT * (myTauP * myTauP);
    Qd = (Constant::BOLTZMANN * myMuTemp) * (myTauD * myTauD);
    Ql = (Constant::BOLTZMANN * myMuTemp) * (myTauL * myTauL);

    initializeForces();
    myEnergies->output(false);

  } // end function initialize

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NPT RATTLE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifier* NfPTIntegrator::createRattleModifier(Real eps, int maxIter){
    return (new ModifierNPTRattle<NfPTIntegrator>(eps, maxIter, this));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NPT SHAKE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifier* NfPTIntegrator::createShakeModifier(Real eps, int maxIter){
    return (new ModifierNPTShake<NfPTIntegrator>(eps, maxIter, this));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Add the modifiers
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::addModifierAfterInitialize() {
    adoptPreStepModifier(new oSGModifierPreForceChemostat(this,3));
    adoptPreStepModifier(new ModifierPreForceBarostat<NfPTIntegrator>(this,2));
    adoptPreStepModifier(new ModifierPreForceThermostat<NfPTIntegrator>(this,1));
    adoptPostStepModifier(new ModifierPostForceThermostat<NfPTIntegrator>(this,3));
    adoptPostStepModifier(new ModifierPostForceBarostat<NfPTIntegrator>(this,2));
    adoptPostStepModifier(new oSGModifierPostForceChemostat(this,1));  
    STSIntegrator::addModifierAfterInitialize();
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // getParameters().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void NfPTIntegrator::getParameters(vector< Parameter >& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("components", Value(myNumComp)));
    parameters.push_back(Parameter("temperature", Value(myTargetTemp)));
    parameters.push_back(Parameter("pressure", Value(myTargetPres)));
    parameters.push_back(Parameter("tauT", Value(myTauT)));
    parameters.push_back(Parameter("tauV", Value(myTauV)));
    parameters.push_back(Parameter("tauP", Value(myTauP)));
    parameters.push_back(Parameter("tauD", Value(myTauD)));
    parameters.push_back(Parameter("tauL", Value(myTauL)));
    parameters.push_back(Parameter("lambdaTemp", Value(myMuTemp)));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // doMake().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  STSIntegrator* NfPTIntegrator::doMake(string&, const vector<Value>& values,ForceGroup* fg)const{
    return new NfPTIntegrator(values[0], values[1], values[2], values[3], values[4],
			      values[5], values[6], values[7], values[8], values[9], fg );
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // finalXSCs
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  XSC & NfPTIntegrator::getXSC() const {

    // create an XSC object to store the final xsc values
    XSC *xsc = new XSC;

    // set the integrator type
    xsc->simType = "NfPTVerlet";

    // for chemostat
    xsc->Lambda = myLambda;
    xsc->Lambda_vel = myLambdaVel;
    xsc->myMolecule = T+1;
    xsc->old_type = OldType;
    if (Insert) xsc->new_type = OldType;
    else xsc->new_type = OldType + 1;

    // for thermostats
    xsc->Eta = myEta;
    xsc->EtaVol = myEtaV;
    xsc->Eta_vel = myEtaVel;
    xsc->EtaVol_vel = myEtaVolVel;
    xsc->EtaLambda = myEtaLambda;
    xsc->EtaLambda_vel = myEtaLambdaVel;

    // for barostat
    xsc->Vol = myVolume;
    xsc->Epsilon_vel = myEpsilonVel;

    return (*xsc);
  }

}
