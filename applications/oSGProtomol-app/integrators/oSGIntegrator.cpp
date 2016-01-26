//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for iSG dynamics                    //
//                                                                           //
//  Unless modified, this integrator uses the molecular virial to control    //
//  the pressure. -- TIM                                                     //
//  -----------------------------------------------------------------------  //

#include "oSGIntegrator.h"
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

#include "mathutilities.h"
#include "Report.h"
#include "Integrator.h"
#include "XSCReader.h"
#include "AtomType.h"
#include "Array.h"
#include "Parallel.h"

//#define DEBUG_MODIFIER_ATOM
//#define DEBUG_CHEMOSTAT

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Keyword.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const string oSGIntegrator::scope("oSGIntegrator");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Default or empty constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  oSGIntegrator::oSGIntegrator(): STSIntegrator(),
				  myNumComp(0),
                                  myNumStages(0),
				  myTargetTemp(0.0),
				  myTargetPres(0.0),
                                  myMuTemp(0.0),
				  myTauT(0.0),
				  myTauV(0.0),
				  myTauP(0.0),
				  myTauD(0.0),
                                  myTauL(0.0),
                                  kbT(0.0),
                                  myLambda(0.0),
                                  NumAtoms(0),
                                  NumMols(0),
                                  myNumFree(0),
				  T(-1),
                                  thisStage(0),
				  NumTransSteps(0),
				  OldType(0),
				  Insert(true),
				  Transformed(true),
				  myVolume(0.0),
				  myEpsilonVel(0.0),
				  Qo(0.0),
				  Qv(0.0),
				  W(0.0),
				  Qd(0.0),
				  Ql(0.0),
                                  myEta(0.0),
                                  myEtaV(0.0),
                                  myEtaLambda(0.0),
                                  myEtaVel(0.0),
                                  myEtaVolVel(0.0),
                                  myLambdaVel(0.0),
                                  myEtaLambdaVel(0.0),
				  AveCQ(0.0),
				  AveCQSq(0.0),
                                  AveDeltaMu(0.0),
                                  LambdaT(0.0),
                                  TotSteps(0),
                                  numAttempts(0),
                                  fugacityFactor(0) {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  oSGIntegrator::oSGIntegrator(Real timestep,
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
    : STSIntegrator(timestep, overloadedForces),
      myNumComp(numComp),
      myNumStages(0),
      myTargetTemp(temperature),
      myTargetPres(pressure),
      myMuTemp(MuTemp),
      myTauT(tauT), 
      myTauV(tauV), 
      myTauP(tauP),
      myTauD(tauD),
      myTauL(tauL),
      kbT(temperature * Constant::BOLTZMANN),
      myLambda(0.0),
      NumAtoms(0),
      NumMols(0),
      myNumFree(0),
      T(-1),
      thisStage(1),
      NumTransSteps(0),
      OldType(0),
      Insert(true),
      Transformed(true),
      myVolume(0.0),
      myEpsilonVel(0.0),
      Qo(0.0),
      Qv(0.0),
      W(0.0),
      Qd(0.0),
      Ql(0.0),
      myEta(0.0),
      myEtaV(0.0),
      myEtaLambda(0.0),
      myEtaVel(0.0),
      myEtaVolVel(0.0),
      myLambdaVel(0.0),
      myEtaLambdaVel(0.0),
      AveCQ(0.0),
      AveCQSq(0.0),
      AveDeltaMu(0.0),
      LambdaT(0.0),
      TotSteps(0),
      numAttempts(0),
      fugacityFactor(1.43836263972E-5) {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Thermostat -- prior to force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::PreForceThermostat() {
    
    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();
    //  Twice the kinetic energy.  Units: (kcal / mol)
    const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);

    //  Advance the particle thermostat variable.
    myEta += myEtaVel * halfDeltaT;
    //  Advance the particle thermostat variable velocity.
    myEtaVel += (twiceKE - myNumFree * kbT) * halfDeltaT / Qo;
  
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Thermostat -- after force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::PostForceThermostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();
    //  Get the new KE.
    const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);

    //  Advance the particle thermostat variable velocity.
    myEtaVel += (twiceKE - myNumFree * kbT) * halfDeltaT / Qo;
    //  Advance the particle thermostat variable.
    myEta += myEtaVel * halfDeltaT;

    //  New particle thermostat kinetic energy.
    Real pEta = (Qo * 0.5) * (myEtaVel * myEtaVel);
    //  New particle thermostat potential energy.
    Real VEta = myEta * myNumFree * kbT;

    //  Add the energy from the extended system thermostat.
    (*myEnergies)[ScalarStructure::INTEGRATOR] += pEta + VEta;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Barostat -- prior to force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::PreForceBarostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();
    //  Twice the kinetic energy.  Units: (kcal / mol)
    const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);
    const Real MolKE = molecularKineticEnergy(myTopo, myVelocities);

    //  Advance the volume thermostat variable.
    myEtaV += myEtaVolVel * halfDeltaT;
    //  Advance the volume thermostat variable velocity.
    myEtaVolVel += (W * myEpsilonVel * myEpsilonVel - kbT) * halfDeltaT / Qv;

    //  Calculate the current molecular pressure.  Units: (bar)
    Real currentPres = computeMolecularPressure(myEnergies, myVolume, MolKE);

    //  Advance the box volume velocity.  Divide by PRESSUREFACTOR to get
    //  correct (fs)^-1 units.
    myEpsilonVel += (3.0 * (myVolume) * (currentPres - myTargetPres)) *
      halfDeltaT / (Constant::PRESSUREFACTOR * W);
    //  Advance the box volume velocity due to thermostat forces.
    myEpsilonVel += (3.0 * twiceKE / myNumFree) * halfDeltaT / W;
    myEpsilonVel *= exp(-myEtaVolVel * halfDeltaT);

    //  Advance the box volume.
    Real fac = exp(3.0 * myEpsilonVel * 2. * halfDeltaT);
    myVolume *= fac;
    myTopo->rescaleVolume(fac);
    uncache();  

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Barostat -- after force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::PostForceBarostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();
    //  Get the new KE.
    const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);
    const Real MolKE = molecularKineticEnergy(myTopo, myVelocities);

    //  Calculate the current pressure.  Units: (bar)
    Real currentPres = computeMolecularPressure(myEnergies, myVolume, MolKE);

    //  Advance the box volume velocity.
    myEpsilonVel *= exp(-myEtaVolVel * halfDeltaT);
    myEpsilonVel += (3.0 * twiceKE / myNumFree) * halfDeltaT / W;
    //  Advance the box volume velocity.  Divide by PRESSUREFACTOR to get
    //  correct (fs)^-1 units.
    myEpsilonVel += (3.0 * (myVolume) * (currentPres - myTargetPres)) *
      halfDeltaT / (Constant::PRESSUREFACTOR * W);

    //  Advance the volume thermostat variable velocity.
    myEtaVolVel += (W * myEpsilonVel * myEpsilonVel - kbT) * halfDeltaT / Qv;
    //  Advance the volume thermostat variable.
    myEtaV += myEtaVolVel * halfDeltaT;

    //  New volume thermostat kinetic energy.
    Real pEtaV = (Qv * 0.5) * (myEtaVolVel * myEtaVolVel);
    //  New volume thermostat potential energy
    Real VEtaV = myEtaV * kbT;
    //  New box volume kinetic energy.
    Real pVol = (W * 0.5) * (myEpsilonVel * myEpsilonVel);
    //  New box volume potential energy.
    Real VVol = (myVolume * myTargetPres) / Constant::PRESSUREFACTOR;

    //  Add the energy from the extended system barostat.
    (*myEnergies)[ScalarStructure::INTEGRATOR] += pEtaV + VEtaV + pVol + VVol;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Chemostat -- prior to force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::PreForceChemostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;
    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;

    //  Advance the chemostat velocity
    myLambdaVel += (myTargetMu - (*myEnergies).deltaMu()) * halfDeltaT / Qd;

    //  Advance the chemostat velocity due to the thermostat
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
  void oSGIntegrator::PostForceChemostat() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    //  Advance the chemostat velocity
    myLambdaVel *= exp(-myEtaLambdaVel * halfDeltaT);
    
    //  Advance the chemostat velocity
    myLambdaVel += (myTargetMu - (*myEnergies).deltaMu()) * halfDeltaT / Qd;

    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;
    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;

#ifdef DEBUG_CHEMOSTAT
    report.precision(6);
    report << hint << "*******************************************************" << endr;
    report << hint << "Molecule being transformed: " << T+1 << ", current step: " << NumTransSteps << endr;
    report << hint << "Current identity is: " << OldType << ", being inserted? " << Insert << endr;
    report << hint << "TargetMu = " << myTargetMu << endr;
    report << hint << "Lambda = " << myLambda << ", " << myTopo->molecules[T].lambda << endr;
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
  //  modifyForces().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::modifyForces() {
  
    // exit if we have not chosen a transforming molecule yet
    if (T == -1) return;

    // if lambda should be larger than the current stage # or smaller than the previous stage # then
    // use a temporary lambda = newStage or lambda = oldStage so that we don't end up multiplying any force
    // parameters by a number that is larger than 1 or less than 0.
    int oldStage = thisStage - 1;
    Real PreModLambda = myLambda;

    // if myLambda is larger than thisStage then temporarily set myLambda to be just less than thisStage;
    // the original myLambda will be restored at the end of this function so that thisStage will be changed
    // by the function checkForTransformation()
    if (myLambda > thisStage) {
      myLambda = thisStage;
      myTopo->molecules[T].lambda = thisStage - 0.000001;}
    // if myLambda is less than oldStage then temporarily set myLambda to be just more than oldStage;
    // the original myLambda will be restored at the end of this function so that thisStage will be changed
    // by the function checkForTransformation()
    else if (myLambda < oldStage) {
      myLambda = oldStage;
      myTopo->molecules[T].lambda = oldStage + 0.000001;}
    // if myLambda is in between thisStage and oldStage then we don't need to adjust myLambda
    else myTopo->molecules[T].lambda = myLambda;

    unsigned int localAtom = 0;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the atomic charges
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all atoms on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
         iter != myTopo->molecules[T].atoms.end(); iter++) {

      // update the atom's charge (TIM -- check this when you try a molecule w/partial charges!!!)
      if (Insert) {
        int s = oldStage;
        myTopo->atoms[*iter].scaledCharge = transformMaps[OldType].atomTypes[localAtom].insert_charge[s]
          * (myLambda - oldStage) + transformMaps[OldType].atomTypes[localAtom].delete_charge[s] * (thisStage - myLambda);}
      else {
        int ds = myNumStages - 1 - oldStage;
        myTopo->atoms[*iter].scaledCharge = transformMaps[OldType].atomTypes[localAtom].delete_charge[ds]
          * (myLambda - oldStage) + transformMaps[OldType].atomTypes[localAtom].insert_charge[ds] * (thisStage - myLambda);}
      myTopo->atoms[*iter].scaledCharge *= Constant::SQRTCOULOMBCONSTANT;

      // increment the atom counter
      localAtom++;
      
#ifdef DEBUG_MODIFIER_ATOM
        report.precision(8);
        report << hint << "____________________________________________________________________" << endr;
        report << hint << "Looking at atom # " << (*iter) << ", type = "
               << myTopo->atomTypes[myTopo->atoms[*iter].type].name << endr;
        report << hint << "Atom's stage number: " << myTopo->atoms[*iter].stageNumber << ", thisStage = "
               << thisStage << ", Lambda = " << myLambda << endr;
        report << hint << "New q = " << myTopo->atoms[*iter].scaledCharge / Constant::SQRTCOULOMBCONSTANT << endr;
        report << hint << "DeltaQ = " << myTopo->atoms[*iter].deltaQ / Constant::SQRTCOULOMBCONSTANT << endr;
        report << hint << "____________________________________________________________________" << endr;
#endif
    } // end loop over atoms

    // restore the original myLambda
    myLambda = PreModLambda;
  } // end function modifyForces()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  checkForStageCompletion().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::checkForStageCompletion() {
  
    Transformed = false;

    // determine which transformation state the molecule is currently in
    int newStage = static_cast<int>(floor(myLambda)) + 1;

    // determine if the molecule has been completely transformed
    if (myLambda > (static_cast<Real>(myNumStages) - 0.002) ) {
      //  add to the # of molecules of OldType if we were inserting the molecule
      if (Insert) {
        N[OldType]++;
        myTopo->iSGNumMols[OldType]++;
      }
      // subtract from the # of molecules of OldType if we were deleting the molecule
      else {
        N[OldType]--;
        myTopo->iSGNumMols[OldType]--;
      }

      // set Lambda for this molecule to be equal to 0
      myTopo->molecules[T].lambda = 0.0;
      Transformed = true;
      Succeeded = true;
    }
    else if (myLambda < 0.002) {
      // set Lambda for this molecule to be equal to 0
      myTopo->molecules[T].lambda = 0.0;
      Transformed = true;
      Succeeded = false;
    }

    else if (newStage != thisStage) {

      // the transformation attempt has not yet completed,
      // but it has reached a new stage so we need to update
      // the LJ and coloumb force parameters which are used for DeltaMu calculations.

      unsigned int localAtom = 0;
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the stomic charges
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // loop over all atoms on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
           iter != myTopo->molecules[T].atoms.end(); iter++) {

        // to make sure we compute the intramolecular Coulomb force correctly, we must
        // assign the proper charge to Qold and Qnew and update deltaQ
        int s = newStage - 1;
        if (Insert) {
          myTopo->atoms[*iter].Qold = transformMaps[T].atomTypes[localAtom].delete_charge[s];
          myTopo->atoms[*iter].Qnew = transformMaps[T].atomTypes[localAtom].insert_charge[s];
        }
        else {
          int ds = myNumStages - 1 - s;
          myTopo->atoms[*iter].Qold = transformMaps[T].atomTypes[localAtom].insert_charge[ds];
          myTopo->atoms[*iter].Qnew = transformMaps[T].atomTypes[localAtom].delete_charge[ds];
        }

        myTopo->atoms[*iter].deltaQ = myTopo->atoms[*iter].Qnew - myTopo->atoms[*iter].Qold;
        myTopo->atoms[*iter].Qold *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[*iter].Qnew *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[*iter].deltaQ *= Constant::SQRTCOULOMBCONSTANT;

        // update the atom counter
        localAtom++;
      } // end loop over atoms

      // Lastly, update the transformation stage indicator
      thisStage = newStage;
    }

    // if the transformation attempt has completed then compute the
    // average deviation in the conserved quantity
    if (Transformed) {
      // add to the # of transformation attempts since the last trajectory output
      numAttempts++;

      // compute the average of the conserved quantity
      AveCQ /= static_cast<Real>(NumTransSteps);
      AveCQSq /= static_cast<Real>(NumTransSteps);

      // compute the average Lambda temperature over the course of
      // this transformation attempt
      Real AveLambdaT = LambdaT / static_cast<Real>(TotSteps);
      (*myEnergies)[ScalarStructure::LAMBDA_TEMPERATURE] = AveLambdaT;

      // compute the fluctuation in the CQ
      (*myEnergies)[ScalarStructure::CQFLUCTUATION] = AveCQSq - power(AveCQ,2);
      (*myEnergies)[ScalarStructure::DELTATIME] = NumTransSteps;

      // output the system properties to the AllEnergiesFile
      myEnergies->output(true);

      // determine if it is time to output the system trajectory
      if (myTrajectoryFreq == numAttempts) {
        myEnergies->trajectory(true);
        numAttempts = 0;
      }
      // reset the transformation stage indicator to stage 1
      thisStage = 1;
    } // end if (Transformed) statement       
  } // end function checkForStageCompletion()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // run
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::run(int numTimesteps) {
    for (int i = 0; i < numTimesteps; i++) {
      
      // randomly pick a molecule to be inserted/deleted
      if (Transformed) pickNewMolecule();

      // integrate the equations of motion
      preStepModify();
      doHalfKick();
      doDriftOrNextIntegrator();
      modifyForces();
      calculateForces();
      do2ndHalfKick();
      postStepModify();

      //  add to the conserved quantity averages
      Real MyCQ = kineticEnergy(myTopo, myVelocities) + myEnergies->potentialEnergy()
        + (*myEnergies)[ScalarStructure::INTEGRATOR];
      AveCQ += MyCQ;
      AveCQSq += power(MyCQ, 2);
      TotSteps++;

      // check to see if a transformation stage has completed
      checkForStageCompletion();

      // if the transformation is complete, do some bookkeeping
      if (Transformed) setForcesAfterTransformation();

    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  initialize().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::initialize(GenericTopology *topo,
                                 Vector3DBlock *positions,
                                 Vector3DBlock *velocities,
                                 ScalarStructure *energies) {

    STSIntegrator::initialize(topo, positions, velocities, energies);

  } // end function initialize

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Add the modifiers
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::addModifierAfterInitialize() {
    adoptPreStepModifier(new oSGModifierPreForceChemostat(this,2));
    adoptPreStepModifier(new ModifierPreForceThermostat<oSGIntegrator>(this,1));
    adoptPostStepModifier(new ModifierPostForceThermostat<oSGIntegrator>(this,2));
    adoptPostStepModifier(new oSGModifierPostForceChemostat(this,1));
    STSIntegrator::addModifierAfterInitialize();
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  addModifierBeforeInitialize().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::addModifierBeforeInitialize() {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // indexTypes()
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::indexTypes(GenericTopology *topo, const STAGE &stage, const PAR& par, bool dihedralMultPSF, int theSeed) {
  
    // store the random number seed
    seed = theSeed;
  
    // determine the # of osmotic components
    numOsmComp = stage.components.size();   
    
    // make sure we didn't get too many fugacities
    if(getIdNoAlias() == "muVTVerlet") {
      if(numOsmComp > myNumComp)           
        report << error << "[muVTIntegrator::indexTypes] oh no, you tried to pass "   
               << numOsmComp << " fugacities, when you cannot pass any more than "
               << myNumComp  << " fugacities... bye!" << endr;
       }
    else {
      if(numOsmComp >= myNumComp)           
        report << error << "[NfPTIntegrator::indexTypes] oh no, you tried to pass "   
               << numOsmComp << " fugacities, when you cannot pass any more than "
               << myNumComp - 1 << " fugacities... bye!" << endr;
    }      
    
    // copy the information in the STAGE file to the local transformMaps object
    transformMaps.resize(numOsmComp);
    
    for (unsigned int i=0; i<numOsmComp; i++) {
      transformMaps[i].fugacity = stage.components[i].fugacity;
      transformMaps[i].NumberOfStages = stage.components[i].NumberOfStages;
      transformMaps[i].myMolecule = stage.components[i].myMolecule;
      for (unsigned int j=0; j<stage.components[i].atomTypes.size(); j++)
        transformMaps[i].atomTypes.push_back(stage.components[i].atomTypes[j]);
      transformMaps[i].myCoordinates = stage.components[i].myCoordinates;

      // compute the chemical potentials --> kT * ln (fugacity/kT)
      Real fugacity = transformMaps[i].fugacity;
      transformMaps[i].fugacity = kbT * log (fugacity / kbT * fugacityFactor);


      /// incorporate the PSF information into the transformMaps
      //~~~~~~~~~~
      // atoms
      //~~~~~~~~~~
      // loop over all atoms in the PSF object
      for (vector<PSF::Atom>::const_iterator atom = stage.components[i].myStructure.atoms.begin();
           atom != stage.components[i].myStructure.atoms.end(); ++atom) {

        // Two data members for AtomType, name and mass
        AtomType tempatomtype;
        tempatomtype.name = atom->atom_type;
        tempatomtype.mass = atom->mass;
        tempatomtype.symbolName = atomTypeToSymbolName(atom->atom_type);
        tempatomtype.charge = atom->charge;

        Atom tempatom;
       
        // First, we need to find the index. (an integer corresponding
        // to the type of the atom)
        for (unsigned int j=0; j<transformMaps[i].atomTypes.size(); j++) {
          std::string myType_name = transformMaps[i].atomTypes[j].type_name;
          bool Found = false;
          
          // loop over all atomTypes in the topology
          for (unsigned int k=0; k<topo->atomTypes.size(); k++) {
            // is this atom type already in the topology?
            if (myType_name == topo->atomTypes[k].name) {
              Found = true;
              tempatom.type = k;
              break;}
          }

          // if this atomType is not in the topology then we must add it
          if (!Found) {
            topo->atomTypes.push_back(tempatomtype);
            tempatom.type = topo->atomTypes.size() - 1;}           
        } // end loop over j

        // Now, the scaled charge.  This is straightforward.
        tempatom.scaledCharge = (atom->charge)*Constant::SQRTCOULOMBCONSTANT;
        tempatom.scaledMass = atom->mass;
        
        // add the atom to the transformMap
        transformMaps[i].atoms.push_back(tempatom);
      } // end loop over PSF atoms


      /// center the COM of each stock osmotic molecule (places COM at the origin)
      /// to do this, first compute the COM of the stock molecule, in case it is not already at (0,0,0)
      /// here I have borrowed some code from the file topologyutilities.cpp -- TIM
      Real sumM = 0.0;
      Real m = transformMaps[i].atoms[0].scaledMass;
      sumM += m;
      Vector3D COM(transformMaps[i].myCoordinates[0] * m);
      Vector3D tempC(0.0,0.0,0.0);
      for (unsigned int a=1; a<transformMaps[i].atoms.size(); a++) {
        m = transformMaps[i].atoms[i].scaledMass;
        sumM += m;
        Vector3D tempX(transformMaps[i].myCoordinates[a] * m);
        Vector3D tempY(tempX - tempC);
        Vector3D tempT(COM + tempY);
        tempC = (tempT - COM) - tempY;
        COM = tempT;
      }
      transformMaps[i].myMolecule.position = COM / sumM;

      /// now subtract the COM coordinates from each atom on the molecule
      for (unsigned int a=0; a<transformMaps[i].atoms.size(); a++)
        transformMaps[i].myCoordinates[a] -= transformMaps[i].myMolecule.position;
      Vector3D Zip(0.0,0.0,0.0);
      transformMaps[i].myMolecule.position = Zip;

/*
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      /// add bonds to transformMaps
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // First create look-up-table
      map<string,vector<PAR::Bond>::const_iterator> bondLookUpTable;
      for (vector<PAR::Bond>::const_iterator bond = par.bonds.begin(); bond != par.bonds.end(); ++bond)
        bondLookUpTable[bond->atom1+","+bond->atom2] = bond;


      // Find the parameters from PAR
      for (vector<PSF::Bond>::const_iterator bond = stage.components[i].myStructure.bonds.begin();
           bond != stage.components[i].myStructure.bonds.end(); ++bond){

        // store the ID numbers of the bonded atoms
        int atom1 = bond->atom1-1;
        int atom2 = bond->atom2-1;

        // store the type names of the bonded atoms
        string bond1(transformMaps[i].atomTypes[atom1].type_name);
        string bond2(transformMaps[i].atomTypes[atom2].type_name);

        map<string,vector<PAR::Bond>::const_iterator>::const_iterator currentbond = bondLookUpTable.find(bond1+","+bond2);
        if(currentbond == bondLookUpTable.end()){
          currentbond = bondLookUpTable.find(bond2+","+bond1);
        }

        // if we still have not found this bond type in the PAR object, report an error
        if(currentbond == bondLookUpTable.end()){
          report << error << "Could not find bond \'"<<bond1<<"\'-\'"<<bond2<<"\' ("<<bond->atom1 <<","<< bond->atom2<<")"<<std::endl;
          for (map<string,vector<PAR::Bond>::const_iterator>::const_iterator i = bondLookUpTable.begin();
               i != bondLookUpTable.end(); i++){
            report << plain << i->first<<std::endl;
          }
          report << endr;
        }

        // if we have found this bond type then copy the bond parameters
        // into the transformMaps topology
        Bond tempbond;
        tempbond.springConstant = currentbond->second->forceConstant;
        tempbond.restLength = currentbond->second->distance;
        tempbond.atom1 = atom1;
        tempbond.atom2 = atom2;
        transformMaps[i].bonds.push_back(tempbond);

        // populate the vector of bonds maintained at each atom
        //report << std::endl <<"Size of Bonds Vector "<<topo->bonds.size() <<endr;
        transformMaps[i].atoms[atom1].mybonds.push_back((transformMaps[i].bonds.size())-1);
        transformMaps[i].atoms[atom2].mybonds.push_back((transformMaps[i].bonds.size())-1);
        // output to screen for testing purposes
        //for (int j = 0; j < transformMaps[i].atoms[atom1].mybonds.size(); j++){
        //report <<"Atom " << atom1 << " bond index = "<< transformMaps[i].atoms[atom1].mybonds[j] <<endr;
        //}
        //for (int k = 0; k < transformMaps[i].atoms[atom2].mybonds.size(); k++){
        //report <<"Atom " << atom2 << " bond index = "<< transformMaps[i].atoms[atom2].mybonds[k] <<endr;
        //}
      } // end loop over PSF bonds section
*/
/*
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      /// add angles to transformMaps
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // First create look-up-table
      map<string,vector<PAR::Angle>::const_iterator> angleLookUpTable;
      for (vector<PAR::Angle>::const_iterator angle = par.angles.begin(); angle != par.angles.end(); ++angle)
        angleLookUpTable[angle->atom1+","+angle->atom2+","+angle->atom3] = angle;

      // Find the parameters from PAR
      // loop over the angle list in the PSF object
      for (vector<PSF::Angle>::const_iterator angle = stage.components[i].myStructure.angles.begin();
           angle != stage.components[i].myStructure.angles.end(); ++angle){

        // store the ID numbers of the atoms in this angle
        int atom1 = angle->atom1-1;
        int atom2 = angle->atom2-1;
        int atom3 = angle->atom3-1;

        // store the type names of the atoms in this angle
        string angle1(transformMaps[i].atomTypes[atom1].type_name);
        string angle2(transformMaps[i].atomTypes[atom2].type_name);
        string angle3(transformMaps[i].atomTypes[atom3].type_name);

        map<string,vector<PAR::Angle>::const_iterator>::const_iterator currentangle = angleLookUpTable.find(angle1+","+angle2+","+angle3);
        if(currentangle == angleLookUpTable.end())
          currentangle = angleLookUpTable.find(angle3+","+angle2+","+angle1);

        // if we still have not found this angle type in the PAR object, report an error
        if(currentangle == angleLookUpTable.end())
          report << error << "Could not find angle \'"<<angle1<<"\'-\'"<<angle2<<"\'-\'"<<angle3<<"\'."<<endr;

        // if we have found this angle type then copy the angle parameters
        // into the transformMaps topology
        Angle tempangle;
        tempangle.atom1 = atom1;
        tempangle.atom2 = atom2;
        tempangle.atom3 = atom3;
        tempangle.forceConstant = currentangle->second->forceConstant;
        tempangle.restAngle = dtor(currentangle->second->angleval);
        if (currentangle->second->ub_flag){ // do we want defaults for these
          tempangle.ureyBradleyConstant = currentangle->second->k_ub;
          tempangle.ureyBradleyRestLength = currentangle->second->r_ub;
        }
        // no Urey-Bradley term specified
        else{
          tempangle.ureyBradleyConstant = 0.0;
          tempangle.ureyBradleyRestLength = 0.0;
        }

        transformMaps[i].angles.push_back(tempangle);
      } // end loop over PSF angle section
*/
/*
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      /// add dihedrals to transformMaps
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // First create look-up-table
      map<string,vector<PAR::Dihedral>::const_iterator> dihedralLookUpTable;
      for (vector<PAR::Dihedral>::const_iterator dihedral = par.dihedrals.begin();
           dihedral != par.dihedrals.end(); ++dihedral)
        dihedralLookUpTable[dihedral->atom1+","+dihedral->atom2+","+dihedral->atom3+","+dihedral->atom4] = dihedral;

      // Find the parameters from PAR
      // loop over the dihedral list in the PSF object
      for(vector<PSF::Dihedral>::const_iterator dihedral = stage.components[i].myStructure.dihedrals.begin();
          dihedral != stage.components[i].myStructure.dihedrals.end(); ++dihedral) {

        // store the ID numbers of the atoms in this dihedral
        int atom1 = dihedral->atom1 - 1;
        int atom2 = dihedral->atom2 - 1;
        int atom3 = dihedral->atom3 - 1;
        int atom4 = dihedral->atom4 - 1;

        // store the type names of the atoms in this dihedral
        string dihedral1 = transformMaps[i].atomTypes[atom1].type_name;
        string dihedral2 = transformMaps[i].atomTypes[atom2].type_name;
        string dihedral3 = transformMaps[i].atomTypes[atom3].type_name;
        string dihedral4 = transformMaps[i].atomTypes[atom4].type_name;

        map<string,vector<PAR::Dihedral>::const_iterator>::const_iterator currentdihedral =
          dihedralLookUpTable.find(dihedral1+","+dihedral2+","+dihedral3+","+dihedral4);

        // if this dihedral type has not been found, try reversing the order of the atom types
        if(currentdihedral == dihedralLookUpTable.end())
          currentdihedral = dihedralLookUpTable.find(dihedral4+","+dihedral3+","+dihedral2+","+dihedral1);

        // Try wildcards if necessary
        if(currentdihedral == dihedralLookUpTable.end()){
          currentdihedral = dihedralLookUpTable.find("X,"+dihedral2+","+dihedral3+",X");
          if(currentdihedral == dihedralLookUpTable.end())
            currentdihedral = dihedralLookUpTable.find("X,"+dihedral3+","+dihedral2+",X");
        }

        // if we still have not found this dihedral type in the PAR object, report an error
        if(currentdihedral == dihedralLookUpTable.end())
          report << error << "Could not find dihedral \'"<<dihedral1<<"\'-\'"<<dihedral2<<"\'-\'"<<dihedral3<<"\'-\'"<<dihedral4<<"\'."<<endr;

        // if we have found this dihedral type then copy the
        // dihedral parameters into the transformMaps topology
        Torsion torsion;
        torsion.atom1         = atom1;
        torsion.atom2         = atom2;
        torsion.atom3         = atom3;
        torsion.atom4         = atom4;
        torsion.periodicity   = currentdihedral->second->periodicity;
        torsion.forceConstant = currentdihedral->second->forceConstant;
        torsion.phaseShift    = dtor(currentdihedral->second->phaseShift);
        torsion.multiplicity  = currentdihedral->second->multiplicity;
        if(topo->dihedrals.empty() ||
           topo->dihedrals[topo->dihedrals.size()-1].atom1 != atom1 ||
           topo->dihedrals[topo->dihedrals.size()-1].atom2 != atom2 ||
           topo->dihedrals[topo->dihedrals.size()-1].atom3 != atom3 ||
           topo->dihedrals[topo->dihedrals.size()-1].atom4 != atom4){
          if(dihedralMultPSF){
            torsion.periodicity.resize(1);
            torsion.forceConstant.resize(1);
            torsion.phaseShift.resize(1);
            torsion.multiplicity = 1;
          }
          transformMaps[i].dihedrals.push_back(torsion);
        }
        else {
          if(dihedralMultPSF){
            Torsion& tmp = topo->dihedrals[topo->dihedrals.size()-1];
            if(tmp.multiplicity > torsion.multiplicity )
              report << error<< "PSF multiplicity definition of dihedral ("
                     << dihedral1 <<","
                     << dihedral2 <<","
                     << dihedral3 <<","
                     << dihedral4 <<") exceeded PAR definition.";
            tmp.periodicity.push_back(torsion.periodicity[tmp.multiplicity]);
            tmp.forceConstant.push_back(torsion.forceConstant[tmp.multiplicity]);
            tmp.phaseShift.push_back(torsion.phaseShift[tmp.multiplicity]);
            tmp.multiplicity++;
          }
          else {
            report << error << "Unexpected PSF multiplicity definition of dihedral ("
                   << dihedral1 <<","
                   << dihedral2 <<","
                   << dihedral3 <<","
                   << dihedral4 <<") occurred, use dihedral multiplicity PSF.";
          }
        }
      } // end loop over PSF dihedral section
*/
/*
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      /// add impropers to transformMaps
      ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // First create look-up-table
      map<string,vector<PAR::Improper>::const_iterator> improperLookUpTable;
      for (vector<PAR::Improper>::const_iterator improper = par.impropers.begin();
           improper != par.impropers.end(); improper++ )
        improperLookUpTable[improper->atom1+","+improper->atom2+","+improper->atom3+","+improper->atom4] = improper;

      // Find the parameters from PAR
      // loop over the improper list in the PSF object
      for(vector<PSF::Improper>::const_iterator improper = stage.components[i].myStructure.impropers.begin();
          improper != stage.components[i].myStructure.impropers.end(); improper++ ) {

        // store the ID numbers of the atoms in this improper
        int atom1 = improper->atom1 - 1;
        int atom2 = improper->atom2 - 1;
        int atom3 = improper->atom3 - 1;
        int atom4 = improper->atom4 - 1;

        // store the type names of the atoms in this improper
        string improper1 = transformMaps[i].atomTypes[atom1].type_name;
        string improper2 = transformMaps[i].atomTypes[atom2].type_name;
        string improper3 = transformMaps[i].atomTypes[atom3].type_name;
        string improper4 = transformMaps[i].atomTypes[atom4].type_name;

        map<string,vector<PAR::Improper>::const_iterator>::const_iterator currentimproper =
           improperLookUpTable.find(improper1+","+improper2+","+improper3+","+improper4);
        if(currentimproper == improperLookUpTable.end())
          currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+","+improper1);

        // Try wildcards if necessary
        // 2) A - X - X - B
        if(currentimproper == improperLookUpTable.end()){
          currentimproper = improperLookUpTable.find(improper1+",X,X,"+improper4);
          if(currentimproper == improperLookUpTable.end())
            currentimproper = improperLookUpTable.find(improper4+",X,X,"+improper1);
        }

        // 3) X - A - B - C
        if(currentimproper == improperLookUpTable.end()){
          currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+","+improper4);
          if(currentimproper == improperLookUpTable.end())
            currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+",X");
        }

        // 4) X - A - B - X
        if(currentimproper == improperLookUpTable.end()){
          currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+",X");
          if(currentimproper == improperLookUpTable.end())
            currentimproper = improperLookUpTable.find("X,"+improper3+","+improper2+",X");
        }

        // 5) X - X - A - B
        if(currentimproper == improperLookUpTable.end()){
          currentimproper = improperLookUpTable.find("X,X,"+improper3+","+improper4);
          if(currentimproper == improperLookUpTable.end())
            currentimproper = improperLookUpTable.find(improper4+","+improper3+",X,X");
        }

        // if we still have not found this improper type in the PAR object, report an error
        if(currentimproper == improperLookUpTable.end())
          report << error << "Could not find improper."<<endr;

        // if we have found this improper type then copy the
        // improper parameters into the transformMaps topology
        Torsion torsion;
        torsion.atom1         = atom1;
        torsion.atom2         = atom2;
        torsion.atom3         = atom3;
        torsion.atom4         = atom4;
        torsion.periodicity.push_back(currentimproper->second->periodicity);
        torsion.forceConstant.push_back(currentimproper->second->forceConstant);
        torsion.phaseShift.push_back(dtor(currentimproper->second->phaseShift));
        torsion.multiplicity  = 1;
        transformMaps[i].impropers.push_back(torsion);
        //     report << plain<< (wildcard ? "#":"")
        //            << improper1 <<","
        //            << improper2 <<","
        //            << improper3 <<","
        //            << improper4 <<","
        //            << torsion.atom1 <<","
        //            << torsion.atom2 <<","
        //            << torsion.atom3 <<","
        //            << torsion.atom4 <<","
        //            << torsion.periodicity[0] <<","
        //            << torsion.periodicity.size() <<","
        //            << torsion.forceConstant[0] <<","
        //            << torsion.phaseShift[0] <<","
        //            << torsion.multiplicity << endr;
      } // end loop over PSF improper section
*/      
    } // end loop over osmotic components

    ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// put the LennardJonesParameters in the topology
    ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // get some array sizes
    unsigned int sizeAtomTypes  = topo->atomTypes.size();
    unsigned int sizeNonbondeds = par.nonbondeds.size();
    unsigned int sizeNbfixs     = par.nbfixs.size();

    topo->lennardJonesParameters.resize(sizeAtomTypes);

    // loop over all atomtypes
    for(unsigned int i=0;i<sizeAtomTypes;i++){
      int ti = 0;
      unsigned int bi = sizeNonbondeds;

      for(unsigned int k=0;k<sizeNonbondeds;++k){
        int ok = equalWildcard(par.nonbondeds[k].atom,topo->atomTypes[i].name);
              
        if(ok > ti){
          bi = k;
          ti = ok;
        }
        if(ti > 1) break;
      } // end loop over k

      if(ti<=0)
        report << error << "Could not find matching parameter nonbonded of atom \'"
               << topo->atomTypes[i].name << "\'." << endr;
       for(unsigned int j=i;j<sizeAtomTypes;j++){
        int tj = 0;
        unsigned int bj = sizeNonbondeds;
        for(unsigned int k=0;k<sizeNonbondeds;++k){
          int ok = equalWildcard(par.nonbondeds[k].atom,topo->atomTypes[j].name);
          if(ok > tj){
            bj = k;
            tj = ok;
          }
          if(tj > 1) break;
        } // end loop over k

        if(tj<=0)
          report << error << "Could not find matching parameter nonbonded of atom \'"
                 << topo->atomTypes[j].name << "\'." << endr;

        LennardJonesParameters paramsij;

        //Charmm28
        Real sigma_i = par.nonbondeds[bi].sigma;
        Real sigma_j = par.nonbondeds[bj].sigma;
        Real sigma14_i = par.nonbondeds[bi].sigma14;
        Real sigma14_j = par.nonbondeds[bj].sigma14;

        Real epsilon_i = par.nonbondeds[bi].epsilon;
        Real epsilon_j = par.nonbondeds[bj].epsilon;
        Real epsilon14_i = par.nonbondeds[bi].epsilon14;
        Real epsilon14_j = par.nonbondeds[bj].epsilon14;

        Real r_ij = sigma_i + sigma_j;
        Real e_ij = sqrt(epsilon_i * epsilon_j);
        Real r14_ij = sigma14_i + sigma14_j;
        Real e14_ij = sqrt(epsilon14_i * epsilon14_j);

        paramsij.A = power<12>(r_ij) * e_ij;
        paramsij.B = 2 * power<6>(r_ij) * e_ij;
        paramsij.A14 = power<12>(r14_ij) * e14_ij;
        paramsij.B14 = 2 * power<6>(r14_ij) * e14_ij;
        //report << topo->atomTypes[i].name<<","<<bi<<","<<topo->atomTypes[j].name<<","<<bj<<" # "<<paramsij.A<<","<<paramsij.B<<","<<paramsij.A14<<","<<paramsij.B14<<endr;
        topo->lennardJonesParameters.set(i, j, paramsij);
      } // wns loop over atom type j
    } // end loop over atom type i

    ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /// put the LennardJonesParameters (NBfix) in the topology
    ///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(unsigned int k=0;k<sizeNbfixs;++k){

      //report << par.nbfixs[k].atom1 << ","<<par.nbfixs[k].atom2<<endr;

      int ti = 0;
      int tj = 0;
      unsigned int bi = sizeNbfixs;
      unsigned int bj = sizeNbfixs;

      for(unsigned int i=0;i<sizeAtomTypes;i++){
        int ok = equalWildcard(par.nbfixs[k].atom1,topo->atomTypes[i].name);
        if(ok > ti){
          bi = i;
          ti = ok;
        }
        if(ti > 1) break;
      }

      if(ti <=0) // ???
        continue;

      for(unsigned int j=0;j<sizeAtomTypes;j++){
        int ok = equalWildcard(par.nbfixs[k].atom2,topo->atomTypes[j].name);
        if(ok > tj){
          bj = j;
          tj = ok;
        }
        if(tj > 1) break;
      }

      if(tj<=0) report << error << "Could not find matching parameter nbfix of atoms \'"
                       << par.nbfixs[k].atom1 << "\' - '" << par.nbfixs[k].atom2 << "\'." << endr;

      LennardJonesParameters paramsij;

      paramsij.A = par.nbfixs[k].a;
      paramsij.B = par.nbfixs[k].b;
      paramsij.A14 = par.nbfixs[k].a14;
      paramsij.B14 = par.nbfixs[k].b14;
      topo->lennardJonesParameters.set(bi, bj, paramsij);

      //report << topo->atomTypes[bi].name<<","<<topo->atomTypes[bj].name<<","<<paramsij.A<<","<<paramsij.B<<","<<paramsij.A14<<","<<paramsij.B14<<endr;
    } // end loop over NbFix types    
  } // end function indexTypes

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  setForcesAfterTransformation().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::setForcesAfterTransformation() {

    unsigned int localAtom = 0;

    /// For the two cases of:
    /// 1) a molecule insertion succeeded, and
    /// 2) a molecule deletion failed, we simply store the permanent atomic
    /// charges but do not make any other changes to the topology
    if ( (Insert && Succeeded) ||
         (!(Insert) && !(Succeeded)) ) {
      // store the molecule´s permanent atomic charges
      // loop over all atoms on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
	  iter != myTopo->molecules[T].atoms.end(); iter++) {

	// update the charge
        myTopo->atoms[*iter].scaledCharge = myTopo->atomTypes[localAtom].charge;
        myTopo->atoms[*iter].scaledCharge *= Constant::SQRTCOULOMBCONSTANT;

	// update deltaQ, Qold, and Qnew
	myTopo->atoms[*iter].deltaQ = 0.0;
	myTopo->atoms[*iter].Qold = 0.0;
	myTopo->atoms[*iter].Qnew = 0.0;
	localAtom++;
      } // end loop over atoms
        
      // make sure that the molecule's type #s are correct
      myTopo->molecules[T].newtype = myTopo->molecules[T].type;
        
    } // end if statement

    /// For the two cases of:
    /// 1) a molecule deletion succeeded, and
    /// 2) a molecule insertion failed, then we have to completely remove the 
    /// transforming molecule/atoms from the topology
    else {
      // remove the transforming molecule´s information from the topology
      // create a series of iterators to remove the proper array elements
      Vector3DBlock::iterator myRemover; 
      std::vector<Atom>::iterator myAtomRemover; 
      std::vector<Molecule>::iterator myMoleculeRemover = myTopo->molecules.begin() + T;
//      std::vector<Bond>::iterator myBondsRemover;
//      std::vector<Angle>::iterator myAnglesRemover;
//      std::vector<Torsion>::iterator myTorsionsRemover;
//      std::vector<Bond::Constraint>::iterator myConstraintRemover;

/*
      // loop over all bonds on molecule myT
      for (std::vector<int>::iterator iter = topo->molecules[myT].bondList.begin();
	   iter != topo->molecules[myT].bondList.end(); iter++) {
        // remove this bond from the exclusion list
        topo->exclusions.remove(topo->bonds[*iter].atom1, topo->bonds[*iter].atom2, EXCLUSION_FULL);

        // remove this bond from the constraint list
        for (unsigned int i=0; i<topo->bondRattleShakeConstraints.size(); i++) {
          if(topo->bondRattleShakeConstraints[i].atom1 == topo->bonds[*iter].atom1 &&
             topo->bondRattleShakeConstraints[i].atom2 == topo->bonds[*iter].atom2) {
            myConstraintRemover = topo->bondRattleShakeConstraints.begin() + i;
            topo->bondRattleShakeConstraints.erase(myConstraintRemover);
            break;
          } // end if constraint == constraint
        } // end loop over constraints

        // remove this bond from the list
        myBondsRemover = topo->bonds.begin() + (*iter);
	topo->bonds.erase(myBondsRemover);
      } // end loop over bonds  /// need to readjust molecule bondlists!!!
*/
/*
      // loop over all angles on molecule myT
      for (std::vector<int>::iterator iter = topo->molecules[myT].angleList.begin();
	   iter != topo->molecules[myT].angleList.end(); iter++) {
        // remove this angle from the exclusion list
        if( topo->exclude != ExclusionType::ONE2 )
          topo->exclusions.remove(topo->angles[*iter].atom1, topo->angles[*iter].atom3, EXCLUSION_FULL);

        // it is possible that there could be a bond constraint associated with this
        // angle if this is an H-X-H angle. Here we check for and remove any such constraint.
        for (unsigned int i=0; i<topo->bondRattleShakeConstraints.size(); i++) {
          if(topo->bondRattleShakeConstraints[i].atom1 == topo->angles[*iter].atom1 &&
             topo->bondRattleShakeConstraints[i].atom2 == topo->angles[*iter].atom3) {
            myConstraintRemover = topo->bondRattleShakeConstraints.begin() + i;
            topo->bondRattleShakeConstraints.erase(myConstraintRemover);
            break;
          } // end if constraint == constraint
        } // end loop over constraints

        // remove this angle from the list
        myAnglesRemover = topo->angles.begin() + (*iter);
	topo->angles.erase(myAnglesRemover);
      } // end loop over angles  /// need to readjust molecule bondlists!!!
*/
/*
      // loop over all dihedrals on molecule myT
      for (std::vector<int>::iterator iter = topo->molecules[myT].dihedralList.begin();
	   iter != topo->molecules[myT].dihedralList.end(); iter++) {
        if( topo->exclude != ExclusionType::ONE4 )
          topo->exclusions.remove(topo->dihedrals[*iter].atom1, topo->dihedrals[*iter].atom4, EXCLUSION_FULL);
        else
          topo->exclusions.remove(topo->dihedrals[*iter].atom1, topo->dihedrals[*iter].atom4, EXCLUSION_MODIFIED);

        // remove this dihedral from the list
        myTorsionsRemover = topo->dihedrals.begin() + (*iter);
	topo->dihedrals.erase(myTorsionsRemover);
      } // end loop over dihedrals  /// need to readjust molecule bondlists!!!
*/
/*
      // loop over all impropers on molecule myT
      for (std::vector<int>::iterator iter = topo->molecules[myT].improperList.begin();
	   iter != topo->molecules[myT].improperList.end(); iter++) {
        // remove this improper from the list
        myTorsionsRemover = topo->impropers.begin() + (*iter);
	topo->impropers.erase(myTorsionsRemover);
      } // end loop over impropers  /// need to readjust molecule bondlists!!!
*/

      /// compute the # of atoms to be removed from the topology
      unsigned int numAtomsToRemove = myTopo->molecules[T].size();

      // loop over all molecules after and including the transformed molecule
      for (unsigned int m=T; m<myTopo->molecules.size(); m++) {
        //  Loop over the atoms on this molecule
        for (unsigned int a=0; a<myTopo->molecules[m].size(); a++) {

          //  Current atom #
          int atom;

          /// remove the elements belonging to the transformed atoms
          if (m == T) {
            atom = myTopo->molecules[T][a];
              
            // remove the Atom element
            myAtomRemover = myTopo->atoms.begin() + atom;
            myTopo->atoms.erase(myAtomRemover);

            // delete the velocity
            myRemover = myVelocities->begin() + atom;
            myVelocities->erase(myRemover);

            // delete the positions
            myRemover = myPositions->begin() + atom;
            myPositions->erase(myRemover);

            // delete the forces
            myRemover = myForces->begin() + atom;
            myForces->erase(myRemover);
          }  // end if m == T statement

          /// for elements of atoms on molecules other than the transforming molecule, we must update the atom ID#s
          else {
            // list of atom ID#s contained in the molecule object
            myTopo->molecules[m].atoms[a] -= numAtomsToRemove;
            atom = myTopo->molecules[m][a];

            // ID# of this atom
            myTopo->atoms[atom].atomNum = atom;

            // the ID# of the molecule that these atoms belong to
            myTopo->atoms[atom].molecule--;
          } // end else statement
        } // end loop over a
      } // end loop over m 
        
      /// remove this molecule from the list
      myTopo->molecules.erase(myMoleculeRemover); 
          
      /// update the # of degrees of freedom and number of molecules
      Real Qofac = 1.0 / myTopo->degreesOfFreedom;
      Real Wfac = 1.0 / (myTopo->degreesOfFreedom + 3);
           
      myTopo->degreesOfFreedom = 3 * myTopo->atoms.size() - 3 - myTopo->bondRattleShakeConstraints.size();
      myNumFree = myTopo->degreesOfFreedom;
      NumMols = myTopo->molecules.size();

      /// compute new thermostat and barostat masses (helps with stability for large changes in # of atoms)
      Qofac *= myTopo->degreesOfFreedom;
      Wfac *= myTopo->degreesOfFreedom + 3;
      Qo *= Qofac;
      W *= Wfac;
        
      /// mark the cell lists as out of date
      myTopo->uncacheCellList();
        
    } // end else statement  
  } // end function setForcesAfterTransformation()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  pickNewMolecule().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::pickNewMolecule() {
    
    unsigned int localAtom = 0;

    // reset lambda to be close to zero
    myLambda = 0.005;
      
    //  randomly select a molecule type to be inserted/deleted
    int myChoice = static_cast<int>(randomNumber(seed)*numOsmComp) % numOsmComp;
    OldType = transformMaps[myChoice].myMolecule.type;

    // update the # of transformation stages in the integrator
    myNumStages = transformMaps[myChoice].NumberOfStages;

    // randomly determine if we should insert or delete a molecule of this type
    Real my2ndChoice = static_cast<int>(randomNumber(seed) * 2) % 2;
    if (my2ndChoice == 0) Insert = false;   /// delete a molecule
    else Insert = true;                     /// insert a molecule

    //  if the # of molecules of the chosen type is zero, then we cannot hoose to delete this type
    if (N[OldType] == 0) {
      my2ndChoice = 1;
      Insert = true;
    }

       
    /// enter this block if we are to delete a molecule
    if (my2ndChoice == 0) {
      // loop over all molecules and make a list of all available molecules of type Old
      std::vector<unsigned int> availableMols;
      for (unsigned int i=0; i<myTopo->molecules.size(); i++) 
        if (myTopo->molecules[i].type == OldType) availableMols.push_back(i);

      //  choose a random molecule Old type molecule
      int my3rdChoice = (static_cast<int>(randomNumber(seed)*availableMols.size())) % availableMols.size();
        
      // let the integrator know that this is the new transforming molecules
      T = availableMols[my3rdChoice];

      // reassign this molecule's Lambda since it is to be transformed
      myTopo->molecules[T].lambda = myLambda;

      // change this molecule's newtype to be different from its type, which is
      // a flag for the Lennard-Jones force function to do a deletion calculation
      myTopo->molecules[T].newtype = myTopo->molecules[T].type + 1;

      // Update the new target chemical potential:
      // myTargetMu = -kT*ln(f/kT * V/N)
      deletionDeltaMu(transformMaps[myChoice].fugacity,OldType);          
        
      // loop over all atoms on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
                                iter != myTopo->molecules[T].atoms.end(); iter++) {

        // the index into the proper element of insert_charge and delete_charge for deletion
        int ds = myNumStages - 1;

        // update stageNumber, deltaQ, Qold, and Qnew, and alphaLJ
        myTopo->atoms[*iter].stageNumber = transformMaps[myChoice].atomTypes[localAtom].stage;
        myTopo->atoms[*iter].Qold        = transformMaps[myChoice].atomTypes[localAtom].insert_charge[ds];
        myTopo->atoms[*iter].Qnew        = transformMaps[myChoice].atomTypes[localAtom].delete_charge[ds];
        myTopo->atoms[*iter].deltaQ      = myTopo->atoms[*iter].Qnew - myTopo->atoms[*iter].Qold;
        myTopo->atoms[*iter].Qold       *= Constant::SQRTCOULOMBCONSTANT;          
        myTopo->atoms[*iter].Qnew       *= Constant::SQRTCOULOMBCONSTANT;          
        myTopo->atoms[*iter].deltaQ     *= Constant::SQRTCOULOMBCONSTANT;          
        myTopo->atoms[*iter].alphaLJ     = transformMaps[myChoice].atomTypes[localAtom].alphaLJ;

        // update the localAtom counter
        localAtom++;
          
      } // end loop over atoms
    } /// end deletion block
      
    /// enter this block if we are to insert a molecule
    else {
      // copy the stock Old type molecule into the incoming molecule's data structures
      Molecule mol = transformMaps[OldType].myMolecule;
      mol.newtype = mol.type;

      // reassign this molecule's Lambda since it is to be transformed
      T = myTopo->molecules.size();
      myTopo->molecules.push_back(mol);
      myTopo->molecules[T].lambda = myLambda;

      // get the current box length
      Real boxLength = pow(myVolume,1.0/3.0);

      // assign a random xyz position to the center of mass of this molecule
      myTopo->molecules[T].position.x = 0.0;
      myTopo->molecules[T].position.y = 0.0;
      myTopo->molecules[T].position.z = 0.0;
      myTopo->molecules[T].position.x = randomNumber(seed) * boxLength - boxLength * 0.5;
      myTopo->molecules[T].position.y = randomNumber(seed) * boxLength - boxLength * 0.5;
      myTopo->molecules[T].position.z = randomNumber(seed) * boxLength - boxLength * 0.5;

      // assign a random orientation to this molecule
      Real Phi = 2 * 3.1415926 * randomNumber(seed);
      Real Theta = 2 * 3.1415926 * randomNumber(seed);
      Real Psi = 2 * 3.1415926 * randomNumber(seed);

      // compute the coordinate tranformation matrix (see Allen & Tildesley)
      Array<Real,2> TransferMatrix;
      TransferMatrix.resize(ArraySizes(3)(3));
      TransferMatrix[0][0] = cos(Phi)*cos(Psi) - sin(Phi)*cos(Theta)*sin(Psi);
      TransferMatrix[0][1] = sin(Phi)*cos(Psi) + cos(Phi)*cos(Theta)*sin(Psi);
      TransferMatrix[0][2] = sin(Theta)*sin(Psi);
      TransferMatrix[1][0] = -cos(Phi)*sin(Psi) - sin(Phi)*cos(Theta)*cos(Psi);
      TransferMatrix[1][1] = -sin(Phi)*sin(Psi) + cos(Phi)*cos(Theta)*cos(Psi);
      TransferMatrix[1][2] = sin(Theta)*cos(Psi);
      TransferMatrix[2][0] = sin(Phi)*sin(Theta);
      TransferMatrix[2][1] = -cos(Phi)*sin(Theta);
      TransferMatrix[2][2] = cos(Theta);

      // loop over all the atoms on this molecule and add them to the topology
      // one by one with an xyz coordinate consistent with the COM coordinates that were
      // just chosen above
      for (unsigned int i=0; i<transformMaps[OldType].myCoordinates.size(); i++) {
        // update the atom id numbers in the molecule object
        myTopo->molecules[T].atoms[i] = myTopo->atoms.size();
          
        // copy in the atom object
        myTopo->atoms.push_back(transformMaps[OldType].atoms[i]);
          
        // set the proper attributes for this new atom
        unsigned int myAtom = myTopo->molecules[T].atoms[i];
        myTopo->atoms[myAtom].molecule    = T;
        myTopo->atoms[myAtom].stageNumber = transformMaps[myChoice].atomTypes[i].stage;
        myTopo->atoms[myAtom].Qold        = transformMaps[myChoice].atomTypes[i].delete_charge[0];
        myTopo->atoms[myAtom].Qnew        = transformMaps[myChoice].atomTypes[i].insert_charge[0];
        myTopo->atoms[myAtom].deltaQ      = myTopo->atoms[myAtom].Qnew - myTopo->atoms[myAtom].Qold;
        myTopo->atoms[myAtom].Qold       *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[myAtom].Qnew       *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[myAtom].deltaQ     *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[myAtom].alphaLJ     = transformMaps[myChoice].atomTypes[i].alphaLJ;
        myTopo->atoms[myAtom].name        = transformMaps[myChoice].atomTypes[localAtom].type_name;
        myTopo->atoms[myAtom].atomNum     = myAtom; 
        myTopo->atoms[myAtom].hvyAtom     = 1;       ///< The size of the Group for Heavy Atom Come First ordering
                                                     ///< for hydrogen hvyAtom=0, else, hvyAtom = 1+#of attached hydrogens

 
        // assign this atom an initial force (must be zero net force to be consistent with a true
        // lambda = 0 state)
        Vector3D NewXYZ(0.0,0.0,0.0);
        myForces->push_back(NewXYZ);

        // apply the new molecule orientation to this atom
        NewXYZ.x = transformMaps[OldType].myCoordinates[i].x * TransferMatrix[0][0]
                 + transformMaps[OldType].myCoordinates[i].y * TransferMatrix[0][1]
                 + transformMaps[OldType].myCoordinates[i].z * TransferMatrix[0][2];
        NewXYZ.y = transformMaps[OldType].myCoordinates[i].x * TransferMatrix[1][0]
                 + transformMaps[OldType].myCoordinates[i].y * TransferMatrix[1][1]
                 + transformMaps[OldType].myCoordinates[i].z * TransferMatrix[1][2];
        NewXYZ.z = transformMaps[OldType].myCoordinates[i].x * TransferMatrix[2][0]
                 + transformMaps[OldType].myCoordinates[i].y * TransferMatrix[2][1]
                 + transformMaps[OldType].myCoordinates[i].z * TransferMatrix[2][2];

        // assign this atom its new xyz coordinate
        myPositions->push_back(myTopo->molecules[T].position + NewXYZ);

        // assign this atom an xyz velocity consistent with the target temperature
        // (i.e. select from a Maxwell-Boltzmann distribution)
        Real kbToverM = sqrt(kbT/transformMaps[OldType].atoms[i].scaledMass);
        NewXYZ.x = kbToverM * randomGaussian(6.0, seed);
        NewXYZ.y = kbToverM * randomGaussian(6.0, seed);
        NewXYZ.z = kbToverM * randomGaussian(6.0, seed);
        myVelocities->push_back(NewXYZ);
      }

      /// must add to the forces list (bonds, angles, constraints, etc.)
      myTopo->molecules[T].bondList.clear();
      myTopo->molecules[T].angleList.clear();
      myTopo->molecules[T].dihedralList.clear();
      myTopo->molecules[T].improperList.clear();
/*
      // loop over all bonds on molecule myT and add them to the topology
      for (unsigned int i=0; i<transformMaps[Old].bonds.size(); i++) {

        // create a temporary bond
        Bond tempBond(transformMaps[Old].bonds[i]);

        // update the ID#s of the two atoms in this bond
        tempBond.atom1 = topo->molecules[myT].atoms[ transformMaps[Old].bonds[i].atom1 ];
        tempBond.atom2 = topo->molecules[myT].atoms[ transformMaps[Old].bonds[i].atom2 ];

        // add this bond to the topology
        topo->bonds.push_back(tempBond);

        // store this bond index number in the molecule's bond list
        topo->molecules[myT].bondList.push_back(topo->bonds.size() - 1);

        // update the vector of bonds maintained at each atom
        topo->atoms[tempBond.atom1].mybonds.push_back((topo->bonds.size())-1);
        topo->atoms[tempBond.atom2].mybonds.push_back((topo->bonds.size())-1);

        // add this bond to the exclusion list
        topo->exclusions.add(tempBond.atom1, tempBond.atom2, EXCLUSION_FULL );

        // add this bond to the constraint list
        topo->bondRattleShakeConstraints.push_back(Bond::Constraint(tempBond.atom1,
                                                                    tempBond.atom2,
                                                                    tempBond.restLength));

      } // end loop over bonds
*/
/*
      // loop over all angles on molecule myT and add them to the topology
      for (unsigned int i=0; i<transformMaps[Old].angles.size(); i++) {

        // create a temporary angle
        Angle tempAngle(transformMaps[Old].angles[i]);

        // update the ID#s of the three atoms in this angle
        tempAngle.atom1 = topo->molecules[myT].atoms[ transformMaps[Old].angles[i].atom1 ];
        tempAngle.atom2 = topo->molecules[myT].atoms[ transformMaps[Old].angles[i].atom2 ];
        tempAngle.atom3 = topo->molecules[myT].atoms[ transformMaps[Old].angles[i].atom3 ];

        // add this angle to the topology
        topo->angles.push_back(tempAngle);

        // store this angle index number in the molecule's angle list
        topo->molecules[myT].angleList.push_back(topo->angles.size() - 1);

        // add this angle to the exclusion list
        if( topo->exclude != ExclusionType::ONE2 )
          topo->exclusions.add(tempAngle.atom1, tempAngle.atom3, EXCLUSION_FULL );

        // it is possible that there could be a bond constraint associated with this
        // angle if this is an H-X-H angle. Here we check for and add any such constraint.
        int a1 = tempAngle.atom1;
        int a2 = tempAngle.atom2;
        int a3 = tempAngle.atom3;

        PairIntSorted p1(a1,a2);
        PairIntSorted p2(a2,a3);

        Real M1 = topo->atoms[a1].scaledMass;
        Real M3 = topo->atoms[a3].scaledMass;
        // For regular water M3 = M1 = 1.008, and for heavy waters,
        // M1 or M3 should be 2 ~ 3 times heavier
        if((M1<5) && (M3 <5)){

	  // ... properly retrieve the right length!
	  int b1 = -1;
	  int b2 = -1;
	  const std::vector< int >& mybonds1 = topo->atoms[a1].mybonds;
	  const std::vector< int >& mybonds3 = topo->atoms[a3].mybonds;

	  for (unsigned int j = 0; j < mybonds1.size(); j++){
	    PairIntSorted s1(topo->bonds[mybonds1[j]].atom1,topo->bonds[mybonds1[j]].atom2);
	    if(p1 == s1){
	      for (unsigned int k = 0; k < mybonds3.size(); k++){
	        PairIntSorted s2(topo->bonds[mybonds3[k]].atom1,topo->bonds[mybonds3[k]].atom2);
	        if(p2 == s2){
		  if(b1 > -1 || b2 > -1)
		    report << warning << "Found already two bonds ("<<b1<<","<<b2<<") for angle "<<i<<"."<<endr;
		  else {
		    b1 = mybonds1[j];
		    b2 = mybonds3[k];
		  }
	        }
	      }
	    }
	  }
	  if(b1 > -1 && b2 > -1){
	    Real b     = topo->bonds[b1].restLength;
	    Real c     = topo->bonds[b2].restLength;
	    Real alpha = topo->angles[i].restAngle;
	    topo->bondRattleShakeConstraints.push_back(Bond::Constraint(a1,a3,sqrt(b*b+c*c-2*b*c*cos(alpha))));
	  }
	  else
	    report << debug(10) << "Could not find two matching bonds for angle "<<i<<"."<<endr;
        }
      } // end loop over angles
*/
/*
      // loop over all dihedrals on molecule myT and add them to the topology
      for (unsigned int i=0; i<transformMaps[Old].dihedrals.size(); i++) {

        // create a temporary dihedral
        Torsion tempTorsion(transformMaps[Old].dihedrals[i]);

        // update the ID#s of the four atoms in this dihedral
        tempTorsion.atom1 = topo->molecules[myT].atoms[ transformMaps[Old].dihedrals[i].atom1 ];
        tempTorsion.atom2 = topo->molecules[myT].atoms[ transformMaps[Old].dihedrals[i].atom2 ];
        tempTorsion.atom3 = topo->molecules[myT].atoms[ transformMaps[Old].dihedrals[i].atom3 ];
        tempTorsion.atom4 = topo->molecules[myT].atoms[ transformMaps[Old].dihedrals[i].atom4 ];

        // add this dihedral to the topology
        topo->dihedrals.push_back(tempTorsion);

        // store this dihedral index number in the molecule's dihedral list
        topo->molecules[myT].dihedralList.push_back(topo->dihedrals.size() - 1);

        // add this dihedral to the exclusion list
        if( topo->exclude != ExclusionType::ONE2 &&
            topo->exclude != ExclusionType::ONE3 ) {
          if ( topo->exclude == ExclusionType::ONE4 )
            topo->exclusions.add(tempTorsion.atom1, tempTorsion.atom4, EXCLUSION_FULL );
          else topo->exclusions.add(tempTorsion.atom1, tempTorsion.atom4, EXCLUSION_MODIFIED );
        }
      } // end loop over dihedrals
*/
/*
      // loop over all impropers on molecule myT and add to the topology
     for (unsigned int i=0; i<transformMaps[Old].impropers.size(); i++) {

        // create a temporary improper
        Torsion tempTorsion(transformMaps[Old].impropers[i]);

        // update the ID#s of the four atoms in this improper
        tempTorsion.atom1 = topo->molecules[myT].atoms[ transformMaps[Old].impropers[i].atom1 ];
        tempTorsion.atom2 = topo->molecules[myT].atoms[ transformMaps[Old].impropers[i].atom2 ];
        tempTorsion.atom3 = topo->molecules[myT].atoms[ transformMaps[Old].impropers[i].atom3 ];
        tempTorsion.atom4 = topo->molecules[myT].atoms[ transformMaps[Old].impropers[i].atom4 ];

        // add this improper to the topology
        topo->impropers.push_back(tempTorsion);

        // store this improper index number in the molecule's improper list
        topo->molecules[myT].improperList.push_back(topo->impropers.size() - 1);

      } // end loop over impropers
*/
      /// update the # of degrees of freedom and number of molecules
      Real Qofac = 1.0 / myTopo->degreesOfFreedom;
      Real Wfac = 1.0 / (myTopo->degreesOfFreedom + 3);
           
      myTopo->degreesOfFreedom = 3 * myTopo->atoms.size() - 3 - myTopo->bondRattleShakeConstraints.size();
      myNumFree = myTopo->degreesOfFreedom;
      NumMols = myTopo->molecules.size();

      /// compute new thermostat and barostat masses (helps with stability for large changes in # of atoms)
      Qofac *= myTopo->degreesOfFreedom;
      Wfac *= myTopo->degreesOfFreedom + 3;
      Qo *= Qofac;
      W *= Wfac;
        
      /// mark the cell lists as out of date
      myTopo->uncacheCellList();
        
      /// Update the new target chemical potential:
      /// myTargetMu = kT*ln(f/kT * V/(N+1))
      insertionDeltaMu(transformMaps[myChoice].fugacity,OldType);
    } /// end insertion block


    // new Lambda velocity chosen from MB distribution
    // width of Maxwell-Boltzmann distribution (in units of [1/fs])
    Real sigma = sqrt(kbT / Qd);

    // assign a new positive lambda velocity
    Real width = randomGaussianNumber(0.0,sigma,seed);
    width *= width;
    myLambdaVel = sqrt(width);

    //  reset the conserved quantity averages
    AveCQ = 0.0;
    AveCQSq = 0.0;
    NumTransSteps = 0;

    // reset the transformation flag
    Transformed = false;

    // reset AllEnergiesFile and trajectoryFile flags
    myEnergies->output(false);
    myEnergies->trajectory(false);
    (*myEnergies).clear();

    /// unCache any precomputed lists or values for all of the forces
    uncache();

  } // end function pickNewMolecule()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // readXSCs (eXstended System Coordinates)
  // read this only if an XSC file is specified in the config file.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void oSGIntegrator::readXSCs(const string myFile, GenericTopology *topo) {

    // create the XSC reader object
    XSCReader xscReader;

    if(!xscReader.open(myFile))
      report << error << "Can't open XSC file \'"<< myFile <<"\'."<<endr;

    // create the storage object for the XSC information
    XSC xsc;

    // read in the XSC file, which will contain the integrator's
    // time-zero extended system variable values
    if(!(xscReader >> xsc))
      report << error << "Could not parse XSC file \'"<< myFile <<"\'."<<endr;
    report << plain << "Using XSC file \'"<< myFile <<"\'."<< endr;

    // store the time-zero values of the extended system variables
    int Mol = xsc.myMolecule-1;
    unsigned int Osmotic = xsc.old_type;
    T = Mol;
    myLambda = xsc.Lambda;
    myLambdaVel = xsc.Lambda_vel;
    OldType = Osmotic;
    unsigned int NewType = xsc.new_type;
    myEta = xsc.Eta;
    myEtaV = xsc.EtaVol;
    myEtaLambda = xsc.EtaLambda;
    myEtaVel = xsc.Eta_vel;
    myEtaVolVel = xsc.EtaVol_vel;
    myEtaLambdaVel = xsc.EtaLambda_vel;
    myVolume = xsc.Vol;
    myEpsilonVel = xsc.Epsilon_vel;

    Transformed = false;
    if (NewType == OldType) Insert = true;
    else Insert = false;
    report << "Initial values of the extended system variables:" << endr;
    report << "Transforming molecule: " << Mol+1 << endr;
    report << "OldType = " << OldType << endr;
    report << "Insertion? " << Insert << endr;
    report.precision(10);
    report << "Lambda = " << myLambda << endr;
    report << "Lambda velocity (1/fs) = " << myLambdaVel << endr;
    report << "Eta = " << myEta << endr;
    report << "Eta velocity (1/fs) = " << myEtaVel << endr;
    report << "EtaVol = " << myEtaV << endr;
    report << "EtaVol velocity (1/fs) = " << myEtaVolVel << endr;
    report << "EtaLambda = " << myEtaLambda << endr;
    report << "EtaLambda velocity (1/fs) = " << myEtaLambdaVel << endr;
    report << "Volume = " << myVolume << endr;
    report << "Epsilon velocity (1/fs) = " << myEpsilonVel << endr;

    // set the # of transformation stages in the integrator
    myNumStages = transformMaps[Osmotic].NumberOfStages;

    // set the target chemical potential in the integrator
    if (Insert) 
      myTargetMu = transformMaps[Osmotic].fugacity;
    else
      myTargetMu = -transformMaps[Osmotic].fugacity;

    // determine which transformation state the molecule is currently in
    thisStage = static_cast<int>(floor(xsc.Lambda)) + 1;
    report << "Current transform stage = " << thisStage << endr;

    // loop over all atoms on molecule T
    unsigned int localAtom = 0;
    for (std::vector<int>::iterator iter = topo->molecules[Mol].atoms.begin();
	 iter != topo->molecules[Mol].atoms.end(); iter++) {
      
      // update stageNumber, deltaQ, Qold, and Qnew, and alphaLJ
      topo->atoms[*iter].stageNumber = transformMaps[Osmotic].atomTypes[localAtom].stage;

      int s = thisStage - 1;
      if (Insert) {

        topo->atoms[*iter].Qold = transformMaps[Osmotic].atomTypes[localAtom].delete_charge[s];
        topo->atoms[*iter].Qnew = transformMaps[Osmotic].atomTypes[localAtom].insert_charge[s];
      }
      else {
        int ds = myNumStages - 1 - s;
        topo->atoms[*iter].Qold = transformMaps[Osmotic].atomTypes[localAtom].insert_charge[ds];
        topo->atoms[*iter].Qnew = transformMaps[Osmotic].atomTypes[localAtom].delete_charge[ds];
      }

      topo->atoms[*iter].deltaQ = topo->atoms[*iter].Qnew - topo->atoms[*iter].Qold;
      topo->atoms[*iter].Qold *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].Qnew *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].deltaQ *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].alphaLJ = transformMaps[Osmotic].atomTypes[localAtom].alphaLJ;
    } // end loop over atoms
  } // end function readXSCs()
  
}
