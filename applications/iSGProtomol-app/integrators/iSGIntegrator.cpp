//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for iSG dynamics                    //
//                                                                           //
//  Unless modified, this integrator uses the molecular virial to control    //
//  the pressure. -- TIM                                                     //
//  -----------------------------------------------------------------------  //

#include "iSGIntegrator.h"
#include "ScalarStructure.h" 
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h" 
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierNPTShake.h"
#include "ModifierNPTRattle.h"
#include "iSGModifyForces.h"
#include "XSCReader.h"

//#define DEBUG_MODIFIER_BOND
//#define DEBUG_MODIFIER_ANGLE
//#define DEBUG_MODIFIER_TORSION
//#define DEBUG_MODIFIER_ATOM
//#define DEBUG_CHEMOSTAT

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Keyword.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const string iSGIntegrator::keyword("iSGVerlet");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Default or empty constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  iSGIntegrator::iSGIntegrator(): STSIntegrator(),
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
				                          NewType(0),
				                          Transformed(true),
				                          myTargetDeltaMu(0.0),
                                  mylnMassRatio(0.0),
                                  myDMuIG(0.0),
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
                                  TotSteps(0) {}
 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  iSGIntegrator::iSGIntegrator(Real timestep, 
                               unsigned int numComp,
			                         Real temperature,
                               Real pressure, 
                               const vector<Real> &fugacityFrac,
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
      NewType(0), 
      Transformed(true),
      myTargetDeltaMu(0.0),
      mylnMassRatio(0.0),
      myDMuIG(0.0),
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
      TotSteps(0) {

    if(numComp != fugacityFrac.size())
      report << error << "[iSGIntegrator::iSGIntegrator] oh no, you tried to pass fugacityFrac("<<fugacityFrac.size()<<") with a different size than number of components("<<numComp<<")... bye!"<<endr;

    // store the fugacity fractions --> kT * ln (fugacityFrac)
    myTargetMu.resize(myNumComp,0.0);
    myFugacityFrac = fugacityFrac;

    Real sum = 0.0;

    for (unsigned int i=0; i<myNumComp; i++) {
      myTargetMu[i] = myTargetTemp * Constant::BOLTZMANN * log(fugacityFrac[i]);
      sum += fugacityFrac[i];
    }

    // Stop if the fugacity fractions do not sum to 1
    if (!equal(sum,1.0,1e-13))
      report << error << "Fugacity fractions do not sum (=1 + "<<sum-1.0<<") to 1."<< endr;   

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Thermostat -- prior to force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  void iSGIntegrator::PreForceThermostat() {

    //report << hint << "PreForceThermostat" << endr;

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
  void iSGIntegrator::PostForceThermostat() {

    //report << hint << "PostForceThermostat" << endr;

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
  void iSGIntegrator::PreForceBarostat() {

    //report << hint << "PreForceBarostat" << endr;

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
  void iSGIntegrator::PostForceBarostat() {

    //report << hint << "PostForceBarostat" << endr;

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
  void iSGIntegrator::PreForceChemostat() {

    //report << hint << "PreForceChemostat" << endr;

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();   


    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;
    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;


    //  Advance the chemostat velocity
    Real deltaM_vsq = 0.0;
    for (unsigned int i=0; i < myTopo->molecules[T].size(); i++) {   
      //  Current atom # and mass
      int atom = myTopo->molecules[T][i];
      deltaM_vsq += myTopo->atoms[atom].deltaM * ((*myVelocities)[atom]).normSquared();
    }

    //  Advance the chemostat velocity due to the fugacity fraction
    myLambdaVel += (myTargetDeltaMu + myDMuIG - (*myEnergies).deltaMu()
      + (deltaM_vsq * 0.5) + mylnMassRatio) * halfDeltaT / Qd;

    //  Advance the chemostat velocity due to the thermostat
    myLambdaVel *= exp(-myEtaLambdaVel * halfDeltaT);

    //  Advance the chemostat
    myLambda += myLambdaVel * 2. * halfDeltaT;

    //  If Lambda should become negative, reset it to be zero,
    //  and if Lambda should become larger than myNumStages, reset it to be myNumStages.
    if (myLambda < 0) myLambda = 0.0;
    else if (myLambda > myNumStages) myLambda = static_cast<Real>(myNumStages);

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Chemostat -- after force calculations
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  void iSGIntegrator::PostForceChemostat() {

    //report << hint << "PostForceChemostat" << endr;

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();    

    
    //  Advance the chemostat velocity
    myLambdaVel *= exp(-myEtaLambdaVel * halfDeltaT);

    Real deltaM_vsq = 0.0;
    for (unsigned int i=0; i < myTopo->molecules[T].size(); i++) {  
      //  Current atom # and mass
      int atom = myTopo->molecules[T][i];
      deltaM_vsq += myTopo->atoms[atom].deltaM * ((*myVelocities)[atom]).normSquared();
    }

    //  Advance the chemostat velocity due to the fugacity fraction
    myLambdaVel += (myTargetDeltaMu + myDMuIG - (*myEnergies).deltaMu()
      + (deltaM_vsq * 0.5) + mylnMassRatio ) * halfDeltaT / Qd;

    //  Advance the lambda thermostat variable velocity.
    myEtaLambdaVel += (Qd * myLambdaVel * myLambdaVel - Constant::BOLTZMANN * myMuTemp) * halfDeltaT / Ql;
    //  Advance the lambda thermostat variable.
    myEtaLambda += myEtaLambdaVel * halfDeltaT;

#ifdef DEBUG_CHEMOSTAT
    report.precision(6);
    report << hint << "*******************************************************" << endr;
    report << hint << "Molecule being transformed: " << T+1 << ", current step: " << NumTransSteps << endr;
    report << hint << "Current identity is: " << OldType << ", New identity is: " << NewType << endr;
    report << hint << "deltaM_vsq = " << deltaM_vsq << endr;
    report << hint << "TargetDeltaMu = " << myTargetDeltaMu << ", mass ratio term = "
                   << mylnMassRatio << ", DeltaMuIG = " << myDMuIG << endr;
    report << hint << "Lambda = " << myLambda << endr;
    report << hint << "Drift Current DeltaMu = " << (*myEnergies).deltaMu() << endr;
    report << hint << "Bond DeltaMu = " << (*myEnergies)[ScalarStructure::BOND_DELTAMU] << endr;
    report << hint << "Angle DeltaMu = " << (*myEnergies)[ScalarStructure::ANGLE_DELTAMU] << endr;
    report << hint << "Dihedral DeltaMu = " << (*myEnergies)[ScalarStructure::DIHEDRAL_DELTAMU] << endr;
    report << hint << "Improper DeltaMu = " << (*myEnergies)[ScalarStructure::IMPROPER_DELTAMU] << endr;
    report << hint << "LennardJones DeltaMu = " << (*myEnergies)[ScalarStructure::LENNARDJONES_DELTAMU] << endr;
    report << hint << "Coulomb DeltaMu = " << (*myEnergies)[ScalarStructure::COULOMB_DELTAMU] << endr;
    report << hint << "LambdaVel = " << myLambdaVel << endr;
    report << hint << "*******************************************************" << endr;
#endif


    //  New chemostat kinetic energy
    Real pLambda = (Qd * 0.5) * (myLambdaVel * myLambdaVel);
    //  New chemostat potential energy
    Real VLambda = -myLambda * (myTargetDeltaMu)
     + (thisStage - 1 - myLambda) * (myDMuIG + mylnMassRatio);
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
  void iSGIntegrator::doHalfKick() {

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
  void iSGIntegrator::doDrift() {

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
  void iSGIntegrator::do2ndHalfKick() {

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
  // run
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::run(int numTimesteps) {
    for (int i = 0; i < numTimesteps; i++) {

      // randomly pick a molecule to be transformed
      if (Transformed) pickNewMolecule();
      
      // integrate the equations of motion
      PreForceThermostat();
      PreForceBarostat();
      PreForceChemostat();
      doHalfKick();
      doDriftOrNextIntegrator();
      calculateForces();
      do2ndHalfKick();
      PostForceChemostat();
      PostForceBarostat();
      PostForceThermostat();

      //  add to the conserved quantity averages
      Real MyCQ = kineticEnergy(myTopo, myVelocities) + myEnergies->potentialEnergy()
        + (*myEnergies)[ScalarStructure::INTEGRATOR];
      AveCQ += MyCQ;
      AveCQSq += power(MyCQ, 2);
      TotSteps++;

      // check to see if a transformation stage has completed
      checkForCompletion();

      // if the transformation is complete, do some bookkeeping
      if (Transformed) updateTopology();
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  initialize().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::initialize(GenericTopology *topo,
                                 Vector3DBlock *positions,
                                 Vector3DBlock *velocities, 
                                 ScalarStructure *energies) {

    STSIntegrator::initialize(topo, positions, velocities, energies);

    //  store the total number of atoms and molecules
    NumAtoms = myTopo->atoms.size();
    NumMols = myTopo->molecules.size();
    myNumFree = myTopo->degreesOfFreedom;

    //  determine the # of molecules of each component in the mixture
    N.resize(myNumComp);
    myTopo->iSGNumMols.resize(myNumComp);
    for (unsigned int C=0; C < myNumComp; C++) {
      N[C] = 0;
      myTopo->iSGNumMols[C] = 0;
    }
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
      myTopo->molecules[T].newtype = NewType;

      // Compute the residual and ideal gas chemical potential differences:
      // ideal gas chemical potential difference (using the precomputed chemical potentials)
      myDMuIG = myDeltaMuIG[OldType][NewType][thisStage-1];

      // residual chemical potential difference
      // DeltaMu = kT*ln(f1/f0 * N0/(N1+1))
      myTargetDeltaMu = (myTargetMu[NewType] - myTargetMu[OldType])
	      + kbT * log (static_cast<Real>(N[OldType]) / static_cast<Real>(N[NewType]+1));
      myTargetDeltaMu /= static_cast<Real>(myNumStages);
    }

    //  Compute the fixed barostat mass and thermostat masses.
    Qo = myNumFree * kbT * (myTauT * myTauT);
    Qv = kbT * (myTauV * myTauV);
    W = (myNumFree + 3) * kbT * (myTauP * myTauP);
    Qd = (Constant::BOLTZMANN * myMuTemp) * (myTauD * myTauD);
    Ql = (Constant::BOLTZMANN * myMuTemp) * (myTauL * myTauL);

    initializeForces();
    myEnergies->output(false);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NPT RATTLE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  Modifier* iSGIntegrator::createRattleModifier(Real eps, int maxIter){
    return (new ModifierNPTRattle<iSGIntegrator>(eps, maxIter, this));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NPT SHAKE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifier* iSGIntegrator::createShakeModifier(Real eps, int maxIter){
    return (new ModifierNPTShake<iSGIntegrator>(eps, maxIter, this));
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Add the modifiers
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::addModifierAfterInitialize() {
    //adoptPreStepModifier(new ModifierPreForceThermostat<iSGIntegrator>(this,1));   // first thermostat update
    //adoptPreStepModifier(new ModifierPreForceBarostat<iSGIntegrator>(this,2));     // first barostat update
    //adoptPreStepModifier(new iSGModifierPreForceChemostat(this,3));                // first chemostat update

    adoptPreForceModifier(new iSGModifyForces(this,1));                            // update the topology with the new lambda value

    //adoptPostStepModifier(new iSGModifierPostForceChemostat(this,1));              // second chemostat update
    //adoptPostStepModifier(new ModifierPostForceBarostat<iSGIntegrator>(this,2));   // second barostat update
    //adoptPostStepModifier(new ModifierPostForceThermostat<iSGIntegrator>(this,3)); // second thermostat update

    STSIntegrator::addModifierAfterInitialize();
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  addModifierBeforeInitialize().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::addModifierBeforeInitialize() {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // getParameters().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::getParameters(vector< Parameter >& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("components", Value(myNumComp)));
    parameters.push_back(Parameter("temperature", Value(myTargetTemp)));
    parameters.push_back(Parameter("pressure", Value(myTargetPres)));
    parameters.push_back(Parameter("fugacityfrac", Value(myFugacityFrac)));
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
  STSIntegrator* iSGIntegrator::doMake(string&, const vector<Value>& values,ForceGroup* fg)const{
    return new iSGIntegrator(values[0], values[1], values[2], values[3], values[4],
			     values[5], values[6], values[7], values[8], values[9], values[10], fg);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // indexTopology()
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::indexTopology(GenericTopology *topo, const TRANS &trans, const PSF& psf,
                                    const iSGPAR& par, bool dihedralMultPSF, int theSeed) {

    // store the random number seed
    seed = theSeed;

    /// first we need to properly size the topology storage elements
    // create an atomtype map
    map< string, int > tempAtomList;

    tempAtomList[psf.atoms[0].atom_type] = 0;
    for(unsigned int i=1;i<psf.atoms.size();i++){
      
      if( tempAtomList.find( psf.atoms[i].atom_type ) == tempAtomList.end() )
        tempAtomList[ psf.atoms[i].atom_type ] = 0;

    }

    //~~~~~~~~~~~~~~~~~~~~~~~
    // Bonds
    //~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds in the par object
    for( unsigned int i = 0; i < par.bonds.size(); i++ ) {

      // determine if this PAR bond is in our system
      if( ( tempAtomList.find( par.bonds[i].atom1 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.bonds[i].atom2 ) != tempAtomList.end() ) ) {
        
        // it is in our system, so store the bond parameters here
        iSGPAR::Bond tempBond = par.bonds[i];
        bonds.push_back(tempBond);
      }
    } // end loop over i

    // now loop over all bonds in the topology
    for(unsigned int i=0; i<topo->bonds.size(); i++) {

      // the ID numbers of the bonded atoms
      int I = topo->bonds[i].atom1;
      int J = topo->bonds[i].atom2;

      // the atom type of the atoms in this bond
      string atom1(topo->atomTypes[topo->atoms[I].type].name);
      string atom2(topo->atomTypes[topo->atoms[J].type].name);

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      unsigned int j;

      // loop over all bonds in the local topology storage element
      for (j=0; j<bonds.size(); j++) {
  
        if (bonds[j].atom1 == atom1 || bonds[j].atom2 == atom1) {found1 = true;}
        if (bonds[j].atom1 == atom2 || bonds[j].atom2 == atom2) {found2 = true;}
  
        if (found1 && found2) {break;}
      } // end loop over j

      // record the local bond index into the topology
      topo->bonds[i].iSGmodifierIndex = j;
      
    } // end loop over i (all bonds in the topology)

    //~~~~~~~~~~~~~~~~~~~~~~~
    // Angles
    //~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles in the par object
    for( unsigned int i = 0; i < par.angles.size(); i++ ) {

      // determine if this PAR angle is in our system
      if( ( tempAtomList.find( par.angles[i].atom1 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.angles[i].atom2 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.angles[i].atom3 ) != tempAtomList.end() ) ) {

        // it is in our system, so store the angle parameters here
        iSGPAR::Angle tempAngle = par.angles[i];
        tempAngle.angleval = dtor(par.angles[i].angleval);
        angles.push_back(tempAngle); 
      }   
    } // end loop over i

    // now loop over all angles in the topology
    for(unsigned int i=0; i<topo->angles.size(); i++) {
     
      // the ID numbers of the atoms in the angle
      int I = topo->angles[i].atom1;
      int J = topo->angles[i].atom2;
      int K = topo->angles[i].atom3;

      // the atom type of the atoms in this angle
      string atom1 = topo->atomTypes[topo->atoms[I].type].name;
      string atom2 = topo->atomTypes[topo->atoms[J].type].name;
      string atom3 = topo->atomTypes[topo->atoms[K].type].name;

      // loop counter
      unsigned int j;

      // loop over all angles in the modifier
      for (j=0; j<angles.size(); j++) {

        // loop sentinels
        bool found1 = false;
        bool found2 = false;
        bool found3 = false;

        if (angles[j].atom1 == atom1 || angles[j].atom2 == atom1 || angles[j].atom3 == atom1) {found1 = true;}
        if (angles[j].atom1 == atom2 || angles[j].atom2 == atom2 || angles[j].atom3 == atom2) {found2 = true;}
        if (angles[j].atom1 == atom3 || angles[j].atom2 == atom3 || angles[j].atom3 == atom3) {found3 = true;}

        if (found1 && found2 && found3) {break;}
      } // end loop over j

      // record the local angle index into the topology
      topo->angles[i].iSGmodifierIndex = j;
      
    } // end loop over i (all angles in the topology)

    //~~~~~~~~~~~~~~~~~~~~~~~
    // Dihedrals
    //~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all dihedrals in the par object
    for( unsigned int i = 0; i < par.dihedrals.size(); i++ ) {

      // determine if this PAR dihedral is in our system
      if( ( tempAtomList.find( par.dihedrals[i].atom1 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.dihedrals[i].atom2 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.dihedrals[i].atom3 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.dihedrals[i].atom4 ) != tempAtomList.end() ) ) {

        // it is in our system, so store the angle parameters here
        iSGPAR::Dihedral tempDihedral = par.dihedrals[i];
        for (unsigned int j=0; j < par.dihedrals[i].multiplicity.size(); j++) {
          for (int k=0; k < par.dihedrals[i].multiplicity[j]; k++)
            tempDihedral.phaseShift[j][k] = dtor(par.dihedrals[i].phaseShift[j][k]);
        }
        dihedrals.push_back(tempDihedral); 
      } // end if statement   
    } // end loop over i

    // now loop over all dihedrals in the topology
    for(unsigned int i=0; i<topo->dihedrals.size(); i++) {
      
      // the ID numbers of the atoms in the dihedral
      int I = topo->dihedrals[i].atom1;
      int J = topo->dihedrals[i].atom2;
      int K = topo->dihedrals[i].atom3;
      int L = topo->dihedrals[i].atom4;

      // the atom type of the atoms in this dihedral
      string atom1 = topo->atomTypes[topo->atoms[I].type].name;
      string atom2 = topo->atomTypes[topo->atoms[J].type].name;
      string atom3 = topo->atomTypes[topo->atoms[K].type].name;
      string atom4 = topo->atomTypes[topo->atoms[L].type].name;

      // loop counter
      unsigned int j;

      // loop over all dihedrals in the modifier
      for (j=0; j<dihedrals.size(); j++) {

        // loop sentinels
        bool found1 = false;
        bool found2 = false;
        bool found3 = false;
        bool found4 = false;
  
        if (dihedrals[j].atom1 == atom1 || dihedrals[j].atom2 == atom1 ||
            dihedrals[j].atom3 == atom1 || dihedrals[j].atom4 == atom1) {found1 = true;}
        if (dihedrals[j].atom1 == atom2 || dihedrals[j].atom2 == atom2 ||
            dihedrals[j].atom3 == atom2 || dihedrals[j].atom4 == atom2) {found2 = true;}
        if (dihedrals[j].atom1 == atom3 || dihedrals[j].atom2 == atom3 ||
            dihedrals[j].atom3 == atom3 || dihedrals[j].atom4 == atom3) {found3 = true;}
        if (dihedrals[j].atom1 == atom4 || dihedrals[j].atom2 == atom4 ||
            dihedrals[j].atom3 == atom4 || dihedrals[j].atom4 == atom4) {found4 = true;}

        if (found1 && found2 && found3 && found4) {break;}
      } // end loop over j

      // record the local dihedral index into the topology
      topo->dihedrals[i].iSGmodifierIndex = j;

    } // end loop over dihedrals

    //~~~~~~~~~~~~~~~~~~~~~~~
    // Impropers
    //~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers in the par object
    for( unsigned int i = 0; i < par.impropers.size(); i++ ) {

      // determine if this PAR improper is in our system
      if( ( tempAtomList.find( par.impropers[i].atom1 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.impropers[i].atom2 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.impropers[i].atom3 ) != tempAtomList.end() ) &&
          ( tempAtomList.find( par.impropers[i].atom4 ) != tempAtomList.end() ) ) {

        // it is in our system, so store the angle parameters here
        iSGPAR::Improper tempImproper = par.impropers[i];
        for (unsigned int k=0; k < par.impropers[i].phaseShift.size(); k++)
          tempImproper.phaseShift[k] = dtor(par.impropers[i].phaseShift[k]);        
        impropers.push_back(tempImproper); 
      }   
    } // end loop over i

    // nowloop over all impropers in the topology
    for(unsigned int i=0; i<topo->impropers.size(); i++) {
      
      // the ID numbers of the atoms in the improper
      int I = topo->impropers[i].atom1;
      int J = topo->impropers[i].atom2;
      int K = topo->impropers[i].atom3;
      int L = topo->impropers[i].atom4;

      // the atom type of the atoms in this improper
      string atom1 = topo->atomTypes[topo->atoms[I].type].name;
      string atom2 = topo->atomTypes[topo->atoms[J].type].name;
      string atom3 = topo->atomTypes[topo->atoms[K].type].name;
      string atom4 = topo->atomTypes[topo->atoms[L].type].name;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      unsigned int j;

      // loop over all impropers in the modifier
      for (j=0; j<impropers.size(); j++) {
  
        if (impropers[j].atom1 == atom1 || impropers[j].atom2 == atom1 ||
            impropers[j].atom3 == atom1 || impropers[j].atom4 == atom1) {found1 = true;}
        if (impropers[j].atom1 == atom2 || impropers[j].atom2 == atom2 ||
            impropers[j].atom3 == atom2 || impropers[j].atom4 == atom2) {found2 = true;}
        if (impropers[j].atom1 == atom3 || impropers[j].atom2 == atom3 ||
            impropers[j].atom3 == atom3 || impropers[j].atom4 == atom3) {found3 = true;}
        if (impropers[j].atom1 == atom4 || impropers[j].atom2 == atom4 ||
            impropers[j].atom3 == atom4 || impropers[j].atom4 == atom4) {found4 = true;}

        if (found1 && found2 && found3 && found4) {break;}
      } // end loop over j

      // record the local improper index into the topology
      topo->impropers[i].iSGmodifierIndex = j;

    } // end loop over dihedrals

    //~~~~~~~~~~~~~~~~~~~~~~~
    // AtomTypes
    //~~~~~~~~~~~~~~~~~~~~~~~
    // set the # of transformation stages in the integrator
    myNumStages = trans.NumberOfStages;

    // loop over all atom types in the topology
    for(unsigned int i=0; i<topo->atomTypes.size(); i++) {

      // copy this atom type into the local atomType structure
      atomTypes.push_back(iSGPAR::AtomType(myNumComp, topo->atomTypes[i].name, topo->atomTypes[i].symbolName));

      // copy this atomtype's set of charges for each transformation stage into the ChargeMap structure
      chargeMaps.push_back(ChargeMap(myNumComp, myNumStages));
    }

    // loop over all atom types in the local atomType structure
    for (unsigned int i=0; i<atomTypes.size(); i++) {
      
      // get the name of this atom type
      string myName = atomTypes[i].name;
      
      // loop over all atom types in the mqf file
      for (unsigned int j=0; j<trans.atomTypes.size(); j++) {

        // if this TRANS atomType is the same as our atom type (myName)
        // then store the masses, charges, and alphaLJs for this atom type
        if (myName == trans.atomTypes[j].type_name) {
          atomTypes[i].mass = trans.atomTypes[j].mass;
          atomTypes[i].charge = trans.atomTypes[j].charge;
          atomTypes[i].stageNumber = trans.atomTypes[j].stage;
          chargeMaps[i].old_charge = trans.atomTypes[j].old_charge;
          chargeMaps[i].new_charge = trans.atomTypes[j].new_charge;
          chargeMaps[i].alphaLJ = trans.atomTypes[j].alphaLJ;
        }
      } // end loop over j
    } // end loop over i

    // fill in the ideal gas state chemical potential difference matrix
    myDeltaMuIG.resize(ArraySizes(myNumComp)(myNumComp)(myNumStages));
    for (unsigned int o=0; o<myNumComp; o++)
      for (unsigned int n=0; n<myNumComp; n++)
        for (unsigned int s=0; s<myNumStages; s++)
          myDeltaMuIG[o][n][s] = trans.DeltaMuIG[o][n][s];

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // readXSCs (eXstended System Coordinates)
  // read this only if an XSC file is specified in the config file.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::readXSCs(const string myFile, GenericTopology *topo) {

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
    T = xsc.myMolecule-1;
    myLambda = xsc.Lambda;
    myLambdaVel = xsc.Lambda_vel;
    OldType = xsc.old_type;
    NewType = xsc.new_type;
    myEta = xsc.Eta;
    myEtaV = xsc.EtaVol;
    myEtaLambda = xsc.EtaLambda;
    myEtaVel = xsc.Eta_vel;
    myEtaVolVel = xsc.EtaVol_vel;
    myEtaLambdaVel = xsc.EtaLambda_vel;
    myVolume = xsc.Vol;
    myEpsilonVel = xsc.Epsilon_vel;

    Transformed = false;

    report << "Initial values of the extended system variables:" << endr;
    report << "Transforming molecule: "                   << T+1 << endr;
    report << "OldType = "                            << OldType << endr;
    report << "NewType = "                            << NewType << endr;
    report.precision(10);
    report << "Lambda = "                            << myLambda << endr;
    report << "Lambda velocity (1/fs) = "         << myLambdaVel << endr;
    report << "Eta = "                                  << myEta << endr;
    report << "Eta velocity (1/fs) = "               << myEtaVel << endr;
    report << "EtaVol = "                              << myEtaV << endr;
    report << "EtaVol velocity (1/fs) = "         << myEtaVolVel << endr;
    report << "EtaLambda = "                      << myEtaLambda << endr;
    report << "EtaLambda velocity (1/fs) = "   << myEtaLambdaVel << endr;
    report << "Volume = "                            << myVolume << endr;
    report << "Epsilon velocity (1/fs) = "       << myEpsilonVel << endr;

    // compute the ideal gas chemical potential
    // difference from changing the molecule's mass
    mylnMassRatio = 0.0;

    // determine which transformation state the molecule is currently in
    thisStage = static_cast<int>(floor(xsc.Lambda)) + 1;

    // loop over all atoms on molecule T
    for (std::vector<int>::iterator iter = topo->molecules[T].atoms.begin();
         iter != topo->molecules[T].atoms.end(); iter++) {

      // the index # into the atomTypes array
      int I = topo->atoms[*iter].type;

      // the atom type name of this atom
      string myName(topo->atomTypes[I].name);

      // find the proper iSGAtomType element
      bool Found = false;
      unsigned int myISGAtomType = 0;
      while (!Found) {
        if (myName == atomTypes[myISGAtomType].name) Found = true;
        else myISGAtomType++;
      }

      // update stageNumber, deltaM, deltaQ, Qold, and Qnew, and alphaLJ
      topo->atoms[*iter].stageNumber = atomTypes[myISGAtomType].stageNumber[OldType];
      topo->atoms[*iter].Qold = chargeMaps[myISGAtomType].old_charge[OldType][NewType][thisStage-1];
      topo->atoms[*iter].Qnew = chargeMaps[myISGAtomType].new_charge[OldType][NewType][thisStage-1];
      topo->atoms[*iter].deltaQ = chargeMaps[myISGAtomType].new_charge[OldType][NewType][thisStage-1]
        - chargeMaps[myISGAtomType].old_charge[OldType][NewType][thisStage-1];
      topo->atoms[*iter].Qold *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].Qnew *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].deltaQ *= Constant::SQRTCOULOMBCONSTANT;
      topo->atoms[*iter].alphaLJ = chargeMaps[myISGAtomType].alphaLJ[OldType][NewType];
      // add the following terms only for the current transformation stage
      if (topo->atoms[*iter].stageNumber == thisStage) {
        topo->atoms[*iter].deltaM = atomTypes[myISGAtomType].mass[NewType]
          - atomTypes[myISGAtomType].mass[OldType];
        
        // update mylnMassRatio for the current stage
        mylnMassRatio += (3 * kbT * 0.5) * ( log(atomTypes[myISGAtomType].mass[OldType])
          - log(atomTypes[myISGAtomType].mass[NewType]) );
      } // end if statement
    } // end loop over atoms

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the BOND parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds on molecule T
    for (std::vector<int>::iterator iter = topo->molecules[T].bondList.begin();
         iter != topo->molecules[T].bondList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = topo->atoms[topo->bonds[*iter].atom1].stageNumber;
      int atom2_stage = topo->atoms[topo->bonds[*iter].atom2].stageNumber;

      // determine which stage this bond should be transformed in
      // (it should be the LARGER of the two stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage);

      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {
        // find the index # of the proper iSGBond element
        int myISGBond = topo->bonds[*iter].iSGmodifierIndex;
        
        // update DeltaK and DeltaR0
        topo->bonds[*iter].DeltaK = bonds[myISGBond].forceConstant[NewType]
          - bonds[myISGBond].forceConstant[OldType];
        topo->bonds[*iter].DeltaR0 = bonds[myISGBond].distance[NewType]
          - bonds[myISGBond].distance[OldType];
      }
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the ANGLE parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles on molecule T
    for (std::vector<int>::iterator iter = topo->molecules[T].angleList.begin();
         iter != topo->molecules[T].angleList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = topo->atoms[topo->angles[*iter].atom1].stageNumber;
      int atom2_stage = topo->atoms[topo->angles[*iter].atom2].stageNumber;
      int atom3_stage = topo->atoms[topo->angles[*iter].atom3].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage);
       
      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {
        // find the index # of the proper iSGAngle element
        int myISGAngle = topo->angles[*iter].iSGmodifierIndex;
        
        // update DeltaK and DeltaTheta0
        topo->angles[*iter].DeltaK = angles[myISGAngle].forceConstant[NewType]
          - angles[myISGAngle].forceConstant[OldType];
        topo->angles[*iter].DeltaTheta0 = angles[myISGAngle].angleval[NewType]
          - angles[myISGAngle].angleval[OldType];
        topo->angles[*iter].Delta_ubK = angles[myISGAngle].k_ub[NewType]
          - angles[myISGAngle].k_ub[OldType];
        topo->angles[*iter].Delta_ubR0 = angles[myISGAngle].r_ub[NewType]
          - angles[myISGAngle].r_ub[OldType];
      }
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the DIHEDRAL parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all dihedrals on molecule T
    for (std::vector<int>::iterator iter = topo->molecules[T].dihedralList.begin();
         iter != topo->molecules[T].dihedralList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = topo->atoms[topo->dihedrals[*iter].atom1].stageNumber;
      int atom2_stage = topo->atoms[topo->dihedrals[*iter].atom2].stageNumber;
      int atom3_stage = topo->atoms[topo->dihedrals[*iter].atom3].stageNumber;
      int atom4_stage = topo->atoms[topo->dihedrals[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);
       
      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {
        // find the index # of the proper iSGDihedral element
        int myISGDihedral = topo->dihedrals[*iter].iSGmodifierIndex;
        
        for (int i=0; i<topo->dihedrals[*iter].multiplicity; i++) {            
          // update DeltaK and DeltaPhase
          topo->dihedrals[*iter].DeltaK[i] = dihedrals[myISGDihedral].forceConstant[NewType][i]
            - dihedrals[myISGDihedral].forceConstant[OldType][i];
          topo->dihedrals[*iter].DeltaPhase[i] = dihedrals[myISGDihedral].phaseShift[NewType][i]
            - dihedrals[myISGDihedral].phaseShift[OldType][i];
        } // end loop over i
      } // end if statement
    } // end loop over dihedrals
        
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the IMPROPER parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers on molecule T
    for (std::vector<int>::iterator iter = topo->molecules[T].improperList.begin();
         iter != topo->molecules[T].improperList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = topo->atoms[topo->impropers[*iter].atom1].stageNumber;
      int atom2_stage = topo->atoms[topo->impropers[*iter].atom2].stageNumber;
      int atom3_stage = topo->atoms[topo->impropers[*iter].atom3].stageNumber;
      int atom4_stage = topo->atoms[topo->impropers[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);

      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {
        // find the index # of the proper iSGImproper element
        int myISGImproper = topo->impropers[*iter].iSGmodifierIndex;
        
        // update DeltaK and DeltaPhase
        topo->impropers[*iter].DeltaK[0] = impropers[myISGImproper].forceConstant[NewType]
          - impropers[myISGImproper].forceConstant[OldType];
        topo->impropers[*iter].DeltaPhase[0] = impropers[myISGImproper].phaseShift[NewType]
          - impropers[myISGImproper].phaseShift[OldType];
      }
    } // end loop over impropers
  } // end function readXSCs()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  modifyForces().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::modifyForces() {

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

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the atomic charges and masses
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // the new mass of the molecule
    Real newMolMass = 0.0;

    // loop over all atoms on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
         iter != myTopo->molecules[T].atoms.end(); iter++) {

      // this is the correct stage for this atom
      // the index # into the atomTypes array
      int I = myTopo->atoms[*iter].type;

      // the atom type name of this atom
      string myName(myTopo->atomTypes[I].name);
      
      // find the proper iSGAtomType element
      bool Found = false;
      unsigned int myISGAtomType = 0;
      while (!Found) {
        if (myName == atomTypes[myISGAtomType].name) Found = true;
        else myISGAtomType++;
      }

      // The mass only needs to be updated if we are at the proper transformation stage for this atom
      if (myTopo->atoms[*iter].stageNumber == thisStage)        
        myTopo->atoms[*iter].scaledMass = atomTypes[myISGAtomType].mass[NewType] * (myLambda - oldStage)
          + atomTypes[myISGAtomType].mass[OldType] * (thisStage - myLambda);

      // update the atom's charge and deltaQ
      int o = OldType;
      int n = NewType;
      int s = oldStage;
      myTopo->atoms[*iter].scaledCharge = chargeMaps[myISGAtomType].new_charge[o][n][s] * (myLambda - oldStage)
        + chargeMaps[myISGAtomType].old_charge[o][n][s] * (thisStage - myLambda);
      myTopo->atoms[*iter].scaledCharge *= Constant::SQRTCOULOMBCONSTANT;

      // update the molecular mass
      newMolMass += myTopo->atoms[*iter].scaledMass;
      
#ifdef DEBUG_MODIFIER_ATOM
        report.precision(8);
        report << hint << "____________________________________________________________________" << endr;
        report << hint << "Looking at atom # " << (*iter) << ", type = "
               << myTopo->atomTypes[myTopo->atoms[*iter].type].name << endr;
        report << hint << "Atom's stage number: " << myTopo->atoms[*iter].stageNumber << ", thisStage = "
               << thisStage << ", Lambda = " << myLambda << endr;
        report << hint << "New m and q = " << myTopo->atoms[*iter].scaledMass << ", "
               << myTopo->atoms[*iter].scaledCharge / Constant::SQRTCOULOMBCONSTANT << endr;
        report << hint << "DeltaM = " << myTopo->atoms[*iter].deltaM << ", DeltaQ = "
               << myTopo->atoms[*iter].deltaQ / Constant::SQRTCOULOMBCONSTANT << endr;
        report << hint << "____________________________________________________________________" << endr;
#endif
    } // end loop over atoms

    // update the mass of molecule T
    myTopo->molecules[T].mass = newMolMass;
        
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the BOND parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].bondList.begin();
      iter != myTopo->molecules[T].bondList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->bonds[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->bonds[*iter].atom2].stageNumber;

      // determine which stage this bond should be transformed in
      // (it should be the LARGER of the two stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage);

      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {      
        // find the index # of the proper iSGBond element
        int myISGBond = myTopo->bonds[*iter].iSGmodifierIndex;
      
        // update the force constant, k, and rest distance, r0
        myTopo->bonds[*iter].springConstant = bonds[myISGBond].forceConstant[NewType] * (myLambda - oldStage)
          + bonds[myISGBond].forceConstant[OldType] * (thisStage - myLambda);
        myTopo->bonds[*iter].restLength = bonds[myISGBond].distance[NewType] * (myLambda - oldStage)
          + bonds[myISGBond].distance[OldType] * (thisStage - myLambda);

        // if we are using rigid bonds then we must loop over the constraint list to find
        // the constraint corresponding to this bond
        if (myTopo->bondRattleShakeConstraints.size() != 0) {

          // get the atom ID#s of the atoms in this bond
          int myAtom1 = myTopo->bonds[*iter].atom1;
          int myAtom2 = myTopo->bonds[*iter].atom2;

          // loop over the constraint list
          for (unsigned int shakeIter = 0; shakeIter < myTopo->bondRattleShakeConstraints.size(); shakeIter++) {

            // get the atom ID#s of the atoms in this constraint
            int shakeAtom1 = myTopo->bondRattleShakeConstraints[shakeIter].atom1;
            int shakeAtom2 = myTopo->bondRattleShakeConstraints[shakeIter].atom2;
            
            // if this is the correct constraint the update the restLength
            if (myAtom1 == shakeAtom1 && myAtom2 == shakeAtom2)
              myTopo->bondRattleShakeConstraints[shakeIter].restLength = myTopo->bonds[*iter].restLength;
          } // end loop over constraints
        } // end if statement

#ifdef DEBUG_MODIFIER_BOND
        report.precision(8);
        report << hint << "Looking at bond #" << (*iter) << endr;
        report << hint << "New k and r0 = " << myTopo->bonds[*iter].springConstant
               << ", " << myTopo->bonds[*iter].restLength << endr;
#endif
      } // end if statement
    } // end loop over bonds

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the ANGLE parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].angleList.begin();
      iter != myTopo->molecules[T].angleList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->angles[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->angles[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->angles[*iter].atom3].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage);
       
      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {          
        // find the index # of the proper iSGAngle element
        int myISGAngle = myTopo->angles[*iter].iSGmodifierIndex;

        // update the force constant, k, and the rest angle, theta0
        myTopo->angles[*iter].forceConstant = angles[myISGAngle].forceConstant[NewType] * (myLambda - oldStage)
          + angles[myISGAngle].forceConstant[OldType] * (thisStage - myLambda);
        myTopo->angles[*iter].ureyBradleyConstant = angles[myISGAngle].k_ub[NewType] * (myLambda - oldStage)
          + angles[myISGAngle].k_ub[OldType] * (thisStage - myLambda);
        myTopo->angles[*iter].restAngle = angles[myISGAngle].angleval[NewType] * (myLambda - oldStage)
          + angles[myISGAngle].angleval[OldType] * (thisStage - myLambda);
        myTopo->angles[*iter].ureyBradleyRestLength = angles[myISGAngle].r_ub[NewType] * (myLambda - oldStage)
          + angles[myISGAngle].r_ub[OldType] * (thisStage - myLambda);

#ifdef DEBUG_MODIFIER_ANGLE
        report.precision(8);
        report << hint << "Looking at angle # " << (*iter) << endr;
        report << hint << "Index # " << myISGAngle << endr;
        report << hint << "New k and theta0 = " << myTopo->angles[*iter].forceConstant << ", "
               << myTopo->angles[*iter].restAngle << endr;
        report << hint << "New k and theta0 for UB = " << myTopo->angles[*iter].ureyBradleyConstant 
               << ", " << myTopo->angles[*iter].ureyBradleyRestLength << endr;
#endif
      } // end if statement
    } // end loop over angles
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the DIHEDRAL parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all dihedrals on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].dihedralList.begin();
       iter != myTopo->molecules[T].dihedralList.end(); iter++) {
      
      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom3].stageNumber;
      int atom4_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);
       
      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {      
        // find the index # of the proper iSGDihedral element
        int myISGDihedral = myTopo->dihedrals[*iter].iSGmodifierIndex;
        
        for (int i=0; i<myTopo->dihedrals[*iter].multiplicity; i++) {
          // update the force constants
          myTopo->dihedrals[*iter].forceConstant[i] =
            dihedrals[myISGDihedral].forceConstant[NewType][i] * (myLambda - oldStage)
            + dihedrals[myISGDihedral].forceConstant[OldType][i] * (thisStage - myLambda);
 
          // update the phase shift angles
          myTopo->dihedrals[*iter].phaseShift[i] =
            dihedrals[myISGDihedral].phaseShift[NewType][i] * (myLambda - oldStage)
            + dihedrals[myISGDihedral].phaseShift[OldType][i] * (thisStage - myLambda);

#ifdef DEBUG_MODIFIER_TORSION
          report.precision(8);
          report << hint << "Looking at dihedral # " << (*iter) << endr;
          report << hint << "New k and P.S. = " << myTopo->dihedrals[*iter].forceConstant[i]
                 << ", " << myTopo->dihedrals[*iter].phaseShift[i] << endr;
          report << hint << "DeltaK = " << myTopo->dihedrals[*iter].DeltaK[i] << endr;
          report << hint << "DeltaPhase = " << myTopo->dihedrals[*iter].DeltaPhase[i] << endr;
#endif
        } // end loop over i
      } // end if statement
    } // end loop over dihedrals

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the IMPROPER parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers on molecule T
    for (vector<int>::iterator iter = myTopo->molecules[T].improperList.begin();
      iter != myTopo->molecules[T].improperList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->impropers[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->impropers[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->impropers[*iter].atom3].stageNumber;
      int atom4_stage = myTopo->atoms[myTopo->impropers[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);

      // update these parameters only if myStage = thisStage
      if (myStage == thisStage) {      
        // find the index # of the proper iSGImproper element
        int myISGImproper = myTopo->impropers[*iter].iSGmodifierIndex;
      
        // update the force constant 
        myTopo->impropers[*iter].forceConstant[0] =
          impropers[myISGImproper].forceConstant[NewType] * (myLambda - oldStage)
          + impropers[myISGImproper].forceConstant[OldType] * (thisStage - myLambda);

        // update the phase shift angle
        myTopo->impropers[*iter].phaseShift[0] =
          impropers[myISGImproper].phaseShift[NewType] * (myLambda - oldStage)
          + impropers[myISGImproper].phaseShift[OldType] * (thisStage - myLambda);

#ifdef DEBUG_MODIFIER_TORSION
        report.precision(8);
        report << hint << "Looking at improper # " << (*iter) << endr;
        report << hint << "New k and P.S. = " << myTopo->impropers[*iter].forceConstant[0]
               << ", " << myTopo->impropers[*iter].phaseShift[0] << endr;
#endif
      } // end if statement
    } // end loop over impropers

    // restore the original myLambda
    myLambda = PreModLambda;
  } // end function modifyForces()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  PickNewMolecule().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::pickNewMolecule() {

    //report << hint << "iSGPickNewMolecule" << endr;

    //  randomly select a molecule type to be transformed
    bool Chosen = false;
    while (!Chosen) {

      OldType = (static_cast<int>(randomNumber(seed)*myNumComp)) % myNumComp;

      //  if the # of molecules of the chosen type is zero, then we cannot
      //  choose to transform this type
      if (N[OldType] != 0) Chosen = true;
    }

    //  choose a random species OldType molecule.
    Chosen = false;
    while (!Chosen) {
      
      // random number
      T = (static_cast<int>(randomNumber(seed)*NumMols)) % NumMols;

      // exit if this molecule type is the same as OldType
      if (myTopo->molecules[T].type == OldType) Chosen = true;
    } // end while loop


    // reassign this molecule's Lambda since it is to be transformed
    myLambda = 0.005;
    myTopo->molecules[T].lambda = myLambda;
      
    // chose a random species for this molecule to be transformed into
    // this species must not the be same species as myChoice
    Chosen = false;
    while (!Chosen) {
      NewType = (static_cast<int>(randomNumber(seed)*myNumComp)) % myNumComp;
      if (NewType != OldType) {Chosen = true;}
    }

    // Update the type and newtype of this molecule
    myTopo->molecules[T].type = OldType;
    myTopo->molecules[T].newtype = NewType;
          
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the atomic mass and charge
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mylnMassRatio = 0.0;
      
    // loop over all atoms on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
         iter != myTopo->molecules[T].atoms.end(); iter++) {
        
      // the index # into the atomTypes array
      int I = myTopo->atoms[*iter].type;
        
      // the atom type name of this atom
      string myName(myTopo->atomTypes[I].name);
        
      // find the proper iSGAtomType element
      bool Found = false;
      unsigned int myISGAtomType = 0;
      while (!Found) {
        if (myName == atomTypes[myISGAtomType].name) Found = true;
        else myISGAtomType++;
      }
        
      // update stageNumber, deltaM, deltaQ, Qold, and Qnew, and alphaLJ
      myTopo->atoms[*iter].stageNumber = atomTypes[myISGAtomType].stageNumber[OldType];
      myTopo->atoms[*iter].Qold = chargeMaps[myISGAtomType].old_charge[OldType][NewType][0];
      myTopo->atoms[*iter].Qnew = chargeMaps[myISGAtomType].new_charge[OldType][NewType][0];
      myTopo->atoms[*iter].deltaQ = chargeMaps[myISGAtomType].new_charge[OldType][NewType][0]
        - chargeMaps[myISGAtomType].old_charge[OldType][NewType][0];
      myTopo->atoms[*iter].Qold *= Constant::SQRTCOULOMBCONSTANT;
      myTopo->atoms[*iter].Qnew *= Constant::SQRTCOULOMBCONSTANT;
      myTopo->atoms[*iter].deltaQ *= Constant::SQRTCOULOMBCONSTANT;
      myTopo->atoms[*iter].alphaLJ = chargeMaps[myISGAtomType].alphaLJ[OldType][NewType];
        
      // update the following terms only for stage 1 atoms
      if (myTopo->atoms[*iter].stageNumber == 1) {
        myTopo->atoms[*iter].deltaM = atomTypes[myISGAtomType].mass[NewType]
          - atomTypes[myISGAtomType].mass[OldType];

        mylnMassRatio += (3 * kbT * 0.5) * ( log(atomTypes[myISGAtomType].mass[OldType])
          - log(atomTypes[myISGAtomType].mass[NewType]) );
      } // end if statement
    } // end loop over atoms
            
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the BOND parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].bondList.begin();
         iter != myTopo->molecules[T].bondList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->bonds[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->bonds[*iter].atom2].stageNumber;

      // determine which stage this bond should be transformed in
      // (it should be the LARGER of the two stage #s, if they are different)
      int myStage =  max(atom1_stage,atom2_stage);

      // update these parameters only if myStage = 1
      if (myStage == 1) {
        // find the index # of the proper iSGBond element
        int myISGBond = myTopo->bonds[*iter].iSGmodifierIndex;
  
        // update DeltaK and DeltaR0
        myTopo->bonds[*iter].DeltaK = bonds[myISGBond].forceConstant[NewType]
          - bonds[myISGBond].forceConstant[OldType];
        myTopo->bonds[*iter].DeltaR0 = bonds[myISGBond].distance[NewType]
          - bonds[myISGBond].distance[OldType];
      }
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the ANGLE parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].angleList.begin();
         iter != myTopo->molecules[T].angleList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->angles[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->angles[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->angles[*iter].atom3].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage =  max(atom1_stage,atom2_stage,atom3_stage);

      // update these parameters only if myStage = 1
      if (myStage == 1) {
        // find the index # of the proper iSGAngle element
        int myISGAngle = myTopo->angles[*iter].iSGmodifierIndex;
  
        // update DeltaK and DeltaTheta0
        myTopo->angles[*iter].DeltaK = angles[myISGAngle].forceConstant[NewType]
          - angles[myISGAngle].forceConstant[OldType];
        myTopo->angles[*iter].DeltaTheta0 = angles[myISGAngle].angleval[NewType]
          - angles[myISGAngle].angleval[OldType];
        myTopo->angles[*iter].Delta_ubK = angles[myISGAngle].k_ub[NewType]
          - angles[myISGAngle].k_ub[OldType];
        myTopo->angles[*iter].Delta_ubR0 = angles[myISGAngle].r_ub[NewType]
          - angles[myISGAngle].r_ub[OldType];
      }
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the DIHEDRAL parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all dihedrals on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].dihedralList.begin();
      iter != myTopo->molecules[T].dihedralList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom3].stageNumber;
      int atom4_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);

      // update these parameters only if myStage = 1
      if (myStage == 1) {
        // find the index # of the proper iSGDihedral element
        int myISGDihedral = myTopo->dihedrals[*iter].iSGmodifierIndex;
  
        for (int i=0; i<myTopo->dihedrals[*iter].multiplicity; i++) {
            
          // update DeltaK and DeltaPhase
          myTopo->dihedrals[*iter].DeltaK[i] = dihedrals[myISGDihedral].forceConstant[NewType][i]
            - dihedrals[myISGDihedral].forceConstant[OldType][i];
          myTopo->dihedrals[*iter].DeltaPhase[i] = dihedrals[myISGDihedral].phaseShift[NewType][i]
            - dihedrals[myISGDihedral].phaseShift[OldType][i];
        } // end loop over i
      } // end if statement
    } // end loop over dihedrals
        
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the IMPROPER parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].improperList.begin();
         iter != myTopo->molecules[T].improperList.end(); iter++) {

      // get the transformation stage #'s for these atoms
      int atom1_stage = myTopo->atoms[myTopo->impropers[*iter].atom1].stageNumber;
      int atom2_stage = myTopo->atoms[myTopo->impropers[*iter].atom2].stageNumber;
      int atom3_stage = myTopo->atoms[myTopo->impropers[*iter].atom3].stageNumber;
      int atom4_stage = myTopo->atoms[myTopo->impropers[*iter].atom4].stageNumber;
        
      // determine which stage this bond should be transformed in
      // (it should be the LARGEST of the stage #s, if they are different)
      int myStage = max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);

      // update these parameters only if myStage = 1
      if (myStage == 1) {
        // find the index # of the proper iSGImproper element
        int myISGImproper = myTopo->impropers[*iter].iSGmodifierIndex;
  
        // update DeltaK and DeltaPhase
        myTopo->impropers[*iter].DeltaK[0] = impropers[myISGImproper].forceConstant[NewType]
          - impropers[myISGImproper].forceConstant[OldType];
        myTopo->impropers[*iter].DeltaPhase[0] = impropers[myISGImproper].phaseShift[NewType]
          - impropers[myISGImproper].phaseShift[OldType];
      }
    }

    // Compute the residual and ideal gas chemical potential differences:
    // ideal gas chemical potential difference (using the precomputed chemical potentials)
    myDMuIG = myDeltaMuIG[OldType][NewType][0];

    // residual chemical potential difference
    // DeltaMu = kT*ln(f1/f0 * N0/(N1+1))
    myTargetDeltaMu = myTargetMu[NewType] - myTargetMu[OldType]
      + kbT * log (static_cast<Real>(N[OldType]) / static_cast<Real>(N[NewType]+1));
    myTargetDeltaMu /= static_cast<Real>(myNumStages);
      
    // new Lambda velocity chosen from MB distribution
    // width of Maxwell-Boltzmann distribution (in units of [1/fs])
    Real sigma = 1.0 / myTauD;
      
    // assign a new positive lambda velocity
    Real width = randomGaussianNumber(0.0,sigma);
    width *= width;
    myLambdaVel = sqrt(width);

    //  reset the conserved quantity averages
    AveCQ = 0.0;
    AveCQSq = 0.0;
    NumTransSteps = 0;

    // reset the transformation flag
    Transformed = false;
      
    // reset AllEnergiesFile flag
    myEnergies->output(false);
    myEnergies->clear();

  } // end PickNewMolecule()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  checkForCompletion().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::checkForCompletion() {

    //report << hint << "CheckForCompletion" << endr;

    Transformed = false;

    // determine which transformation state the molecule is currently in
    int newStage = static_cast<int>(floor(myLambda)) + 1;

    // determine if the molecule has been completely transformed 
    if (myLambda > (myNumStages - 0.002) ) {
      //  add to the # of molecules of NewType and subtract from the 
      //  # of molecules of OldType
      N[NewType]++;
      myTopo->iSGNumMols[NewType]++;
      N[OldType]--;
      myTopo->iSGNumMols[OldType]--;
      
      // set Lambda for this molecule to be equal to 0 and update
      // the molecule's type
      myTopo->molecules[T].lambda = 0.0;
      myTopo->molecules[T].type = NewType;
      myTopo->molecules[T].newtype = NewType;
      FinalType = NewType;
      Transformed = true;
    }
    else if (myLambda < 0.002) {
      // set Lambda for this molecule to be equal to 0 and update
      // the molecule's type
      myTopo->molecules[T].lambda = 0.0;
      myTopo->molecules[T].type = OldType;
      myTopo->molecules[T].newtype = OldType;
      FinalType = OldType;
      Transformed = true;
    }   
    else if (newStage != thisStage) {

      // the transformation attempt has not yet completed,
      // but it has reached a new stage so we need to update
      // the force constants which are used for DeltaMu calculations

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the atomic charges and masses
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mylnMassRatio = 0.0;

      // loop over all atoms on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
           iter != myTopo->molecules[T].atoms.end(); iter++) {
      
        // the index # into the atomTypes array
        int I = myTopo->atoms[*iter].type;
        
        // the atom type name of this atom
        string myName(myTopo->atomTypes[I].name);
        
        // find the proper iSGAtomType element
        bool Found = false;
        unsigned int myISGAtomType = 0;
        while (!Found) {
          if (myName == atomTypes[myISGAtomType].name) Found = true;
          else myISGAtomType++;
        }
        
        // update deltaM, deltaQ, Qold, and Qnew
        if (myTopo->atoms[*iter].stageNumber == thisStage) myTopo->atoms[*iter].deltaM = 0.0;

        // to make sure we compute the intramolecular Coulomb force correctly, we must
        // assign the proper charge to Qold and Qnew
        myTopo->atoms[*iter].Qold = chargeMaps[myISGAtomType].old_charge[OldType][NewType][newStage-1];
        myTopo->atoms[*iter].Qnew = chargeMaps[myISGAtomType].new_charge[OldType][NewType][newStage-1];
        myTopo->atoms[*iter].deltaQ = chargeMaps[myISGAtomType].new_charge[OldType][NewType][newStage-1]
          - chargeMaps[myISGAtomType].old_charge[OldType][NewType][newStage-1];
        myTopo->atoms[*iter].Qold *= Constant::SQRTCOULOMBCONSTANT;
        myTopo->atoms[*iter].Qnew *= Constant::SQRTCOULOMBCONSTANT;          
        myTopo->atoms[*iter].deltaQ *= Constant::SQRTCOULOMBCONSTANT;

        if (myTopo->atoms[*iter].stageNumber == newStage) {
          myTopo->atoms[*iter].deltaM = atomTypes[myISGAtomType].mass[NewType]
            - atomTypes[myISGAtomType].mass[OldType];

          mylnMassRatio += (3 * kbT * 0.5) * (log(atomTypes[myISGAtomType].mass[OldType])
            - log(atomTypes[myISGAtomType].mass[NewType]));
        } // end if-else statements
      } // end loop over atoms
            
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the BOND parameters
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      // loop over all bonds on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].bondList.begin();
           iter != myTopo->molecules[T].bondList.end(); iter++) {

        // get the transformation stage #'s for these atoms
        int atom1_stage = myTopo->atoms[myTopo->bonds[*iter].atom1].stageNumber;
        int atom2_stage = myTopo->atoms[myTopo->bonds[*iter].atom2].stageNumber;

        // determine which stage this bond should be transformed in
        // (it should be the LARGER of the two stage #s, if they are different)
        int myStage =  max(atom1_stage,atom2_stage);
                
        // find the index # of the proper iSGBond element
        int myISGBond = myTopo->bonds[*iter].iSGmodifierIndex;
        
        // update DeltaK and DeltaR0
        if (myStage == thisStage) { 
          // this stage is complete, so we set the parameters to zero
          myTopo->bonds[*iter].DeltaK = 0.0;
          myTopo->bonds[*iter].DeltaR0 = 0.0;
        }
        else if (myStage == newStage) {
          // we are now starting this new stage, so these parameters must be updated
          myTopo->bonds[*iter].DeltaK = bonds[myISGBond].forceConstant[NewType]
            - bonds[myISGBond].forceConstant[OldType];
          myTopo->bonds[*iter].DeltaR0 = bonds[myISGBond].distance[NewType]
            - bonds[myISGBond].distance[OldType];
        } // end if-else statements
      } // end loop over bonds
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the ANGLE parameters
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      // loop over all angles on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].angleList.begin();
           iter != myTopo->molecules[T].angleList.end(); iter++) {

        // get the transformation stage #'s for these atoms
        int atom1_stage = myTopo->atoms[myTopo->angles[*iter].atom1].stageNumber;
        int atom2_stage = myTopo->atoms[myTopo->angles[*iter].atom2].stageNumber;
        int atom3_stage = myTopo->atoms[myTopo->angles[*iter].atom3].stageNumber;
        
        // determine which stage this bond should be transformed in
        // (it should be the LARGEST of the stage #s, if they are different)
        int myStage =  max(atom1_stage,atom2_stage,atom3_stage);
                           
        // find the index # of the proper iSGAngle element
        int myISGAngle = myTopo->angles[*iter].iSGmodifierIndex;
        
        // update DeltaK and DeltaTheta0
        if (myStage == thisStage) {
          // this stage is complete, so we set the parameters to zero
          myTopo->angles[*iter].DeltaK = 0.0;
          myTopo->angles[*iter].DeltaTheta0 = 0.0;
          myTopo->angles[*iter].Delta_ubK = 0.0;
          myTopo->angles[*iter].Delta_ubR0 = 0.0;
        }
        else if (myStage == newStage) {
          // we are now starting this new stage, so these parameters must be updated
          myTopo->angles[*iter].DeltaK = angles[myISGAngle].forceConstant[NewType]
            - angles[myISGAngle].forceConstant[OldType];
          myTopo->angles[*iter].DeltaTheta0 = angles[myISGAngle].angleval[NewType]
            - angles[myISGAngle].angleval[OldType];
          myTopo->angles[*iter].Delta_ubK = angles[myISGAngle].k_ub[NewType]
            - angles[myISGAngle].k_ub[OldType];
          myTopo->angles[*iter].Delta_ubR0 = angles[myISGAngle].r_ub[NewType]
            - angles[myISGAngle].r_ub[OldType];
        } // end if-else statements
      } // end loop over angles
      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the DIHEDRAL parameters
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      // loop over all dihedrals on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].dihedralList.begin();
           iter != myTopo->molecules[T].dihedralList.end(); iter++) {
        
        // get the transformation stage #'s for these atoms
        int atom1_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom1].stageNumber;
        int atom2_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom2].stageNumber;
        int atom3_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom3].stageNumber;
        int atom4_stage = myTopo->atoms[myTopo->dihedrals[*iter].atom4].stageNumber;
        
        // determine which stage this bond should be transformed in
        // (it should be the LARGEST of the stage #s, if they are different)
        int myStage =  max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);
                
        // find the index # of the proper iSGDihedral element
        int myISGDihedral = myTopo->dihedrals[*iter].iSGmodifierIndex;
        
        for (int i=0; i<myTopo->dihedrals[*iter].multiplicity; i++) {

          // update DeltaK and DeltaPhase
          if (myStage == thisStage) {
            // this stage is complete, so we set the parameters to zero
            myTopo->dihedrals[*iter].DeltaK[i] = 0.0;
            myTopo->dihedrals[*iter].DeltaPhase[i] = 0.0;
          }
          else if (myStage == newStage) {
            // we are now starting this new stage, so these parameters must be updated
            myTopo->dihedrals[*iter].DeltaK[i] = dihedrals[myISGDihedral].forceConstant[NewType][i]
              - dihedrals[myISGDihedral].forceConstant[OldType][i];
            myTopo->dihedrals[*iter].DeltaPhase[i] = dihedrals[myISGDihedral].phaseShift[NewType][i]
              - dihedrals[myISGDihedral].phaseShift[OldType][i];
          } // end if-else statements          
        } // end loop over i  
      } // end loop over dihedrals
    
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // update the IMPROPER parameters
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      // loop over all impropers on molecule T
      for (std::vector<int>::iterator iter = myTopo->molecules[T].improperList.begin();
           iter != myTopo->molecules[T].improperList.end(); iter++) {

        // get the transformation stage #'s for these atoms
        int atom1_stage = myTopo->atoms[myTopo->impropers[*iter].atom1].stageNumber;
        int atom2_stage = myTopo->atoms[myTopo->impropers[*iter].atom2].stageNumber;
        int atom3_stage = myTopo->atoms[myTopo->impropers[*iter].atom3].stageNumber;
        int atom4_stage = myTopo->atoms[myTopo->impropers[*iter].atom4].stageNumber;
        
        // determine which stage this bond should be transformed in
        // (it should be the LARGEST of the stage #s, if they are different)
        int myStage =  max(atom1_stage,atom2_stage,atom3_stage,atom4_stage);
                           
        // find the index # of the proper iSGImproper element
        int myISGImproper = myTopo->impropers[*iter].iSGmodifierIndex;
        

        // update DeltaK and DeltaPhase
        if (myStage == thisStage) {
          // this stage is complete, so we set the parameters to zero        
          myTopo->impropers[*iter].DeltaK[0] = 0.0;
          myTopo->impropers[*iter].DeltaPhase[0] = 0.0;
        }
        else if (myStage == newStage) {
          // we are now starting this new stage, so these parameters must be updated
          myTopo->impropers[*iter].DeltaK[0] = impropers[myISGImproper].forceConstant[NewType]
            - impropers[myISGImproper].forceConstant[OldType];
          myTopo->impropers[*iter].DeltaPhase[0] = impropers[myISGImproper].phaseShift[NewType]
            - impropers[myISGImproper].phaseShift[OldType];
        } // end if-else statements
      } // end loop over impropers
      
      // Lastly, update the transformation stage indicator and the ideal gas DeltaMu
      thisStage = newStage;
      myDMuIG = myDeltaMuIG[OldType][NewType][thisStage-1];
    }

    // if the transformation attempt has completed then compute the
    // average deviation in the conserved quantity
    if (Transformed) {     
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

      // reset the transformation stage indicator to stage 1
      thisStage = 1;

    } // end if (Transformed) statement
  } // end function

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // UpdateTopology
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void iSGIntegrator::updateTopology() {

    //report << hint << "iSGUpdateTopology" << endr;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the BOND parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].bondList.begin();
         iter != myTopo->molecules[T].bondList.end(); iter++) {
  
      // find the index # of the proper iSGBond element
      int myISGBond = myTopo->bonds[*iter].iSGmodifierIndex;
 
      // update the force constant, k
      myTopo->bonds[*iter].springConstant = bonds[myISGBond].forceConstant[FinalType];
  
      // update the rest distance, r0
      myTopo->bonds[*iter].restLength = bonds[myISGBond].distance[FinalType];
  
      // update DeltaK and DeltaR0
      myTopo->bonds[*iter].DeltaK = 0.0;
      myTopo->bonds[*iter].DeltaR0 = 0.0;
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the ANGLE parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].angleList.begin();
         iter != myTopo->molecules[T].angleList.end(); iter++) {
  
      // find the index # of the proper iSGAngle element
      int myISGAngle = myTopo->angles[*iter].iSGmodifierIndex;
  
      // update the force constant, k
      myTopo->angles[*iter].forceConstant = angles[myISGAngle].forceConstant[FinalType];
      myTopo->angles[*iter].ureyBradleyConstant = angles[myISGAngle].k_ub[FinalType];
  
      // update the rest distance, r0
      myTopo->angles[*iter].restAngle = angles[myISGAngle].angleval[FinalType];
      myTopo->angles[*iter].ureyBradleyRestLength = angles[myISGAngle].r_ub[FinalType];
  
      // update DeltaK and DeltaTheta0
      myTopo->angles[*iter].DeltaK = 0.0;
      myTopo->angles[*iter].DeltaTheta0 = 0.0;
      myTopo->angles[*iter].Delta_ubK = 0.0;
      myTopo->angles[*iter].Delta_ubR0 = 0.0;
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the DIHEDRAL parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all dihedrals on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].dihedralList.begin();
         iter != myTopo->molecules[T].dihedralList.end(); iter++) {
  
      // find the index # of the proper iSGDihedral element
      int myISGDihedral = myTopo->dihedrals[*iter].iSGmodifierIndex;
  
      for (int i=0; i<myTopo->dihedrals[*iter].multiplicity; i++) {
        // update the force constants
        myTopo->dihedrals[*iter].forceConstant[i] = dihedrals[myISGDihedral].forceConstant[FinalType][i];
   
        // update the phase shift angles
        myTopo->dihedrals[*iter].phaseShift[i] = dihedrals[myISGDihedral].phaseShift[FinalType][i];
    
        // update DeltaK and DeltaPhase
        myTopo->dihedrals[*iter].DeltaK[i] = 0.0;
        myTopo->dihedrals[*iter].DeltaPhase[i] = 0.0;
      } // end loop over i
    } // end loop over dihedrals
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the IMPROPER parameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].improperList.begin();
         iter != myTopo->molecules[T].improperList.end(); iter++) {
  
      // find the index # of the proper iSGImproper element
      int myISGImproper = myTopo->impropers[*iter].iSGmodifierIndex;
  
      // update the force constant
      myTopo->impropers[*iter].forceConstant[0] = impropers[myISGImproper].forceConstant[FinalType];
      
      // update the phase shift angle
      myTopo->impropers[*iter].phaseShift[0] = impropers[myISGImproper].phaseShift[FinalType];
  
      // update DeltaK and DeltaPhase
      myTopo->impropers[*iter].DeltaK[0] = 0.0;
      myTopo->impropers[*iter].DeltaPhase[0] = 0.0;
    }
      
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // update the atomic charges and masses
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // the new mass of the molecule
    Real newMolMass = 0.0;

    // loop over all atoms on molecule T
    for (std::vector<int>::iterator iter = myTopo->molecules[T].atoms.begin();
         iter != myTopo->molecules[T].atoms.end(); iter++) {

      // the index # into the atomTypes array
      int I = myTopo->atoms[*iter].type;

      // the atom type name of this atom
      string myName(myTopo->atomTypes[I].name);

      // find the proper iSGAtomType element
      bool Found = false;
      unsigned int myISGAtomType = 0;
      while (!Found) {
        if (myName == atomTypes[myISGAtomType].name) Found = true;
        else myISGAtomType++;
      }

      // update the mass
      myTopo->atoms[*iter].scaledMass = atomTypes[myISGAtomType].mass[FinalType];
      newMolMass += myTopo->atoms[*iter].scaledMass;

      // update the charge
      myTopo->atoms[*iter].scaledCharge = atomTypes[myISGAtomType].charge[FinalType];
      myTopo->atoms[*iter].scaledCharge *= Constant::SQRTCOULOMBCONSTANT;

      // update deltaM, deltaQ, Qold, and Qnew
      myTopo->atoms[*iter].deltaM = 0.0;
      myTopo->atoms[*iter].deltaQ = 0.0;
      myTopo->atoms[*iter].Qold = 0.0;
      myTopo->atoms[*iter].Qnew = 0.0;
    } // end loop over atoms

    // update the mass of molecule T
    myTopo->molecules[T].mass = newMolMass;
      
  } // end UpdateTopology()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // finalXSCs
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  XSC & iSGIntegrator::getXSC() const {

    // create an XSC object to store the final xsc values
    XSC *xsc = new XSC;

    // set the integrator type
    xsc->simType = "iSGVerlet";
      
    // for chemostat
    xsc->Lambda = myLambda;
    xsc->Lambda_vel = myLambdaVel;
    xsc->myMolecule = T+1;
    xsc->old_type = OldType;
    xsc->new_type = NewType;

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
