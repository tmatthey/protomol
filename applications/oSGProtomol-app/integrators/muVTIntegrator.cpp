//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for muVT or muVT dynamics           //
//  -----------------------------------------------------------------------  //

#include "muVTIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierNVTShake.h"
#include "ModifierNVTRattle.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Keyword.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const string muVTIntegrator::keyword("muVTVerlet");

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Default or empty constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  muVTIntegrator::muVTIntegrator(): oSGIntegrator() {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  Constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  muVTIntegrator::muVTIntegrator(Real timestep,
                                 unsigned int numComp,
			         Real temperature,
                                 Real tauT,
			         Real tauD,
                                 Real tauL,
                                 Real MuTemp,
                                 ForceGroup *overloadedForces)
: oSGIntegrator(timestep, numComp, temperature, 0.0, tauT, 0.0, 0.0, tauD, tauL, MuTemp, overloadedForces) {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  doHalfKick().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void muVTIntegrator::doHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // ---------------------------------------------------------------------
    //  Do the first update of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {
      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];

        //  Advance the velocities due to the thermostat force.
        (*myVelocities)[atom] *= exp(-myEtaVel * halfDeltaT);

      } //  end loop over atoms
    } //  end loop over molecules

    // -------------------------------------------------------------------------
    //  Do the second update of the atom velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to atomic forces.
        (*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;

      }  // end loop over atoms
    } //  end loop over molecules   
  }  //  End doHalfKick().

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  doDrift().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void muVTIntegrator::doDrift() {

    //  Timestep.  Units: (fs)
    const Real deltaT = getTimestep();

    // -------------------------------------------------------------------------
    //  Do the update of the positions
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {
      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];

        //  Advance the positions due to velocity.
        (*myPositions)[atom] += (*myVelocities)[atom] * deltaT * Constant::INV_TIMEFACTOR;

      }  //  end loop over atoms
    } //  end loop over molecules

    // calculate the new COM and momentum of each molecule
    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);
  
  }  // end doDrift()

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  do2ndHalfkick()
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void muVTIntegrator::do2ndHalfKick() {

    //  Timestep.  Units: (fs)
    const Real halfDeltaT = 0.5 * getTimestep();

    // -----------------------------------------------------------------------------
    //  Do the first update of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {
      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to atomic forces.
        (*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;

      }  //  end loop over atoms
    } //  end loop over molecules

    // ------------------------------------------------------------------------------
    //  Do the second update of the atomic velocities
    //  Loop over all molecules
    for (unsigned int i=0; i < NumMols; i++) {

      //  Temporary storage element for the updated molecular momentum
      Vector3D Momentum(0.0,0.0,0.0);

      //  Loop over the atoms on this molecule
      for (unsigned int a=0; a < myTopo->molecules[i].size(); a++) {

        //  Current atom # and mass
        int atom = myTopo->molecules[i][a];
        Real mass = myTopo->atoms[atom].scaledMass;

        //  Advance the velocities due to the thermostat force.
        (*myVelocities)[atom] *= exp(-myEtaVel * halfDeltaT);

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
  void muVTIntegrator::initialize(GenericTopology *topo,
                                 Vector3DBlock *positions,
                                 Vector3DBlock *velocities,
                                 ScalarStructure *energies) {

    oSGIntegrator::initialize(topo, positions, velocities, energies);


    //  store the total number of atoms and molecules
    NumAtoms = myTopo->atoms.size();
    NumMols = myTopo->molecules.size();
    myNumFree = myTopo->degreesOfFreedom;

    //  Initialize myVolume
    myVolume = myTopo->getVolume(*positions);

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
      myEta        = 0.;
      myEtaLambda  = 0.;
      myEtaVel     = 0.;
      myLambdaVel  = 0.;
      myEtaLambdaVel = 0.;
    }
    else {
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
      // so here we just add the V/N or V/N+1 term
      if (Insert) myTargetMu = kbT * log (myVolume / static_cast<Real>(N[OldType]+1));
      else myTargetMu -= kbT * log (myVolume / static_cast<Real>(N[OldType]));
      myTargetMu /= static_cast<Real>(myNumStages);      
    }

    //  Compute the fixed barostat mass and thermostat masses.
    Qo = myNumFree * kbT * (myTauT * myTauT);
    Qd = (Constant::BOLTZMANN * myMuTemp) * (myTauD * myTauD);
    Ql = (Constant::BOLTZMANN * myMuTemp) * (myTauL * myTauL);

    initializeForces();   
    myEnergies->output(false);

  } // end function initialize

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NVT RATTLE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifier* muVTIntegrator::createRattleModifier(Real eps, int maxIter){
    return (new ModifierNVTRattle<muVTIntegrator>(eps, maxIter, this));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //  NVT SHAKE modifier
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Modifier* muVTIntegrator::createShakeModifier(Real eps, int maxIter){
    return (new ModifierNVTShake<muVTIntegrator>(eps, maxIter, this));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // getParameters().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void muVTIntegrator::getParameters(vector< Parameter >& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("components", Value(myNumComp)));
    parameters.push_back(Parameter("temperature", Value(myTargetTemp)));
    parameters.push_back(Parameter("tauT", Value(myTauT)));
    parameters.push_back(Parameter("tauD", Value(myTauD)));
    parameters.push_back(Parameter("tauL", Value(myTauL)));
    parameters.push_back(Parameter("lambdaTemp", Value(myMuTemp)));
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // doMake().
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  STSIntegrator* muVTIntegrator::doMake(string&, const vector<Value>& values,ForceGroup* fg)const{
    return new muVTIntegrator(values[0], values[1], values[2], values[3], values[4], values[5],
                              values[6], fg );
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // finalXSCs
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  XSC & muVTIntegrator::getXSC() const {

    // create an XSC object to store the final xsc values
    XSC *xsc = new XSC;

    // set the integrator type
    xsc->simType = "muVTVerlet";

    xsc->Vol = myVolume;

    // for chemostat
    xsc->Lambda = myLambda;
    xsc->Lambda_vel = myLambdaVel;
    xsc->myMolecule = T+1;
    xsc->old_type = OldType;
    if (Insert) xsc->new_type = OldType;
    else xsc->new_type = OldType + 1;

    // for thermostats
    xsc->Eta = myEta;
    xsc->Eta_vel = myEtaVel;
    xsc->EtaLambda = myEtaLambda;
    xsc->EtaLambda_vel = myEtaLambdaVel;

    return (*xsc);
  }
  
}
