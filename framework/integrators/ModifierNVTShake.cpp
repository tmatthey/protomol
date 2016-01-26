#include "ModifierNVTShake.h"
#include "Integrator.h"
#include "Topology.h"
#include "ScalarStructure.h"
#include "pmconstants.h"

using namespace ProtoMol::Report;

namespace ProtoMol {


  //__________________________________________________ ModifierNVTShakeDetails
  ModifierNVTShakeDetails::ModifierNVTShakeDetails(Real eps, int maxIter, int order):ModifierMetaShake(eps,maxIter,order){}

  void ModifierNVTShakeDetails::doExecute(){
    // estimate the current error in all bond constraints
    Real error = calcError();

    // delta_t
    Real dt = getTimestep() / Constant::TIMEFACTOR;
    Real dtsq = dt * dt;

    // multiplicative constant for the thermostat velocity
    const Real pEta_term = 0.5 * getTimestep() * getEtaVel();

    int iter = 0;
    while(error > myEpsilon) {
      
      for(unsigned int i=0;i<myListOfConstraints->size();i++) {

        // find the ID#s of the two atoms in the current constraint
        int a1 = (*myListOfConstraints)[i].atom1;
        int a2 = (*myListOfConstraints)[i].atom2;

        // reciprocal atomic masses
        Real rM1 = 1/myTopology->atoms[a1].scaledMass;
        Real rM2 = 1/myTopology->atoms[a2].scaledMass;

        // multiplicative constant due to the thermostat effect
        // (for constant temperature simulations -- see Equation 56 of G. Kalibaeva,
        //  M. Ferrario, and G. Ciccotti, "Constant pressure-constant temperature molecular
        //  dynamics: a correct constrained NPT ensemble using the molecular virial", Mol. Phys.
        //  101(6), 2003, p. 765-778)
        Real Exp2_atom1 = exp(- pEta_term );
        Real Exp2_atom2 = exp(- pEta_term );

        // get the target bond distance for this constraint
        Real restLength = (*myListOfConstraints)[i].restLength;
        
        // now lets compute the lambdas.
        // compute the current bond vector
        Vector3D pab = (*myPositions)[a1] - (*myPositions)[a2];   
        Real pabsq = pab.normSquared();
        Real rabsq = restLength*restLength;
        
        // compute the difference between the target bond length and
        // the actual bond length
        Real diffsq = rabsq - pabsq; //-g^k()

        // compute the bond vector from the previous timestep
        Vector3D rab = myLastPositions[a1] - myLastPositions[a2];
        Real rpab = rab * pab;

        // calculate the constraint force, or multiplier
        Real gab = diffsq / (2 * dtsq * (rM1 + rM2) * rpab);
        Vector3D dp = rab * gab;

        // move the positions based upon the multiplier
        (*myPositions)[a1] += dp * dtsq * rM1;
        (*myPositions)[a2] -= dp * dtsq * rM2;

        // move the velocities based upon the multiplier
        (*myVelocities)[a1] += dp * dt * rM1 * Exp2_atom1;
        (*myVelocities)[a2] -= dp * dt * rM2 * Exp2_atom2;

        // the constraint adds a force to each atom since their positions
        // had to be changed.  This constraint force therefore contributes
        // to the atomic virial.  Note that the molecular virial is independent of
        // any intramolecular constraint forces.
        if (myEnergies->virial()) myEnergies->addVirial(dp*2,rab);
      }
      
      // compute the error in all the bond constraints after this SHAKE iteration
      error = calcError();
      iter ++;
      if(iter > myMaxIter) {
	report << warning << "Shake maxIter = " << myMaxIter 
	       << " reached, but still not converged ... error is "<<error<<endr;
	break;
      }      
    }

    // store the old positions
    myLastPositions = (*myPositions);
  }

}
