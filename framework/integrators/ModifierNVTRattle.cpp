#include "ModifierNVTRattle.h"
#include "Integrator.h"
#include "Topology.h"
#include "ScalarStructure.h"
#include "pmconstants.h"

using namespace ProtoMol::Report;

namespace ProtoMol {


  //__________________________________________________ ModifierNVTRattleDetails
  ModifierNVTRattleDetails::ModifierNVTRattleDetails(Real eps, int maxIter, int order):ModifierMetaRattle(eps,maxIter,order){}

  void ModifierNVTRattleDetails::doExecute(){

    // estimate the current error in all velocity constraints
    Real error = calcError();

    // delta_t
    Real dt = getTimestep() / Constant::TIMEFACTOR;

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

        // multiplicative constant due to thermostat effect
        // (for constant temperature simulations -- see Equation 56 of G. Kalibaeva,
        //  M. Ferrario, and G. Ciccotti, "Constant pressure-constant temperature molecular
        //  dynamics: a correct constrained NPT ensemble using the molecular virial", Mol. Phys.
        //  101(6), 2003, p. 765-778)
        Real Exp2_atom1 = exp(- pEta_term );
        Real Exp2_atom2 = exp(- pEta_term );
        
        // now lets compute the lambdas.
        // compute the current bond vector
        Vector3D rab = (*myPositions)[a1] - (*myPositions)[a2];
        Real rabsq = rab.normSquared();

        // compute the current velocity vector
        Vector3D vab = (*myVelocities)[a1] - (*myVelocities)[a2];   

        // dot product of distance and velocity vectors
        Real rvab = rab * vab;

        // compute the change in lambda
        Real gab = -rvab / (dt * (rM1 * Exp2_atom1 + rM2 * Exp2_atom2) * rabsq);
        Vector3D dp = rab * gab;

        // move the velocities based upon the multiplier
        (*myVelocities)[a1] += dp * dt * rM1 * Exp2_atom1;
        (*myVelocities)[a2] -= dp * dt * rM2 * Exp2_atom2;

        // the constraint adds a force to each atom since their positions
        // had to be changed.  This constraint force therefore contributes
        // to the atomic virial.  Note that the molecular virial is independent of
        // any intramolecular constraint forces.
        if (myEnergies->virial()) myEnergies->addVirial(dp*2,rab);
      }

      // compute the error in all the velocity constraints after this RATTLE iteration
      error = calcError();
      iter ++;
      if(iter > myMaxIter) {
        report << warning << "Rattle maxIter = " << myMaxIter
               << " reached, but still not converged ... error is "<<error<<endr;
        break;
      }
    }
  }

}
