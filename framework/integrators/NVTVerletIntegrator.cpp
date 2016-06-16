//  -----------------------------------------------------------------------  //
//  explicit, time-reversible integrator for NVT dynamics                    //
//  -----------------------------------------------------------------------  //

#include "NVTVerletIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierPreForceThermostat.h"
#include "ModifierPostForceThermostat.h"
#include "ModifierNVTShake.h"
#include "ModifierNVTRattle.h"

using std::vector;
using std::string;

namespace ProtoMol
{
	//  -----------------------------------------------------------------------  //
	//  Keyword.
	const string NVTVerletIntegrator::keyword("NVTVerlet");

	//  -----------------------------------------------------------------------  //
	//  Default or empty constructor
	NVTVerletIntegrator::NVTVerletIntegrator(): STSIntegrator(),
	                                            myTargetTemp(0.0),
	                                            myTauT(0.0),
	                                            kbT(0.0)
	{
	}

	//  -----------------------------------------------------------------------  //
	//  Constructor
	NVTVerletIntegrator::NVTVerletIntegrator(Real timestep,
	                                         Real temperature,
	                                         Real tauT,
	                                         ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces),
		  myTargetTemp(temperature),
		  myTauT(tauT),
		  kbT(temperature * Constant::BOLTZMANN)
	{
	}

	//  -----------------------------------------------------------------------  //
	//  Thermostat -- prior to force calculations
	void NVTVerletIntegrator::PreForceThermostat()
	{
		//  Timestep.  Units: (fs)
		const Real halfDeltaT = 0.5 * getTimestep();
		//  Twice the kinetic energy.  Units: (kcal / mol)
		const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);


		//  Advance the particle thermostat variable.
		myEta += myEtaVel * halfDeltaT;
		//  Advance the particle thermostat variable velocity.
		myEtaVel += (twiceKE - myTopo->degreesOfFreedom * kbT) * halfDeltaT / Qo;
	}

	//  -----------------------------------------------------------------------  //
	//  Thermostat -- after force calculations
	void NVTVerletIntegrator::PostForceThermostat()
	{
		//  Timestep.  Units: (fs)
		const Real halfDeltaT = 0.5 * getTimestep();
		//  Get the new KE.
		const Real twiceKE = 2. * kineticEnergy(myTopo, myVelocities);


		//  Advance the particle thermostat variable velocity.
		myEtaVel += (twiceKE - myTopo->degreesOfFreedom * kbT) * halfDeltaT / Qo;
		//  Advance the particle thermostat variable.
		myEta += myEtaVel * halfDeltaT;

		//  New particle thermostat kinetic energy.
		Real pEta = (Qo * 0.5) * (myEtaVel * myEtaVel);
		//  New particle thermostat potential energy.
		Real VEta = myEta * myTopo->degreesOfFreedom * kbT;

		//  Add the energy from the extended system thermostat.
		(*myEnergies)[ScalarStructure::INTEGRATOR] += pEta + VEta;
	}

	//  -----------------------------------------------------------------------  //
	//  doHalfKick().
	void NVTVerletIntegrator::doHalfKick()
	{
		//  Timestep.  Units: (fs)
		const Real halfDeltaT = 0.5 * getTimestep();

		// ---------------------------------------------------------------------
		//  Do the first update of the atom velocities
		//  Loop over all molecules
		for (unsigned int i = 0; i < myTopo->molecules.size(); i++)
		{
			//  Temporary storage element for the updated molecular momentum
			//Vector3D Momentum(0.0,0.0,0.0);

			//  Loop over the atoms on this molecule
			for (unsigned int a = 0; a < myTopo->molecules[i].size(); a++)
			{
				//  Current atom # and mass
				int atom = myTopo->molecules[i][a];

				//  Advance the velocities due to the thermostat and barostat forces.
				(*myVelocities)[atom] *= exp(- myEtaVel * halfDeltaT);
			} //  end loop over atoms
		} //  end loop over molecules

		// -------------------------------------------------------------------------
		//  Do the second update of the atom velocities
		//  Loop over all molecules
		for (unsigned int i = 0; i < myTopo->molecules.size(); i++)
		{
			//  Loop over the atoms on this molecule
			for (unsigned int a = 0; a < myTopo->molecules[i].size(); a++)
			{
				//  Current atom # and mass
				int atom = myTopo->molecules[i][a];
				Real mass = myTopo->atoms[atom].scaledMass;

				//  Advance the velocities due to atomic forces.
				(*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;
			} // end loop over atoms
		} //  end loop over molecules
	} //  End doHalfKick().

	//  -----------------------------------------------------------------------  //
	//  doDrift().
	void NVTVerletIntegrator::doDrift()
	{
		//  Timestep.  Units: (fs)
		const Real deltaT = getTimestep();

		// -------------------------------------------------------------------------
		//  Update of the positions
		//  Loop over all molecules
		for (unsigned int i = 0; i < myTopo->molecules.size(); i++)
		{
			//  Temporary storage element for the updated molecular COM
			Vector3D COM(0.0, 0.0, 0.0);

			//  Loop over the atoms on this molecule
			for (unsigned int a = 0; a < myTopo->molecules[i].size(); a++)
			{
				//  Current atom # and mass
				int atom = myTopo->molecules[i][a];
				Real mass = myTopo->atoms[atom].scaledMass;

				//  Advance the positions due to velocity.
				(*myPositions)[atom] += (*myVelocities)[atom] * deltaT * Constant::INV_TIMEFACTOR;

				//  Add to the new COM of this molecule
				COM += (*myPositions)[atom] * mass;
			} //  end loop over atoms

			//  Store the updated molecular COM
			myTopo->molecules[i].position = COM / (myTopo->molecules[i].mass);
		} //  end loop over molecules

		//  update the COM of each molecule
		buildMolecularCenterOfMass(myPositions, myTopo);
	} // end doDrift()

	//  -------------------------------------------------------------------  //
	//  do2ndHalfkick()
	void NVTVerletIntegrator::do2ndHalfKick()
	{
		//  Timestep.  Units: (fs)
		const Real halfDeltaT = 0.5 * getTimestep();

		// -----------------------------------------------------------------------------
		//  Do the first update of the atomic velocities
		//  Loop over all molecules
		for (unsigned int i = 0; i < myTopo->molecules.size(); i++)
		{
			//  Temporary storage element for the updated molecular momentum
			Vector3D Momentum(0.0, 0.0, 0.0);

			//  Loop over the atoms on this molecule
			for (unsigned int a = 0; a < myTopo->molecules[i].size(); a++)
			{
				//  Current atom # and mass
				int atom = myTopo->molecules[i][a];
				Real mass = myTopo->atoms[atom].scaledMass;

				//  Advance the velocities due to atomic forces.
				(*myVelocities)[atom] += (*myForces)[atom] * halfDeltaT * Constant::INV_TIMEFACTOR / mass;

				//  Add to the new momentum of this molecule
				Momentum += (*myVelocities)[atom] * mass;
			} //  end loop over atoms

			//  Store the updated molecular momentum
			myTopo->molecules[i].momentum = Momentum;
		} //  end loop over molecules

		// ------------------------------------------------------------------------------
		//  Do the second and third updates of the atomic velocities
		//  Loop over all molecules
		for (unsigned int i = 0; i < myTopo->molecules.size(); i++)
		{
			//  Temporary storage element for the updated molecular momentum
			Vector3D Momentum(0.0, 0.0, 0.0);

			//  Loop over the atoms on this molecule
			for (unsigned int a = 0; a < myTopo->molecules[i].size(); a++)
			{
				//  Current atom # and mass
				int atom = myTopo->molecules[i][a];
				Real mass = myTopo->atoms[atom].scaledMass;

				//  Advance the velocities due to the thermostat forces.
				(*myVelocities)[atom] *= exp(- myEtaVel * halfDeltaT);

				//  Add to the new momentum of this molecule
				Momentum += (*myVelocities)[atom] * mass;
			} // end loop over atoms

			//  Store the updated molecular momentum
			myTopo->molecules[i].momentum = Momentum;
		} //  end loop over molecules

		//  Add the new kinetic and potential energy of the system.
		Real PE = myEnergies->potentialEnergy();
		Real KE = kineticEnergy(myTopo, myVelocities);
		(*myEnergies)[ScalarStructure::INTEGRATOR] += KE + PE;
	} // end do2ndHalfKick

	//  -----------------------------------------------------------------------  //
	//  run().
	void NVTVerletIntegrator::run(int numTimesteps)
	{
		for (int i = 0; i < numTimesteps; i++)
		{
			preStepModify();
			doHalfKick();
			doDriftOrNextIntegrator();
			calculateForces();
			do2ndHalfKick();
			postStepModify();
		}
	}

	//  -----------------------------------------------------------------------  //
	//  initialize().
	void NVTVerletIntegrator::initialize(GenericTopology* topo,
	                                     Vector3DBlock* positions,
	                                     Vector3DBlock* velocities,
	                                     ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);

		buildMolecularCenterOfMass(myPositions, myTopo);

		// initialize all forces and modifiers
		myForces->zero(positions->size());
		initializeForces();

		//  Initialize variables to something sane.
		myEta = 0.;
		myEtaVel = 0.;

		//  Compute the fixed thermostat masse.
		Qo = myTopo->degreesOfFreedom * kbT * (myTauT * myTauT);
	}

	//  -----------------------------------------------------------------------  //
	//  createRattleModifier().
	Modifier* NVTVerletIntegrator::createRattleModifier(Real eps, int maxIter)
	{
		return (new ModifierNVTRattle<NVTVerletIntegrator>(eps, maxIter, this));
	}

	//  -----------------------------------------------------------------------  //
	//  createShakeModifier().
	Modifier* NVTVerletIntegrator::createShakeModifier(Real eps, int maxIter)
	{
		return (new ModifierNVTShake<NVTVerletIntegrator>(eps, maxIter, this));
	}

	//  -----------------------------------------------------------------------  //
	//  Add the modifiers
	void NVTVerletIntegrator::addModifierAfterInitialize()
	{
		adoptPreStepModifier(new ModifierPreForceThermostat<NVTVerletIntegrator>(this, 1));
		adoptPostStepModifier(new ModifierPostForceThermostat<NVTVerletIntegrator>(this, 1));
		STSIntegrator::addModifierAfterInitialize();
	}

	//  -----------------------------------------------------------------------  //
	//  getParameters().
	void NVTVerletIntegrator::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("temperature", Value(myTargetTemp)));
		parameters.push_back(Parameter("tauT", Value(myTauT)));
	}

	//  -----------------------------------------------------------------------  //
	//  doMake().
	STSIntegrator* NVTVerletIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new NVTVerletIntegrator(values[0], values[1], values[2], fg);
	}
}
