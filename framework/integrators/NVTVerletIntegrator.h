/*  -*- c++ -*-  */
#ifndef NVTVERLETINTEGRATOR_H
#define NVTVERLETINTEGRATOR_H

#include "STSIntegrator.h"

namespace ProtoMol
{
	//_______________________________________________________________ NVTVerletIntegrator
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class ForceGroup;
	template <class TIntegrator>
	class ModifierPreForceThermostat;
	template <class TIntegrator>
	class ModifierPostForceThermostat;
	class Modifer;

	class NVTVerletIntegrator: public STSIntegrator
	{
		template <class TIntegrator>
		friend class ModifierPreForceThermostat;
		template <class TIntegrator>
		friend class ModifierPostForceThermostat;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NVTVerletIntegrator();
		NVTVerletIntegrator(Real timestep,
		                    Real temperature,
		                    Real tauT,
		                    ForceGroup* overloadedForces);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NVTVerletIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		void do2ndHalfKick();
		void PreForceThermostat();
		void PostForceThermostat();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const;

		virtual unsigned int getParameterSize() const
		{
			return 3;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);
		virtual void run(int numTimesteps);

		/// Create a Rattle modifier
		virtual Modifier* createRattleModifier(Real eps, int maxIter);
		/// Create a Shake modifier 
		virtual Modifier* createShakeModifier(Real eps, int maxIter);

	protected:
		virtual void addModifierAfterInitialize();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class STSIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg) const;
	protected:
		virtual void doDrift();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class StandardIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual void doHalfKick();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;

		Real getEtaVel() const
		{
			return myEtaVel;
		}

		Real getNumAtoms() const
		{
			return NumAtoms;
		}

	private:
		const Real myTargetTemp; //  Target temperature.  Units: (K)
		const Real myTauT; //  thermostat oscillation time period.  Units: (fs)
		const Real kbT; //  Target temperature multiplied by Boltzmann's constant.  Units: (kcal/mol)
		unsigned int NumAtoms; //  Total # of atoms in the system.
		unsigned int myNumFree; //  Total # of degrees of freedom = (3*Natoms - 3) - NumConstraints
		Real Qo; //  Particle thermostat mass.  Units: (kcal fs^2 / mol)
		Real myEta; //  Nose-Hoover particle thermostat variable.  Units: (dimensionless)
		Real myEtaVel; //  Velocity of the thermostat variable.  Units: (fs)^-1
	};
}
#endif
