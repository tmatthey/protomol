/*  -*- c++ -*-  */
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Makeable.h"
#include "IntegratorDefinition.h"
#include "Real.h"
#include "XSC.h"

#include <string>
#include <set>

namespace ProtoMol
{
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class ForceGroup;
	class MTSIntegrator;
	class Modifier;

	//_________________________________________________________________ Integrator
	/*
	  Base class of all integrators. 
	*/

	class Integrator : public Makeable
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Integrator();
		Integrator(ForceGroup* forceGroup);
		virtual ~Integrator();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		///  Run the integrator for the certain number of timesteps.  It can
		///  be assumed that the topology and forces have not changed since 
		///  the last time the integrator was initialized.                  
		virtual void run(int numTimesteps) = 0;

		///  Set the integrator up.  This method can be called at any time, 
		///  and should perform any starting force evaluations the          
		///  integrator needs in order to run correctly.  The simulation    
		///  data should be kept in the given structures.                  
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);

		//  Needed for calculating shadow Hamiltonian.
		virtual void updateBeta(Real /*timestep*/)
		{
		}

	public:
		/// Returns the time step of the actual integrator.                 
		virtual Real getTimestep() const = 0;
		/// Returns the smallest time step of the integrator hierarchy.    
		Real getBottomTimestep() const;
		/// Forward time-integration
		void forward();
		/// Backward time-integration
		void backward();

		bool isForward() const
		{
			return myForward;
		}

	public:
		///  Store the values of the extended system coordinates at
		///  the end of a simulation into an XSC object
		virtual XSC& getXSC() const
		{
			XSC* xsc = new XSC;
			return (*xsc);
		}

	public:
		//  Returns the pointer to forces.
		Vector3DBlock* getForces() const;
		void saveForces();
		void restoreForces();

		ForceGroup* getForceGroup() const
		{
			return myForcesToEvaluate;
		}

	public:
		// Modifier objects are inserted and executed before, in-between and after
		// the force calculation and/or drift. Modifiers are automatically 
		// initialized if possible (myTopo != NULL)

		/// Add a modifier obejct at the begining of a step
		void adoptPreStepModifier(Modifier* modifier);
		/// Add a modifier object before doDriftOrNextIntegrator
		void adoptPreDriftOrNextModifier(Modifier* modifier);
		/// Add a modifier object after doDriftOrNextIntegrator
		void adoptPostDriftOrNextModifier(Modifier* modifier);
		/// Add a modifier object before the force calculation
		void adoptPreForceModifier(Modifier* modifier);
		/// Add a modifier object between system and extended force calculation
		void adoptMediForceModifier(Modifier* modifier);
		/// Add a modifier obejct after the force calculation
		void adoptPostForceModifier(Modifier* modifier);
		/// Add a modifier obejct at the end of a step
		void adoptPostStepModifier(Modifier* modifier);

		/// Delete all external modifiers
		void deleteExternalModifiers();
		/// Remove one specific modifier from the modifier list(s)
		bool removeModifier(const Modifier* modifier);
		bool removeModifier(const std::string modifierName);


		/// Create a Rattle modifier
		virtual Modifier* createRattleModifier(Real eps, int maxIter);
		/// Create a Shake modifier 
		virtual Modifier* createShakeModifier(Real eps, int maxIter);
		/// Create a Shadow modifier 
		virtual Modifier* createShadowModifier(int order2k, int freq);

	public:
		// Access methods of integrator at levels
		Integrator* top();
		const Integrator* top() const;
		Integrator* bottom();
		const Integrator* bottom() const;
		virtual Integrator* next() =0;
		virtual const Integrator* next() const =0;
		virtual Integrator* previous() =0;
		virtual const Integrator* previous() const =0;
		///  Returns the actual level of the integrator.                 
		int level() const;
		///  Returns the number of levels
		int size() const;

	public:
		/// Retrieves the integrator definition of the actual level
		IntegratorDefinition getIntegratorDefinition() const;
		/// Retrieves the complete integrator definition
		std::vector<IntegratorDefinition> getIntegratorDefinitionAll() const;

	public:
		/// Forces all integrators and their associated forces
		/// to clear the cache and pre-computed values
		void uncache();

	protected:
		// excute the modifiers
		void preStepModify();
		void preDriftOrNextModify();
		void postDriftOrNextModify();
		void preForceModify();
		void mediForceModify();
		void postForceModify();
		void postStepModify();

		/// Initialize all modifiers
		void initializeModifiers();
		/// Delete all internal modifiers
		void deleteInternalModifiers();

		/// Add modifiers which should modify during initialize
		virtual void addModifierBeforeInitialize()
		{
		}

		/// Add modifiers which should not modify during initialize
		virtual void addModifierAfterInitialize()
		{
		}

	public:
		// any test of modifiers
		bool anyPreStepModify() const
		{
			return !myPreStepModifiers.empty();
		};

		bool anyPreDriftOrNextModify() const
		{
			return !myPreDriftOrNextModifiers.empty();
		};

		bool anyPostDriftOrNextModify() const
		{
			return !myPostDriftOrNextModifiers.empty();
		};

		bool anyPreForceModify() const
		{
			return !myPreForceModifiers.empty();
		};

		bool anyMediForceModify() const
		{
			return !myMediForceModifiers.empty();
		};

		bool anyPostForceModify() const
		{
			return !myPostForceModifiers.empty();
		};

		bool anyPostStepModify() const
		{
			return !myPostStepModifiers.empty();
		};

	private:
		void addModifier(Modifier* modifier);
		void deleteModifier(Modifier* modifier);

	private:
		/// Integrator specific details of uncache
		virtual void doUncache()
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getScope() const
		{
			return scope;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string scope;
		static Real myBeta;
		Real myPotEnergy;

		//  For use in NormModeInt.  This is a hack and should be fixed.
		Real* mhQ;
		Real maxEigval;
		unsigned int numEigvects;

	protected:
		GenericTopology* myTopo;

		Vector3DBlock* myPositions;
		Vector3DBlock* myVelocities;
		Vector3DBlock* myForces;
		ScalarStructure* myEnergies;

		ForceGroup* myForcesToEvaluate;
		bool myForward;

	private:
		Vector3DBlock* myOldForces;

		std::set<Modifier*> myPreStepModifiers;
		std::set<Modifier*> myPreDriftOrNextModifiers;
		std::set<Modifier*> myPostDriftOrNextModifiers;
		std::set<Modifier*> myPreForceModifiers;
		std::set<Modifier*> myMediForceModifiers;
		std::set<Modifier*> myPostForceModifiers;
		std::set<Modifier*> myPostStepModifiers;
		std::set<Modifier*> myListModifiers;
	};

	//______________________________________________________________________ INLINES

	inline Vector3DBlock* Integrator::getForces() const
	{
		return myForces;
	}

	inline Real Integrator::getBottomTimestep() const
	{
		return (bottom()->getTimestep());
	}

	inline int Integrator::size() const
	{
		return (top()->level() + 1);
	}
}

#endif
