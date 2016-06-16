/*  -*- c++ -*-  */
#ifndef MODIFIER_H
#define MODIFIER_H
#include <vector>

#include "Report.h"
#include "GenericTopology.h"

namespace ProtoMol
{
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class ForceGroup;
	class Integrator;

	//________________________________________________________ Modifier
	/**
  
	   Base class for all kind of modifier implementation, based on the Observer Pattern.
	   A modifier
	   object can be added at five different stages of a single time step:@n
	   - before doDriftOrNextIntegrator (pre) @n
	   - after doDriftOrNextIntegrator (post) @n
	   - before the force calculation (pre) @n
	   - between system and extended force calculation (medi) @n
	   - after the force calculation (post) @n @n
	   Each modifier object performs it changes to the system via the method execute(...) 
	   with the implementation details in doExecute(). The execution is defined 1. by an order number
	   and 2. by the pointer of the object. 
	   An implementation can be either internal -- added by an integrator -- or
	   external -- added by the user or at application level. Integrator provides all methods
	   add, remove and delete modifiers. Integrators provides also interface to add internal
	   modifiers before/after the initialization of forces.
	   Internal modifiers (isInternal() is true) are removed, if any, by Intergrator under
	   initialize and each integrator has to add its modifiers, if any. To add internal
	   modifier you should override addModifierBeforeInitialize() and/or
	   addModifierAfterInitialize() in order the specify if the modifications should be
	   considered during the (force) initialization or not.
	   Furthermore, it is possible to disable and enable a modifier.
	*/
	class Modifier
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Modifier(int order = 0): myOrder(order),
		                         myEnable(true),
		                         myTopology(NULL),
		                         myPositions(NULL),
		                         myVelocities(NULL),
		                         myForces(NULL),
		                         myEnergies(NULL)
		{
		}

		virtual ~Modifier()
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Modifier
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		/// The method, which calls the implemenation
		void execute()
		{
			Report::report << Report::debug(10) << "[Modifier::execute] " << print() << "(" << (long)(this) << ") (enable=" << myEnable << ") at " << myTopology->time << Report::endr;
			if (myEnable)
				doExecute();
		}

		/// If the modifier is internal (added by an integrator) or
		/// external (added by the user)
		virtual bool isInternal() const =0;

		/// Returns order of execution
		int order() const
		{
			return myOrder;
		}

		/// Activate modifier
		void enable() const
		{
			myEnable = true;
		}

		/// Deactivate modifier
		void disable() const
		{
			myEnable = false;
		}

		/// If the modifier is active (doExecute() is called)
		bool isEnabled() const
		{
			return myEnable;
		}

		/// Strict weak order using first order and than pointer to use set<>
		bool operator<(const Modifier& m) const;

		/// Initialize
		void initialize(GenericTopology* topo,
		                Vector3DBlock* positions,
		                Vector3DBlock* velocities,
		                Vector3DBlock* forces,
		                ScalarStructure* energies);

		/// print/debug
		std::string print() const
		{
			return doPrint();
		}

	private:
		/// The method, which does the actual modification
		virtual void doExecute() =0;

		/// Implemenation of initialize
		virtual void doInitialize()
		{
		}

		/// Implemenation print/debug
		virtual std::string doPrint() const =0;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		int myOrder;
		mutable bool myEnable;
	protected:
		GenericTopology* myTopology;
		Vector3DBlock* myPositions;
		Vector3DBlock* myVelocities;
		Vector3DBlock* myForces;
		ScalarStructure* myEnergies;
	};

	//______________________________________________________________________ INLINES

	inline bool Modifier::operator<(const Modifier& m) const
	{
		if (myOrder < m.myOrder)
		{
			return true;
		}
		else if (myOrder > m.myOrder)
		{
			return false;
		}
		return (this < &m);
	}

	inline void Modifier::initialize(GenericTopology* topo,
	                                 Vector3DBlock* positions,
	                                 Vector3DBlock* velocities,
	                                 Vector3DBlock* forces,
	                                 ScalarStructure* energies)
	{
		myTopology = topo;
		myPositions = positions;
		myVelocities = velocities;
		myForces = forces;
		myEnergies = energies;
		Report::report << Report::debug(10) << "[Modifier::initialize] " << print() << "(" << (long)(this) << ") at " << myTopology->time << Report::endr;
		doInitialize();
	}
}
#endif /* MODIFIER_H */
