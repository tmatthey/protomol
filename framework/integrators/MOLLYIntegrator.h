/*  -*- c++ -*-  */
#ifndef MOLLYINTEGRATOR_H
#define MOLLYINTEGRATOR_H

#include "MTSIntegrator.h"

namespace ProtoMol
{
	class GenericTopology;
	class ScalarStructure;
	class ForceGroup;
	class StandardIntegrator;
	class Vector3DBlock;

	//_________________________________________________________ MOLLYIntegrator

	class MOLLYIntegrator : public MTSIntegrator
	{
		friend class ModifierAveraging;
		friend class ModifierMollification;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		// totally overload the function.  Use no defaults (unless the lists passed
		//   are empty lists of forces.  Then the default forces will be assigned);
		//   however, the lists still must be passed nonetheless.
		MOLLYIntegrator();
		MOLLYIntegrator(int cycles, ForceGroup* overloadedForces, StandardIntegrator* nextIntegrator);
		virtual ~MOLLYIntegrator();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class MOLLYIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	public:
		void averagingPositions();
		void mollification();
	private:
		virtual Vector3DBlock* doAveragingPositions() = 0;
		virtual void doMollification(Vector3DBlock* doPreprocessedPositions) = 0;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);
	protected:
		/// Added modifier to call averagingPositions() and mollification()
		virtual void addModifierBeforeInitialize();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Vector3DBlock* mySwapPositions;
	};
}
#endif
