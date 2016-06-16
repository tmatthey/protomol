/*  -*- c++ -*-  */
#ifndef OUTPUT_H
#define OUTPUT_H

#include "Makeable.h"

namespace ProtoMol
{
	class OutputCollection;
	class OutputCache;
	class Configuration;
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class Integrator;

	//________________________________________________________ Output
	/**
	   Base class of all Output classes to dump data at a given frequency. 
	   The actual output frequency is defined by global output frequency times the
	   the frequency of the concrete class. The global output frequency is also used
	   to define the number of steps to run the integrator.
	   It keeps pointers to the topology, positions, velocities, integrator and to a 
	   cache object. The cache object provides method the retrieve values such
	   that other output objects can reuse the values rather to recompute them again.
	   Output objects are aggregated in OutputCollection, which invokes them one by one at
	   application level.
	   If you only need to print some
	   output to a file you should rather inherit from OutputFile, then inherit directly from
	   Output. 
	  */
	class Output : public Makeable
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Output();
		Output(int freq);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
		                const Vector3DBlock* pos, const Vector3DBlock* vel, const ScalarStructure* energies);
		///< To initialize the object, before the simulation starts.

		void run(int step);
		//< Called at each step (e.g., printing total energy on the screen),
		//< takes care of the output frequency.

		void finalize(int step);
		//< At the end of the simulation (e.g., writing final positions), and
		//< calls first run() to ensure that run is called for the last 
		//< step, if needed.

		Output* make(std::string& errMsg, const std::vector<Value>& values) const;
		///< Factory method to create a complete output object from its prototype

		virtual bool isIdDefined(const Configuration* config) const;

		///< Should return true if the concrete object is defined/specified in 
		///< Configuration by the user. Normally if gedId() has a valid value 
		///< in Configuration.

		virtual bool addDoKeyword() const
		{
			return true;
		}

		///< Defines if the output object supports do<getId()> to enable or disable
		///< the output.

		void setCache(const OutputCache* cache);

		///< Set the cache object, such that not each output object
		///< has to re-compute same values of interest.
		///< the cache object is shared among all output objects in OutputCollection.

		int getFirstStep() const
		{
			return myFirstStep;
		}

		int getLastStep() const
		{
			return myLastStep;
		}

		int getOutputFreq() const
		{
			return myOutputFreq;
		}

		int getNext() const
		{
			return myNextStep;
		}

		bool first() const
		{
			return myFirst;
		}

		void updateNextStep(int step);

	private:
		virtual void doInitialize() =0;
		///< Hook method of initialize, implemented in the concrete class
		virtual void doRun(int step) =0;
		///< Hook method of run, implemented in the concrete class
		virtual void doFinalize(int step) =0;
		///< Hook method of finalize, implemented in the concrete class
		virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const =0;

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
	private:
		int myFirstStep;
		int myLastStep;
		int myNextStep;
		bool myFirst;
	protected:
		int myOutputFreq; ///< Output freqeuncy
		const Configuration* myConfig;
		const GenericTopology* myTopology;
		const Integrator* myIntegrator;
		const ScalarStructure* myEnergies;
		const Vector3DBlock* myPositions;
		const Vector3DBlock* myVelocities;
		const OutputCache* myCache; ///< Pointer to the shared cache object 
	};
}
#endif
