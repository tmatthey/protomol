/*  -*- c++ -*-  */
#ifndef PAULTRAPINTEGRATOR_H
#define PAULTRAPINTEGRATOR_H


#include "AbstractEnumType.h"
#include "STSIntegrator.h"
#include "pmconstants.h"
#include <vector>

namespace ProtoMol
{
	//_____________________________________________________ ThermostatEnum

	class ThermostatEnum
	{
	public:
		enum Enum
		{
			FIRST = 0, // Only internal purpose
			UNDEFINED = 0, // Value returned when no string matches
			NVT,
			NVT_ZERO,
			NVT_IND,
			NVT_SHELL,
			NVT_GLOBAL,
			BERENDSEN,
			BERENDSEN_ZERO,
			BERENDSEN_IND,
			BERENDSEN_SHELL,
			BERENDSEN_GLOBAL,
			LAST // Only internal purpose
		};

		static const std::string str[];
	};

	//_____________________________________________________ ThermostatType

	typedef AbstractEnumType<ThermostatEnum> ThermostatType;


	//_________________________________________________________________ PaulTrapIntegrator

	class ScalarStructure;
	class ForceGroup;

	class PaulTrapIntegrator: public STSIntegrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PaulTrapIntegrator();
		PaulTrapIntegrator(Real timestep,
		                   Real temperature,
		                   Real thermalInertia,
		                   Real bathPosition,
		                   Real bathVelocity,
		                   ThermostatType nvttype,
		                   Real part,
		                   const std::vector<Real>& time,
		                   const std::vector<Real>& t,
		                   ForceGroup* overloadedForces);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class PaulTrapIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		void init();

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
			return 9;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);


	private:
		virtual void doUncache();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class STSIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg) const;
	protected:
		virtual void doDrift();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		bool myCached;
		const Real myTemperature;
		const Real myThermalInertia;
		Real myBathPosition;
		Real myBathVelocity;
		ThermostatType myThermostatType;
		Real myPart;
		Real myPartReal;
		unsigned int myCount;
		bool myRemoveAngularMotion;
		bool myRemoveCommonMotion;
		std::vector<Real> myTime;
		std::vector<Real> myT;
		std::vector<int> myKeep;
		std::vector<std::vector<int>> myLayer;
	};

	//______________________________________________________________________ INLINES
}

#endif
