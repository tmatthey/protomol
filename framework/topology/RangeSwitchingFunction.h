/* -*- c++ -*- */
#ifndef RANGESWITCHINGFUNCTION_H
#define RANGESWITCHINGFUNCTION_H

#include <vector>

#include "Parameter.h"
#include "RangeSwitchingFunctionBase.h"

namespace ProtoMol
{
	//_________________________________________________________________ RangeSwitchingFunction

	/**
	 * Defines a range of a given from a switching function.
	 */
	template <class TOriginalSwitchingFunction>
	class RangeSwitchingFunction : private RangeSwitchingFunctionBase
	{
	public:
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Types and Enums
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		enum
		{
			USE=1
		};

		enum
		{
			MODIFY=1
		};

		enum
		{
			CUTOFF=TOriginalSwitchingFunction::CUTOFF
		};

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	public:
		RangeSwitchingFunction(): myR0(0.0), myR1(0.0)
		{
		}

		RangeSwitchingFunction(const TOriginalSwitchingFunction sw, Real r0, Real r1): myOrigFunc(sw),
		                                                                               myR0Squared(r0 * r0),
		                                                                               myR1Squared(r1 * r1),
		                                                                               myR0(r0),
		                                                                               myR1(r1)
		{
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class RangeSwitchingFunction
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		/// simple and fast test if we should apply the switching function
		bool roughTest(Real distSquared) const;

		Real cutoffSquared() const
		{
			return std::max(myR1Squared, myOrigFunc.cutoffSquared());
		}

		Real cutoff() const
		{
			return std::max(myR1, myOrigFunc.cutoff());
		}

		void operator()(Real& value, Real& derivOverD, Real distSquared) const;
		Matrix3by3 hessian(const Vector3D& rij, Real distSquared) const;

		static std::string getId()
		{
			return (keyword + TOriginalSwitchingFunction::getId());
		}

		void getParameters(std::vector<Parameter>& parameters) const;

		static unsigned int getParameterSize()
		{
			return 2 + TOriginalSwitchingFunction::getParameterSize();
		}

		static RangeSwitchingFunction make(std::string& errMsg, std::vector<Value> values);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		TOriginalSwitchingFunction myOrigFunc;
		Real myR0Squared;
		Real myR1Squared;
		Real myR0;
		Real myR1;
	};

	//______________________________________________________________________ INLINES


	template <class TOriginalSwitchingFunction>
	inline bool RangeSwitchingFunction<TOriginalSwitchingFunction>::roughTest(Real distSquared) const
	{
		return ((distSquared >= myR0Squared && distSquared < myR1Squared) && myOrigFunc.roughTest(distSquared));
	}

	template <class TOriginalSwitchingFunction>
	inline void RangeSwitchingFunction<TOriginalSwitchingFunction>::operator()(Real& value, Real& derivOverD, Real distSquared) const
	{
		if (!roughTest(distSquared))
		{
			value = 0.0;
			derivOverD = 0.0;
			return;
		}
		myOrigFunc(value, derivOverD, distSquared);
	}


	template <class TOriginalSwitchingFunction>
	inline Matrix3by3 RangeSwitchingFunction<TOriginalSwitchingFunction>::hessian(const Vector3D& rij, Real distSquared) const
	{
		return myOrigFunc.hessian(rij, distSquared);
	}


	template <class TOriginalSwitchingFunction>
	void RangeSwitchingFunction<TOriginalSwitchingFunction>::getParameters(std::vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-r0", Value(myR0, ConstraintValueType::NotNegative()), Text("range swf from")));
		parameters.push_back(Parameter("-r1", Value(myR1, ConstraintValueType::NotNegative()), Text("range swf to")));
		myOrigFunc.getParameters(parameters);
	}

	template <class TOriginalSwitchingFunction>
	RangeSwitchingFunction<TOriginalSwitchingFunction> RangeSwitchingFunction<TOriginalSwitchingFunction>::make(std::string& errMsg, std::vector<Value> values)
	{
		Real r0, r1;
		values[0].get(r0);
		values[1].get(r1);
		if (!values[0].valid() || !values[0].valid() || r1 < 0.0 || r0 < 0.0 || r0 > r1)
		{
			errMsg += keyword + " switching function: 0 <= r0 (=" + values[0].getString() + ") <= r1 (=" + values[1].getString() + ").";
			return RangeSwitchingFunction();
		}

		return RangeSwitchingFunction(TOriginalSwitchingFunction::make(errMsg, std::vector<Value>(values.begin() + 2, values.end())), r0, r1);
	}
}
#endif /* RANGESWITCHINGFUNCTION_H */
