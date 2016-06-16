/* -*- c++ -*- */
#ifndef TIMEFORCE_H
#define TIMEFORCE_H

#include "Force.h"
#include "Timer.h"

namespace ProtoMol
{
	//_________________________________________________________________ TimeForce

	class TimeForce : virtual public Force
	{
		// This class contains the definition of one force

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		TimeForce(Force* actualForce);
		virtual ~TimeForce();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class TimeForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		unsigned int getIdNumber() const
		{
			return myIdNumber;
		}

	protected:
		void preprocess(unsigned int numAtoms);
		void postprocess(const GenericTopology* topo, Vector3DBlock* forces, ScalarStructure* energies);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Force
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getKeyword() const
		{
			return keyword;
		}

		virtual unsigned int numberOfBlocks(const GenericTopology* topo,
		                                    const Vector3DBlock* pos);
		virtual void uncache();
	private:
		virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const;
		virtual void getParameters(std::vector<Parameter>& parameters) const;
		virtual unsigned int getParameterSize() const;
	private:
		virtual void doSetParameters(std::string& errMsg, std::vector<Value> values);


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
	public:
		static const std::string keyword;
	protected:
		Force* myActualForce;
	private:
		unsigned int myIdNumber;
		std::string myForcename;
		static unsigned int myCounter;

		Timer myTimer;
		std::vector<TimeRep> myTimeList;
	};

	//______________________________________________________________________ INLINES
}
#endif /* TIMEFORCE_H */
