/*  -*- c++ -*-  */
#ifndef OUTPUTXYZTRAJECTORYVEL_H
#define OUTPUTXYZTRAJECTORYVEL_H

#include "Output.h"


namespace ProtoMol
{
	class XYZTrajectoryWriter;

	//________________________________________________________ OutputXYZTrajectoryVel
	class OutputXYZTrajectoryVel : public Output
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OutputXYZTrajectoryVel();
		OutputXYZTrajectoryVel(const std::string& filename, int freq);
		virtual ~OutputXYZTrajectoryVel();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const;
		virtual void doInitialize();
		virtual void doRun(int step);
		virtual void doFinalize(int step);


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		// Returns the identification string

		virtual unsigned int getParameterSize() const
		{
			return 2;
		}

		virtual void getParameters(std::vector<Parameter>& parameter) const;
		virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		XYZTrajectoryWriter* myXYZ;
	};
}
#endif
