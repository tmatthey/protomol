/*  -*- c++ -*-  */
#ifndef OUTPUTFINALXYZBINPOS_H
#define OUTPUTFINALXYZBINPOS_H

#include "Output.h"

namespace ProtoMol
{
	class Configuration;

	//________________________________________________________ OutputFinalXYZBinPos
	class OutputFinalXYZBinPos : public Output
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OutputFinalXYZBinPos();
		OutputFinalXYZBinPos(const std::string& filename, bool minimal);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const;

		virtual void doInitialize()
		{
		};

		virtual void doRun(int)
		{
		};

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
		std::string myFilename;
		bool myMinimalImage;
	};
}
#endif
