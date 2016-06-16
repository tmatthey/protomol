/*  -*- c++ -*-  */
#ifndef OUTPUTDIHEDRALS_H
#define OUTPUTDIHEDRALS_H

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "OutputFile.h"
#include "Vector3DBlock.h"


namespace ProtoMol
{
	class DCDTrajectoryWriter;

	//________________________________________________________ OutputDihedrals
	/** 
	    Writes the dihedral values to a file at given freqeuncy.
	*/
	class OutputDihedrals : public OutputFile
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OutputDihedrals();
		OutputDihedrals(const std::string& filename, int freq, int cacheFreq, int cacheSize,
		                Real closeTime, bool minimal, int index, bool dihset, std::string dsetfile);
		~OutputDihedrals();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class OutputDihedrals
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		/**
		 * PRB Feb05
		 * This function takes in a set of atom indices for select protein structures (backbone, alpha helix, etc...)
		 * and returns the Protomol indices for the matching dihedrals !the method may return extra dihedrals!
		 */
		//void getdihedralsfromatomset();

		std::vector<std::vector<Real>> dihInc(std::vector<int>);
		/**
		 * PRB Feb05
		 * This function was added to calculate a standard increment for measuring dihedral variance and assigning wells
		 */

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class OutputFile
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual void doRunCached(int step);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const;
		virtual void doInitialize();
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
			return 9;
		}

		virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;
		virtual void getParameters(std::vector<Parameter>& parameter) const;
	private:
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;

	private:
		int myDihedralIndex;
		bool myDihedralsSet;
		std::string myDihedralsSetfile;
		std::vector<int> myDihedrals;
		std::vector<std::vector<Real>> myMaximas;
		std::vector<std::vector<Real>> myMinimas;
		std::set<std::string> myConfstrings;
		DCDTrajectoryWriter* myDCD;
		std::vector<std::string> confstrings;
		std::vector<int> confstringsCounter;
		bool myMinimalImage;
		std::string oldconformstring;
		std::vector<std::vector<int>> flipcounts;
	};
}
#endif
