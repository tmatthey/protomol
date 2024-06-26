#include "RMTIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierFriction.h"

using std::string;
using std::vector;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ RMTIntegrator

	const string RMTIntegrator::keyword("RMT");

	RMTIntegrator::RMTIntegrator(): STSIntegrator(),
	                                myTemp(0.0),
	                                myQ1(0.0), myQ2(0.0), myQ3(0.0), myQ4(0.0), myQ5(0.0),
	                                myC2(0.0), myC3(0.0), myC4(0.0), myC5(0.0),
	                                myNumStats(1), tstatfile(""), filerate(0), incTdof(true)
	{
		myQ[0] = myQ1;
		myQ[1] = myQ2;
		myQ[2] = myQ3;
		myQ[3] = myQ4;
		myQ[4] = myQ5;
		myC[1] = myC2;
		myC[2] = myC3;
		myC[3] = myC4;
		myC[4] = myC5;
		if (myNumStats > 5 || myNumStats < 1) myNumStats = 1;
	}

	RMTIntegrator::RMTIntegrator(Real timestep,
	                             Real temp,
	                             Real q1, Real q2, Real q3, Real q4, Real q5,
	                             Real c2, Real c3, Real c4, Real c5,
	                             int numst, std::string tstf, int opr, bool tdof,
	                             ForceGroup* overloadedForces)

		: STSIntegrator(timestep, overloadedForces),
		  myTemp(temp),
		  myQ1(q1), myQ2(q2), myQ3(q3), myQ4(q4), myQ5(q5),
		  myC2(c2), myC3(c3), myC4(c4), myC5(c5),
		  myNumStats(numst), tstatfile(tstf), filerate(opr), incTdof(tdof)
	{
		myQ[0] = myQ1;
		myQ[1] = myQ2;
		myQ[2] = myQ3;
		myQ[3] = myQ4;
		myQ[4] = myQ5;
		myC[1] = myC2;
		myC[2] = myC3;
		myC[3] = myC4;
		myC[4] = myC5;
		if (myNumStats > 5 || myNumStats < 1) myNumStats = 1;
	}

	RMTIntegrator::~RMTIntegrator()
	{
		if (tstatfile != "") myFile.close();
	}

	void RMTIntegrator::initialize(GenericTopology* topo,
	                               Vector3DBlock* positions,
	                               Vector3DBlock* velocities,
	                               ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		//
		initializeForces();
		for (int i = 0; i < 5; i++)
		{
			myS[i] = 1.0;
			myPs[i] = 0.0;
			myAvTKE[i] = 0.0;
			myAvS[i] = 0.0;
		}
		myOldProdS = mySn = 1.0;
		//myNf = myTopo->degreesOfFreedom;
		myNf = (myPositions->size() - 2) * 3;
		if (myNf < 1) myNf = 1;
		myTopo->degreesOfFreedom = (int)myNf;
		//if(myNf >= 3) myNf += 3;
		mykT = Constant::BOLTZMANN * myTemp;
		calculateForces();
		myh0 = myEnergies->potentialEnergy() + kineticEnergy(myTopo, myVelocities);
		report << hint << "[RMTIntegrator::initialize] mykt= " << mykT << " myNf= " << myNf << " MyTemp= " << myTemp << " myH0= " << myh0 << endl;
		if (tstatfile != "")
		{ //ofstream myFile;
			myFile.open(tstatfile.c_str(), ofstream::out);
			myFile.precision(10);
		}
		totStep = 0;
		avKE = avKEsq = 0.0;
	}

	void RMTIntegrator::run(int numTimesteps)
	{
		int step;

		if (numTimesteps < 1)
			return;
		//
		preStepModify();
		//
		step = 1;
		while (step <= numTimesteps)
		{
			totStep++;
			halfUpdtH2(0);
			halfUpdtH31();
			halfUpdtH3j(0);
			UpdtH1();
			calculateForces();
			halfUpdtH3j(1);
			halfUpdtH31();
			halfUpdtH2(1);
			step += 1;
			for (int i = 0; i < myNumStats; i++)
			{
				myAvTKE[i] += 0.5 * myPs[i] * myPs[i] / myQ[i];
				myAvS[i] += myS[i];
			}
			Real ke = kineticEnergy(myTopo, myVelocities);
			avKE += ke;
			avKEsq += ke * ke;
			if (tstatfile != "" && !(totStep % filerate))
			{
				for (int i = 0; i < myNumStats; i++) myFile << myS[i] << " ";
				for (int i = 0; i < myNumStats; i++) myFile << myPs[i] << " ";
				for (int i = 0; i < myNumStats; i++) myFile << myAvS[i] / (Real)totStep << " ";
				for (int i = 0; i < myNumStats; i++) myFile << myAvTKE[i] / (Real)totStep << " ";
				Real avKE2 = avKE / (Real)totStep;
				Real varKE = avKEsq / (Real)totStep - avKE2 * avKE2;
				myFile << totalEnergy(0) << " " << avKE2 << " " << varKE << endl;
			}
		}
		//
		postStepModify();
	}

	void RMTIntegrator::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("temperature", Value(myTemp, ConstraintValueType::NotNegative()), Text("preferred system temperature")));
		parameters.push_back(Parameter("Q1", Value(myQ1), 1.0, Text("heat bath coupling 1")));
		parameters.push_back(Parameter("Q2", Value(myQ2), 1.0, Text("heat bath coupling 2")));
		parameters.push_back(Parameter("Q3", Value(myQ3), 1.0, Text("heat bath coupling 3")));
		parameters.push_back(Parameter("Q4", Value(myQ4), 1.0, Text("heat bath coupling 4")));
		parameters.push_back(Parameter("Q5", Value(myQ5), 1.0, Text("heat bath coupling 5")));
		parameters.push_back(Parameter("C2", Value(myC2), 1.0, Text("auxiliary coefficient 2")));
		parameters.push_back(Parameter("C3", Value(myC3), 1.0, Text("auxiliary coefficient 3")));
		parameters.push_back(Parameter("C4", Value(myC4), 1.0, Text("auxiliary coefficient 4")));
		parameters.push_back(Parameter("C5", Value(myC5), 1.0, Text("auxiliary coefficient 5")));
		parameters.push_back(Parameter("NumStats", Value(myNumStats), 1.0, Text("Length of chain 1-5")));
		parameters.push_back(Parameter("tstatfile", Value(tstatfile, ConstraintValueType::NoConstraints()), std::string(""), Text("Thermostat variable filename")));
		parameters.push_back(Parameter("filerate", Value(filerate), 1.0, Text("File output frequency")));
		parameters.push_back(Parameter("incTdof", Value(incTdof, ConstraintValueType::NoConstraints()), true, Text("Increment Tstat dof in method")));
	}

	STSIntegrator* RMTIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new RMTIntegrator(values[0], values[1], values[2], values[3], values[4], values[5], values[6],
		                         values[7], values[8], values[9], values[10], values[11], values[12], values[13], values[14], fg);
	}

	void RMTIntegrator::halfUpdtH2(int typ)
	{
		//H2 half step, typ=0 first cycle, typ=1 second cycle
		Real prodS, sumSpot, tempS;
		int ii, tdof;

		Real timestep = getTimestep() * Constant::INV_TIMEFACTOR;
		//  Calc Ps 
		prodS = prodSs(1, myNumStats);
		ii = 1;
		sumSpot = 0.0;
		while (ii < myNumStats)
		{
			if (incTdof) tdof = ii;
			else tdof = 0;
			sumSpot += (myNf + tdof) * mykT * log(myS[ii]) + 0.5 * (1.0 - myS[ii]) * (1.0 - myS[ii]) / myC[ii];
			ii += 1;
		}
		myPs[0] -= 0.5 * timestep * prodS * (myEnergies->potentialEnergy() + sumSpot);
		ii = 1;
		while (ii < myNumStats)
		{
			if (incTdof) tdof = ii;
			else tdof = 0;
			myPs[ii] -= 0.5 * timestep * myS[0] * (prodS / myS[ii]) * (myEnergies->potentialEnergy() + sumSpot + (myNf + tdof) * mykT - myS[ii] * (1.0 - myS[ii]) / myC[ii]);
			ii += 1;
		}
		// Calculate Momenta, test for half of cycle as using SCALED velocities
		if (typ < 1)
		{
			for (unsigned int i = 0; i < myPositions->size(); i++)
				(*myVelocities)[i] += (*myForces)[i] * timestep * 0.5 / myTopo->atoms[i].scaledMass;
			//report <<hint<<"[RMTIntegrator::initialize] timest= "<<timestep<<endl;
			// Calc current product of s'
			myOldProdS = prodS * myS[0];
		}
		else
		{
			tempS = myOldProdS / (prodS * myS[0]);
			for (unsigned int i = 0; i < myPositions->size(); i++)
				(*myVelocities)[i] = (*myVelocities)[i] * tempS + (*myForces)[i] * timestep * 0.5 / myTopo->atoms[i].scaledMass;
			//report <<hint<<"[RMTIntegrator::initialize] tempS= "<<tempS<<" timest= "<<timestep<<endl;
		}
	}

	void RMTIntegrator::halfUpdtH3j(int dir)
	{
		// H3j half step, dir=0 direction j=1...M, dir=1 direction j=M...1
		//def halfUpdtH3j(self,dir):
		int ii, kk, jj;
		Real prodSlj, prodSgj, a, c;

		Real timestep = getTimestep() * Constant::INV_TIMEFACTOR;
		kk = 1;
		while (kk < myNumStats)
		{
			if (dir > 0)
			{
				jj = myNumStats - kk;
			}
			else
			{
				jj = kk;
			}
			// 1/4 step for Ps1
			prodSlj = prodSs(0, jj);
			prodSgj = prodSs(jj + 1, myNumStats);
			a = timestep * prodSlj / (8.0 * myQ[jj] * prodSgj);
			c = -myPs[jj];
			myPs[jj] = -2.0 * c / (1.0 + sqrt(1.0 - 4.0 * a * c));
			// 1/2 step for Sj
			mySn = myS[jj];
			a = 0.25 * timestep * prodSlj * myPs[jj] / (myQ[jj] * prodSgj);
			myS[jj] *= (1.0 + a) / (1.0 - a);
			// 1/2 step for Ps2-PsM
			ii = 0;
			while (ii < jj)
			{
				myPs[ii] -= 0.25 * timestep * ((mySn + myS[jj]) / myS[ii]) * (0.5 * prodSlj * myPs[jj] * myPs[jj] / (myQ[jj] * prodSgj));
				ii += 1;
			}
			ii = jj + 1;
			while (ii < myNumStats)
			{
				myPs[ii] += 0.25 * timestep * ((mySn + myS[jj]) / myS[ii]) * (0.5 * prodSlj * myPs[jj] * myPs[jj] / (myQ[jj] * prodSgj));
				ii += 1;
			}
			//  1/4 step for Ps1
			myPs[jj] += 0.25 * timestep * (-0.5 * prodSlj * myPs[jj] * myPs[jj] / (myQ[jj] * prodSgj));
			//
			kk += 1;
		}
	}

	void RMTIntegrator::halfUpdtH31()
	{
		//H31 half step
		//def halfUpdtH31(self):
		Real prodS, a, c;
		int ii;

		Real timestep = getTimestep() * Constant::INV_TIMEFACTOR;
		// 1/4 step for Ps1
		prodS = prodSs(1, myNumStats);
		a = timestep / (8.0 * myQ[0] * prodS);
		c = -myPs[0] - 0.25 * timestep * prodS * myh0;
		myPs[0] = -2.0 * c / (1.0 + sqrt(1.0 - 4.0 * a * c));
		// 1/2 step in S1
		mySn = myS[0];
		a = 0.25 * timestep * myPs[0] / (myQ[0] * prodS);
		myS[0] *= (1.0 + a) / (1.0 - a);
		//report <<hint<<"[RMTIntegrator::initialize] S= "<<myS[0]<<" Q= "<<myQ[0]<<" numS= "<<myNumStats<<endl;
		// 1/2 step for Ps2-PsM
		ii = 1;
		while (ii < myNumStats)
		{
			myPs[ii] += 0.25 * timestep * ((mySn + myS[0]) / myS[ii]) * (0.5 * myPs[0] * myPs[0] / (myQ[0] * prodS) + prodS * myh0);
			ii += 1;
		}
		// 1/4 step for Ps1
		myPs[0] += 0.25 * timestep * (prodS * myh0 - 0.5 * myPs[0] * myPs[0] / (myQ[0] * prodS));
	}

	void RMTIntegrator::UpdtH1()
	{
		// H1 full step
		// def UpdtH1(self):
		Real prodS, tempS, sKinetic;
		int ii;

		Real timestep = getTimestep() * Constant::INV_TIMEFACTOR;
		// Positions
		prodS = prodSs(0, myNumStats);
		tempS = myOldProdS / prodS;
		for (unsigned int i = 0; i < myPositions->size(); i++)
			(*myPositions)[i] += (*myVelocities)[i] * timestep * tempS;
		// Ps'
		sKinetic = kineticEnergy(myTopo, myVelocities) * myOldProdS * myOldProdS / prodS;
		myPs[0] += (timestep / myS[0]) * (sKinetic - prodS * mykT * myNf * (1.0 + log(myS[0])));
		ii = 1;
		while (ii < myNumStats)
		{
			myPs[ii] += (timestep / myS[ii]) * (sKinetic - prodS * mykT * myNf * log(myS[0]));
			ii += 1;
		}
	}

	Real RMTIntegrator::totalEnergy(int typ)
	{
		// Calc total energy of system, typ=0 calc h0, typ=1 calc non time reparam Energy, typ=2 calc total
		// def totalEnergy(self,typ):
		Real prodS, resProdS, tempH;
		int ii, tdof;

		prodS = prodSs(0, myNumStats);
		resProdS = prodS / myS[0];
		tempH = kineticEnergy(myTopo, myVelocities) +
			myEnergies->potentialEnergy() + mykT * myNf * log(myS[0]) + 0.5 * myPs[0] * myPs[0] / (myQ[0] * resProdS * resProdS);
		ii = 1;
		while (ii < myNumStats)
		{
			if (incTdof) tdof = ii;
			else tdof = 0;
			tempH += (myNf + tdof) * mykT * log(myS[ii]) + 0.5 * (1.0 - myS[ii]) * (1.0 - myS[ii]) / myC[ii];
			resProdS /= myS[ii];
			tempH += 0.5 * myPs[ii] * myPs[ii] / (myQ[ii] * resProdS * resProdS);
			ii += 1;
		}
		if (typ > 0) tempH -= myh0;
		if (typ > 1) tempH *= prodS;
		return (tempH);
	}

	Real RMTIntegrator::prodSs(int start, int end)
	{
		int ii = start;
		Real prodS = 1.0;
		while (ii < end)
		{
			prodS *= myS[ii];
			ii += 1;
		}
		return (prodS);
	}
}
