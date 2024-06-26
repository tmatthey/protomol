#include "NumericalDifferentiation.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"

using std::string;
using std::vector;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//__________________________________________________ NumericalDifferentiation

	const string NumericalDifferentiation::keyword("NumericalDifferentiation");

	NumericalDifferentiation::NumericalDifferentiation() : STSIntegrator()
	{
	}

	NumericalDifferentiation::NumericalDifferentiation(Real timestep, Real epsil,
	                                                   ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces), epsilon(epsil)
	{
		hsn.findForces(overloadedForces); //find forces and parameters
	}


	NumericalDifferentiation::~NumericalDifferentiation()
	{
	}


	void NumericalDifferentiation::initialize(GenericTopology* topo,
	                                          Vector3DBlock* positions,
	                                          Vector3DBlock* velocities,
	                                          ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
		//	
		_N = myPositions->size();
		_3N = 3 * _N;
		hsn.initialData(_3N);
		hsn.clear();
		//
	}

	void NumericalDifferentiation::run(int numTimesteps)
	{
		Vector3DBlock pmolForces, numForces, num2ndDeriv;
		Real oldPE, forcePE1, forcePE2, hessPE1, hessPE2;
		double* numHess;
		Real maxForceError = 0.0, maxHessError = 0.0;
		char coor[3] = {'x', 'y', 'z'};

		if (numTimesteps < 1)
			return;

		preStepModify();
		calculateForces();
		//Forces
		oldPE = myEnergies->potentialEnergy();
		pmolForces = *myForces;
		numForces.resize(_N);
		num2ndDeriv.resize(_N);
		//Hessian
		hsn.evaluate(myPositions, myTopo, false); //true for mass re-weight;
		numHess = new double[_3N * _3N];
		for (int k = 0; k < numTimesteps; k++)
		{
			maxForceError = 0.0;
			maxHessError = 0.0;
			if (k) epsilon /= 2.0;
			//Forces and 2nd derivative
			for (unsigned int i = 0; i < _3N; i++)
			{
				(*myPositions)[i / 3][i % 3] += epsilon;
				calculateForces();
				forcePE1 = myEnergies->potentialEnergy();
				(*myPositions)[i / 3][i % 3] -= 2.0 * epsilon;
				calculateForces();
				forcePE2 = myEnergies->potentialEnergy();
				//forces
				numForces[i / 3][i % 3] = -(forcePE1 - forcePE2) / (2.0 * epsilon);
				//2nd derivative
				num2ndDeriv[i / 3][i % 3] = (forcePE1 + forcePE2) / (epsilon * epsilon) - 2.0 * oldPE / (epsilon * epsilon);
				//
				report.precision(15);
				report << debug(2) << "[NumericalDifferentiation::run] Atom " << i / 3 << ":" << coor[i % 3] << ", force= " << pmolForces[i / 3][i % 3] <<
					", num forces= " << numForces[i / 3][i % 3] << ", epsilon= " << epsilon << "." << endr;
				(*myPositions)[i / 3][i % 3] += epsilon; //restore force positions
				//
				Real tempErr = fabs(pmolForces[i / 3][i % 3] - numForces[i / 3][i % 3]);
				if (tempErr > maxForceError) maxForceError = tempErr;
			}
			//Hessian
			for (unsigned int i = 0; i < _3N; i++)
			{
				for (unsigned int j = 0; j < _3N; j++)
				{
					(*myPositions)[i / 3][i % 3] += epsilon;
					(*myPositions)[j / 3][j % 3] += epsilon;
					calculateForces();
					hessPE1 = myEnergies->potentialEnergy();
					(*myPositions)[i / 3][i % 3] -= 2.0 * epsilon;
					(*myPositions)[j / 3][j % 3] -= 2.0 * epsilon;
					calculateForces();
					hessPE2 = myEnergies->potentialEnergy();
					numHess[i * _3N + j] = ((hessPE1 + hessPE2) / (epsilon * epsilon) - 2.0 * oldPE / (epsilon * epsilon) - num2ndDeriv[i / 3][i % 3] - num2ndDeriv[j / 3][j % 3]) * 0.5;
					report << debug(3) << "[NumericalDifferentiation::run] Atom1 " << i / 3 << ":" << coor[i % 3] << ", Atom2 " << j / 3 << ":" << coor[j % 3] << ", Hess.= " <<
						hsn.hessM[i * _3N + j] << ", num. Hess.= " << numHess[i * _3N + j] << ", epsilon= " << epsilon << "." << endr;
					(*myPositions)[i / 3][i % 3] += epsilon; //restore hessian positions
					(*myPositions)[j / 3][j % 3] += epsilon;
					Real tempErr = fabs(numHess[i * _3N + j] - hsn.hessM[i * _3N + j]);
					if (tempErr > maxForceError) maxHessError = tempErr;
				}
			}
			report << debug(1) << "[NumericalDifferentiation::run] force error = " << maxForceError << ", Hessian error = " << maxHessError << ", epsilon = " << epsilon << endr;
		}
		report << hint << "[NumericalDifferentiation::run] force error = " << maxForceError << ", Hessian error = " << maxHessError << ", epsilon = " << epsilon << endr;
		//remove storage
		delete [] numHess;
		//
		postStepModify();
	}

	void NumericalDifferentiation::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("epsilon", Value(epsilon, ConstraintValueType::Positive()), true, Text("epsilon")));
	}

	STSIntegrator* NumericalDifferentiation::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new NumericalDifferentiation(values[0], values[1], fg);
	}
}
