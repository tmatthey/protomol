#include "PaulTrapIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"

using std::string;
using std::vector;
using std::pair;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//________________________________________________________ ThermostatType

	const string ThermostatEnum::str[static_cast<int>(LAST) - static_cast<int>(FIRST)] = {
		// Order is essential, must be in relation to Enum
		string("undefined"), // Returned when no enum matches
		string("NVT"),
		string("NVT_zero"),
		string("NVT_ind"),
		string("NVT_shell"),
		string("NVT_global"),
		string("berendsen"),
		string("berendsen_zero"),
		string("berendsen_ind"),
		string("berendsen_shell"),
		string("berendsen_global")
	};

	//_________________________________________________________________ PaulTrapIntegrator

	const string PaulTrapIntegrator::keyword("PaulTrap");

	PaulTrapIntegrator::PaulTrapIntegrator(): STSIntegrator(),
	                                          myCached(false),
	                                          myTemperature(0.0),
	                                          myThermalInertia(0.0),
	                                          myBathPosition(0.0),
	                                          myBathVelocity(1.0)
	{
	}

	PaulTrapIntegrator::PaulTrapIntegrator(Real timestep,
	                                       Real temperature,
	                                       Real thermalInertia,
	                                       Real bathPosition,
	                                       Real bathVelocity,
	                                       ThermostatType nvttype,
	                                       Real part,
	                                       const vector<Real>& time,
	                                       const vector<Real>& t,
	                                       ForceGroup* overloadedForces)

		: STSIntegrator(timestep, overloadedForces),
		  myCached(false),
		  myTemperature(temperature),
		  myThermalInertia(thermalInertia),
		  myBathPosition(bathPosition),
		  myBathVelocity(bathVelocity),
		  myThermostatType(nvttype),
		  myPart(std::max(std::min(part, 1.0), 0.0)),
		  myPartReal(std::max(std::min(part, 1.0), 0.0)),
		  myCount(0)
	{
		// Test
		const unsigned int count = time.size();
		if (count != t.size())
			report << error << "[PaulTrapIntegrator::PaulTrapIntegrator] size of time(" << count << ") and temperature(" << t.size() << ") vectors differ!" << endr;

		// Sort
		vector<pair<Real, Real>> tmp(count);
		for (unsigned int i = 0; i < count; ++i)
		{
			tmp[i].first = time[i];
			tmp[i].second = t[i];
		}
		std::sort(tmp.begin(), tmp.end());
		myT.resize(count);
		myTime.resize(count);
		for (unsigned int i = 0; i < count; ++i)
		{
			myTime[i] = tmp[i].first;
			myT[i] = tmp[i].second;
		}
	}

	void PaulTrapIntegrator::doDrift(void)
	{
		if (!myCached)
			init();

		const unsigned int numberOfAtoms = myTopo->atoms.size();
		Real deltaT = getTimestep() / Constant::TIMEFACTOR;

		Real theTemperature = myTemperature;
		for (int i = 0; i < static_cast<int>(myT.size()) - 1; ++i)
		{
			if (myTopo->time >= myTime[i] && myTopo->time < myTime[i + 1])
			{
				theTemperature = myT[i] + (myTopo->time - myTime[i]) / (myTime[i + 1] - myTime[i]) * (myT[i + 1] - myT[i]);
				break;
			}
		}

		// Scaling of velocities 
		unsigned int countZero = 0;
		switch (myThermostatType)
		{
		case ThermostatType::NVT_SHELL:
		case ThermostatType::BERENDSEN_SHELL:
			{
				vector<Real> vavg(myLayer.size(), 0.0);
				vector<int> count(myLayer.size(), 0);
				for (unsigned int i = 0; i < myLayer.size(); i++)
				{
					for (unsigned int j = 0; j < myLayer[i].size(); j++)
					{
						unsigned int k = myLayer[i][j];
						Real rx = (*myPositions)[k].norm();
						if (rx > 1e-10)
						{
							count[i]++;
							vavg[i] += (*myVelocities)[k].dot((*myPositions)[k]) / rx;
						}
					}
				}
				for (unsigned int i = 0; i < vavg.size(); i++)
				{
					Real v = 0.0;
					if (count[i] > 0)
						v = vavg[i] / static_cast<Real>(count[i]);

					//report << hint << "Vavg["<<i<<"] ("<<myLayer[i].size()<<") ="<<v<<"."<<endr;
					for (unsigned int j = 0; j < myLayer[i].size(); j++)
					{
						unsigned int k = myLayer[i][j];
						Real rx = (*myPositions)[k].norm();
						if (rx > 1e-10)
						{
							(*myVelocities)[k] = (*myPositions)[k] * (v / rx);
						}
						else
						{
							(*myVelocities)[k] = Vector3D(0.0, 0.0, 0.0);
							countZero++;
						}
					}
				}
			}
			break;
		case ThermostatType::NVT_GLOBAL:
		case ThermostatType::BERENDSEN_GLOBAL:
			{
				Real vavg = 0.0;
				int count = 0;
				for (unsigned int i = 0; i < numberOfAtoms; i++)
				{
					if (myKeep[i] != 0)
					{
						Real rx = (*myPositions)[i].norm();
						if (rx > 1e-10)
						{
							count++;
							vavg += (*myVelocities)[i].dot((*myPositions)[i]) / rx;
						}
					}
				}
				if (count > 0)
				{
					vavg /= static_cast<Real>(count);
					//report << hint << "Vavg="<<vavg<<"."<<endr;
					for (unsigned int i = 0; i < numberOfAtoms; i++)
					{
						if (myKeep[i] != 0)
						{
							Real rx = (*myPositions)[i].norm();
							if (rx > 1e-10)
							{
								(*myVelocities)[i] = (*myPositions)[i] * (vavg / rx);
							}
							else
							{
								(*myVelocities)[i] = Vector3D(0.0, 0.0, 0.0);
								countZero++;
							}
						}
					}
				}
			}
			break;
		case ThermostatType::BERENDSEN_IND:
		case ThermostatType::NVT_IND:
			{
				for (unsigned int i = 0; i < numberOfAtoms; i++)
				{
					if (myKeep[i] != 0)
					{
						Real rx = (*myPositions)[i].normSquared();
						if (rx > 1e-10)
						{
							(*myVelocities)[i] = (*myPositions)[i] * ((*myVelocities)[i].dot((*myPositions)[i]) / rx);
						}
						else
						{
							(*myVelocities)[i] = Vector3D(0.0, 0.0, 0.0);
							countZero++;
						}
					}
				}
			}
			break;
		case ThermostatType::BERENDSEN_ZERO:
		case ThermostatType::NVT_ZERO:
			{
				for (unsigned int i = 0; i < numberOfAtoms; i++)
				{
					if (myKeep[i] != 0)
					{
						(*myVelocities)[i] = Vector3D(0.0, 0.0, 0.0);
						countZero++;
					}
				}
			}
			break;
		default:
			break;
		}


		// Integration
		switch (myThermostatType)
		{
		case ThermostatType::BERENDSEN:
		case ThermostatType::BERENDSEN_IND:
		case ThermostatType::BERENDSEN_ZERO:
		case ThermostatType::BERENDSEN_GLOBAL:
		case ThermostatType::BERENDSEN_SHELL:
			{
				//report.precision(13);
				//report << hint << myBathVelocity;
				Real t = temperature(myTopo, myVelocities);
				if (t > 0.0 && countZero < numberOfAtoms)
					t *= static_cast<Real>(numberOfAtoms) / static_cast<Real>(numberOfAtoms - countZero);
				myBathVelocity = t > 0.0 ? (myBathVelocity * (1 - myThermalInertia) + theTemperature / t * myThermalInertia) : 1.0;
				Real a = sqrt(myBathVelocity);
				//if(t < theTemperature ){
				//  a = 1.0;
				//  myBathVelocity = 1.0;
				//}
				//report <<  "a="<<a-1<<", T="<<t<<endr;
				for (unsigned int i = 0; i < numberOfAtoms; i++)
				{
					(*myVelocities)[i] *= a;
					(*myPositions)[i] += (*myVelocities)[i] * deltaT;
				}
			}
			break;
		default:
			{
				Real kineticEnergy = 0.0;
				//report.precision(13);
				//report << hint << "a="<<-myBathPosition<<", T="<<temperature(myTopo,myVelocities)<<endr;
				for (unsigned int i = 0; i < numberOfAtoms; i++)
				{
					Real m = myTopo->atoms[i].scaledMass;
					(*myVelocities)[i] -= (*myVelocities)[i] * myBathPosition;
					(*myPositions)[i] += (*myVelocities)[i] * deltaT;
					kineticEnergy += ((*myVelocities)[i]).normSquared() * m;
				}
				kineticEnergy *= 0.5;
				Real e = (3.0 / 2.0 * (Real)numberOfAtoms * Constant::BOLTZMANN * theTemperature);

				myBathPosition += (kineticEnergy - e) * myThermalInertia / (Real)numberOfAtoms;
			}
			break;
		}
		buildMolecularCenterOfMass(myPositions, myTopo);
		buildMolecularMomentum(myVelocities, myTopo);
	}

	void PaulTrapIntegrator::initialize(GenericTopology* topo,
	                                    Vector3DBlock* positions,
	                                    Vector3DBlock* velocities,
	                                    ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
		myCached = false;
		init();

		report.precision(13);
		report << plain
			<< "Paul Trap I : thermostate = " << myThermostatType << "\n"
			<< "            : part        = " << myPart << "\n";
		for (unsigned int i = 0; i < myT.size(); ++i)
			report << "            : t" << i << "          = (" << myTime[i] << "," << myT[i] << ")\n";
		report << "            : layers      = " << myLayer.size() << "\n"
			<< "            : holding     = " << myCount << " (" << 100.0 * myPartReal << "%, r=" << pow(myPartReal, 1.0 / 3.0) << ")";
		report << endr;
	}


	void PaulTrapIntegrator::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("temperature", Value(myTemperature, ConstraintValueType::NotNegative())));
		parameters.push_back(Parameter("thermal", Value(myThermalInertia)));
		parameters.push_back(Parameter("bathPos", Value(myBathPosition), 0.0));
		parameters.push_back(Parameter("bathVel", Value(myBathVelocity), 1.0));
		parameters.push_back(Parameter("scheme", Value(myThermostatType.getString(), ConstraintValueType::NotEmpty()), std::string("NVT"), Text(std::string("thermostat scheme (") + ThermostatType::getPossibleValues() + std::string(")"))));
		parameters.push_back(Parameter("part", Value(myPart, ConstraintValueType::NotNegative()), 0.0));
		parameters.push_back(Parameter("time", Value(myTime), vector<Real>()));
		parameters.push_back(Parameter("t", Value(myT, ConstraintValueType::NotNegative()), vector<Real>()));
	}

	STSIntegrator* PaulTrapIntegrator::doMake(string& errMsg, const vector<Value>& values, ForceGroup* fg) const
	{
		ThermostatType nvttype = values[5].getString();
		if (!nvttype.valid())
		{
			errMsg += " ThermostatType \'" + values[5].getString() +
				"\' not recognized, possible values are: " + ThermostatType::getPossibleValues(",") + ".";
			return NULL;
		}
		if (values[7].size() != values[8].size())
		{
			errMsg += " size of time(" + toString(values[7].size()) + ") and temperature(" + toString(values[8].size()) + ") vectors differ.";
		}
		return new PaulTrapIntegrator(values[0], values[1], values[2], values[3], values[4], nvttype, values[6], values[7], values[8], fg);
	}

	void PaulTrapIntegrator::init()
	{
		const unsigned int numberOfAtoms = myTopo->atoms.size();

		myKeep.clear();
		myKeep.resize(numberOfAtoms, 0);
		myLayer.clear();

		myCount = 0;
		if (myThermostatType != "nvt" && myThermostatType != "berendsen" && myPart > 0.0)
		{
			vector<pair<Real, unsigned int>> order(numberOfAtoms);
			for (unsigned int i = 0; i < numberOfAtoms; i++)
			{
				order[i].first = (*myPositions)[i].norm();
				order[i].second = i;
			}
			sort(order.begin(), order.end());
			Real maxr = -1.0;
			Real last = -1.0;
			for (unsigned int i = 0; i < numberOfAtoms; i++)
			{
				if (maxr >= 0.0 && fabs(order[i].first - maxr) > 1e-2)
					break;
				if (fabs(order[i].first - last) > 1e-2)
				{
					myLayer.resize(myLayer.size() + 1);
					last = order[i].first;
				}
				myLayer[myLayer.size() - 1].push_back(i);
				if (static_cast<unsigned int>(numberOfAtoms * myPart) <= i && maxr < 0.0)
				{
					maxr = order[i].first;
				}
				myKeep[order[i].second] = 1;
				myCount++;
			}
		}

		myPartReal = (numberOfAtoms > 0 ? (Real)myCount / (Real)numberOfAtoms : 0.0);
		myCached = true;
	}

	void PaulTrapIntegrator::doUncache()
	{
		myCached = false;
	}
}
