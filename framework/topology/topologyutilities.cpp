#include "topologyutilities.h"

#include "Vector3DBlock.h"
#include "GenericTopology.h"
#include "Report.h"
#include "Topology.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "stringutilities.h"
#include "ScalarStructure.h"

#include <vector>
#include <algorithm>
#include <set>

using namespace ProtoMol::Report;
using std::string;
using std::vector;
using std::set;

namespace ProtoMol
{
	//________________________________________________________________randomVelocity
	void randomVelocity(Real temperature,
	                    const GenericTopology* topology,
	                    Vector3DBlock* velocities,
	                    unsigned int seed)
	{
		// Argument tests
		if (temperature < 0)
		{
			report << error << "[randomVelocity] : "
				<< "Kelvin temperature cannot be negative. "
				<< "Aborting ...." << endr;
			return;
		}

		// Short cuts
		const Real kbT = temperature * Constant::BOLTZMANN;
		const unsigned int nAtoms = topology->atoms.size();

		Real kbToverM;


		// Make sure that the velocitie array has the right size ...
		velocities->resize(nAtoms);

		// Assign the random velocity to each atom
		for (unsigned int i = 0; i < nAtoms; i++)
		{
			kbToverM = sqrt(kbT / topology->atoms[i].scaledMass);

			// Generates a random Gaussian number with mean 0, std dev kbToverM.
			(*velocities)[i].x = kbToverM * randomGaussianNumber(seed);
			(*velocities)[i].y = kbToverM * randomGaussianNumber(seed);
			(*velocities)[i].z = kbToverM * randomGaussianNumber(seed);
		}
	}

	//________________________________________________________________randomVelocity
	void randomVelocity(Real temperatureFrom,
	                    Real temperatureTo,
	                    const GenericTopology* topology,
	                    const Vector3DBlock* positions,
	                    Vector3DBlock* velocities,
	                    bool removeLinear,
	                    bool removeAngular,
	                    unsigned int seed)
	{
		// Make sure temperatureFrom <= temperatureTo
		if (temperatureFrom > temperatureTo)
			std::swap(temperatureFrom, temperatureTo);

		// Argument tests	
		if (temperatureFrom < 0)
		{
			report << error << "[randomVelocity] : "
				<< "Kelvin temperature cannot be negative. "
				<< "Aborting ...." << endr;
			return;
		}

		// Random velocities
		randomVelocity((temperatureFrom + temperatureTo) * 0.5, topology, velocities, seed);

		if (removeLinear)
			removeLinearMomentum(velocities, topology);

		if (removeAngular)
			removeAngularMomentum(positions, velocities, topology);

		Real actual = temperature(topology, velocities);

		if (actual <= 0.0 && temperatureFrom > 0.0)
		{
			report << recoverable << "[randomVelocity] : "
				<< "Actual temperature zero, can not rescale." << endr;
			return;
		}

		// rescale if possible
		if (actual > 0.0)
		{
			Real target = temperatureFrom + (temperatureTo - temperatureFrom) * randomNumber();
			Real scale = sqrt(target / actual);
			const unsigned int nAtoms = topology->atoms.size();
			for (unsigned int i = 0; i < nAtoms; i++)
			{
				(*velocities)[i] *= scale;
			}
		}
	}

	//___________________________________________________________temperature
	Real temperature(const GenericTopology* topology,
	                 const Vector3DBlock* velocities)
	{
		return temperature(kineticEnergy(topology, velocities), topology->degreesOfFreedom);
	}

	//___________________________________________________________temperature
	Real temperature(Real kineticEnergy,
	                 unsigned int degreesOfFreedom)
	{
		return (2.0 * kineticEnergy / (Constant::BOLTZMANN * degreesOfFreedom));
	}

	//___________________________________________________________temperatureForAtomType
	Real temperatureForAtomType(const GenericTopology* topology,
	                            const Vector3DBlock* velocities,
	                            int atomType,
	                            waterOption option)
	{
		// The number of atoms
		int count;
		// The total KE
		Real KE = kineticEnergyForAtomType(topology, velocities, atomType, option, count);

		return temperature(KE, 3 * count);
	}

	//___________________________________________________________temperatureForWater
	Real temperatureForWater(const GenericTopology* topology,
	                         const Vector3DBlock* velocities)
	{
		// The number of waters
		int count;
		// The total water KE
		Real KE = kineticEnergyForWater(topology, velocities, count);

		return temperature(KE, 3 * count);
	}

	//___________________________________________________________temperatureForNonWater
	Real temperatureForNonWater(const GenericTopology* topology,
	                            const Vector3DBlock* velocities)
	{
		// The number of non-water molecules
		int count;
		// The total non-water KE
		Real KE = kineticEnergyForNonWater(topology, velocities, count);

		return temperature(KE, 3 * count);
	}

	//___________________________________________________________getNonWaterAtoms
	int getNonWaterAtoms(const GenericTopology* topology)
	{
		// The number of non-water molecules
		int atomcount = 0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->molecules[topology->atoms[i].molecule].water != true)
			{
				atomcount++;
			}
		}
		return atomcount;
	}

	//___________________________________________________________kineticEnergy
	Real kineticEnergy(const GenericTopology* topology,
	                   const Vector3DBlock* velocities)
	{
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
		}
		return (kineticEnergy * 0.5);
	}

	//___________________________________________________________kineticEnergyForAtomType
	Real kineticEnergyForAtomType(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int atomType,
	                              waterOption option)
	{
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->atoms[i].type == atomType)
			{
				bool good = false;
				if ((option == IGNORE_WATER) && (topology->molecules[topology->atoms[i].molecule].water == false))
				{
					good = true;
				}
				else if ((option == ONLY_WATER) && (topology->molecules[topology->atoms[i].molecule].water == true))
				{
					good = true;
				}
				else if (option == ALL)
				{
					good = true;
				}

				if (good == true)
				{
					kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
				}
			}
		}
		return (kineticEnergy * 0.5);
	}

	Real kineticEnergyForAtomType(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int atomType,
	                              waterOption option,
	                              int& atomCount)
	{
		Real kineticEnergy = 0.0;
		atomCount = 0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->atoms[i].type == atomType)
			{
				bool good = false;
				if ((option == IGNORE_WATER) && (topology->molecules[topology->atoms[i].molecule].water == false))
				{
					good = true;
				}
				else if ((option == ONLY_WATER) && (topology->molecules[topology->atoms[i].molecule].water == true))
				{
					good = true;
				}
				else if (option == ALL)
				{
					good = true;
				}

				if (good == true)
				{
					atomCount++;
					kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
				}
			}
		}
		return (kineticEnergy * 0.5);
	}

	//___________________________________________________________kineticEnergyForWater
	Real kineticEnergyForWater(const GenericTopology* topology,
	                           const Vector3DBlock* velocities)
	{
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->molecules[topology->atoms[i].molecule].water == true)
			{
				kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
			}
		}
		return (kineticEnergy * 0.5);
	}

	Real kineticEnergyForWater(const GenericTopology* topology,
	                           const Vector3DBlock* velocities,
	                           int& waterCount)
	{
		Real kineticEnergy = 0.0;
		waterCount = 0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->molecules[topology->atoms[i].molecule].water == true)
			{
				waterCount++;
				kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
			}
		}
		return (kineticEnergy * 0.5);
	}

	//___________________________________________________________kineticEnergyForNonWater
	Real kineticEnergyForNonWater(const GenericTopology* topology,
	                              const Vector3DBlock* velocities)
	{
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->molecules[topology->atoms[i].molecule].water == false)
			{
				kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
			}
		}
		return (kineticEnergy * 0.5);
	}

	Real kineticEnergyForNonWater(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int& nonWaterCount)
	{
		Real kineticEnergy = 0.0;
		nonWaterCount = 0;
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
		{
			if (topology->molecules[topology->atoms[i].molecule].water == false)
			{
				nonWaterCount++;
				kineticEnergy += topology->atoms[i].scaledMass * ((*velocities)[i]).normSquared();
			}
		}
		return (kineticEnergy * 0.5);
	}

	//___________________________________________________________molecularKineticEnergy
	Real molecularKineticEnergy(const GenericTopology* topology,
	                            const Vector3DBlock* velocities)
	{
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < topology->molecules.size(); i++)
		{
			// compute the COM momentum
			Vector3D tempMolMom = molecularMomentum(topology->molecules[i].atoms, velocities, topology);

			// compute twice the kinetic energy
			kineticEnergy += tempMolMom.dot(tempMolMom) / topology->molecules[i].mass;
		}

		// multiply by 1/2 to get the kinetic energy
		return (kineticEnergy * 0.5);
	}

	//___________________________________________________________velocityVirial
	ScalarStructure velocityVirial(const GenericTopology* topology,
	                               const Vector3DBlock* velocities)
	{
		ScalarStructure res;
		res.clear();
		addVelocityVirial(&res, topology, velocities);
		return res;
	}

	//___________________________________________________________addVelocityVirial
	void addVelocityVirial(ScalarStructure* energies,
	                       const GenericTopology* topology,
	                       const Vector3DBlock* velocities)
	{
		for (unsigned int i = 0; i < topology->atoms.size(); i++)
			energies->addVirial((*velocities)[i], (*velocities)[i] * topology->atoms[i].scaledMass);
	}

	//___________________________________________________________atomTypeToSymbolName
	string atomTypeToSymbolName(const string& type)
	{
		if (equalBeginNocase(type, "SOD"))
			return string("NA");
		else if (equalBeginNocase(type, "POT"))
			return string("K");
		else if (equalBeginNocase(type, "CLA"))
			return string("CL");
		else if (equalBeginNocase(type, "CAL"))
			return string("CA");
		else if (equalBeginNocase(type, "CES"))
			return string("CS");
		else if (equalBeginNocase(type, "FE"))
			return string("FE");
		else if (equalBeginNocase(type, "HE"))
			return string("HE");
		else if (equalBeginNocase(type, "NE"))
			return string("NE");
		else if (equalBeginNocase(type, "H"))
			return string("H");
		else if (equalBeginNocase(type, "C"))
			return string("C");
		else if (equalBeginNocase(type, "O"))
			return string("O");
		else if (equalBeginNocase(type, "F"))
			return string("F");
		else if (equalBeginNocase(type, "N"))
			return string("N");
		else if (equalBeginNocase(type, "P"))
			return string("P");
		else if (equalBeginNocase(type, "S"))
			return string("S");
		else if (equalBeginNocase(type, "MG"))
			return string("MG");
		else if (equalBeginNocase(type, "ZN"))
			return string("ZN");
		else
			return type;
	}


	//___________________________________________________________linearMomentum
	Vector3D linearMomentum(const Vector3DBlock* velocities,
	                        const GenericTopology* topo)
	{
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		//  Loop through all the atoms and remove the motion of center of mass
		//  using an advanced method -- Kahan's magic addition algorithm to get
		//  rid of round-off errors: Scientific Computing pp34.

		Vector3D sumMomentum((*velocities)[0] * topo->atoms[0].scaledMass);
		Vector3D tempC(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < topo->atoms.size(); i++)
		{
			Vector3D tempX((*velocities)[i] * topo->atoms[i].scaledMass);
			Vector3D tempY(tempX - tempC);
			Vector3D tempT(sumMomentum + tempY);
			tempC = (tempT - sumMomentum) - tempY;
			sumMomentum = tempT;
		}
		return sumMomentum;
	}

	//___________________________________________________________linearMomentumSolute
	Vector3D linearMomentumSolute(const Vector3DBlock* velocities,
	                              const GenericTopology* topo)
	{
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		//  Loop through all the solute atoms and remove the motion of center of mass
		//  using an advanced method -- Kahan's magic addition algorithm to get
		//  rid of round-off errors: Scientific Computing pp34.

		Vector3D sumMomentum((*velocities)[0] * topo->atoms[0].scaledMass);
		Vector3D tempC(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < topo->atoms.size(); i++)
		{
			//cout << "The atom type: " << topo->atoms[i].type  << endl;
			//cout << "The molecule number: " << topo->atoms[i].molecule  << endl;
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				//cout << "The atom type (not water): " << topo->atoms[i].type  << endl;
				//cout << "The molecule number: " << topo->atoms[i].molecule  << endl;
				Vector3D tempX((*velocities)[i] * topo->atoms[i].scaledMass);
				Vector3D tempY(tempX - tempC);
				Vector3D tempT(sumMomentum + tempY);
				tempC = (tempT - sumMomentum) - tempY;
				sumMomentum = tempT;
			}
		}
		return sumMomentum;
	}


	//___________________________________________________________removeLinearMomentum
	Vector3D removeLinearMomentum(Vector3DBlock* velocities,
	                              const GenericTopology* topo)
	{
		// UPDATE FRI APR 13th 2007 by PRB
		// ONLY REMOVES LINEAR MOMENTUM OF THE SOLUTE
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		//Vector3D momentum     = linearMomentum(velocities,topo);
		Vector3D momentum = linearMomentumSolute(velocities, topo);

		// Count number of solute atoms
		int soluteAtoms = 0;
		for (unsigned int i = 0; i < topo->atoms.size(); i++)
		{
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				soluteAtoms += 1;
			}
		}

		//Vector3D avgMomentum  = momentum/velocities->size();
		Vector3D avgMomentum = momentum / soluteAtoms;
		for (unsigned int i = 0; i < topo->atoms.size(); i++)
		{
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				(*velocities)[i] -= avgMomentum / topo->atoms[i].scaledMass;
			}
		}

		return momentum;
	}

	//___________________________________________________________centerOfMass
	Vector3D centerOfMass(const Vector3DBlock* positions,
	                      const GenericTopology* topo)
	{
		if (positions->empty())
			return Vector3D(0.0, 0.0, 0.0);

		// Loop through all the atoms and compute center of mass
		// using an advanced method -- Kahan's magic addition algorithm to get
		// rid of round-off errors: Scientific Computing pp34.

		Real sumM = 0.0;
		Real m = topo->atoms[0].scaledMass;
		sumM += m;
		Vector3D centerOfMass((*positions)[0] * m);
		Vector3D tempC(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < topo->atoms.size(); i++)
		{
			m = topo->atoms[i].scaledMass;
			sumM += m;
			Vector3D tempX((*positions)[i] * m);
			Vector3D tempY(tempX - tempC);
			Vector3D tempT(centerOfMass + tempY);
			tempC = (tempT - centerOfMass) - tempY;
			centerOfMass = tempT;
		}
		return centerOfMass / sumM;
	}

	//___________________________________________________________angularMomentum
	Vector3D angularMomentum(const Vector3DBlock* positions,
	                         const Vector3DBlock* velocities,
	                         const GenericTopology* topo)
	{
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		return angularMomentum(positions, velocities, topo, centerOfMass(positions, topo));
	}

	//___________________________________________________________angularMomentum
	Vector3D angularMomentum(const Vector3DBlock* positions,
	                         const Vector3DBlock* velocities,
	                         const GenericTopology* topo,
	                         const Vector3D& centerOfMass)
	{
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		// Loop through all the atoms and compute angular momentum
		// using an advanced method -- Kahan's magic addition algorithm to get
		// rid of round-off errors: Scientific Computing pp34.

		Vector3D momentum(((*positions)[0] - centerOfMass).cross((*velocities)[0]) * topo->atoms[0].scaledMass);
		Vector3D tempC = Vector3D(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < topo->atoms.size(); i++)
		{
			Vector3D tempX(((*positions)[i] - centerOfMass).cross((*velocities)[i]) * topo->atoms[i].scaledMass);
			Vector3D tempY(tempX - tempC);
			Vector3D tempT(momentum + tempY);
			tempC = (tempT - momentum) - tempY;
			momentum = tempT;
		}
		return momentum;
	}

	//___________________________________________________________angularMomentumSolute
	Vector3D angularMomentumSolute(const Vector3DBlock* positions,
	                               const Vector3DBlock* velocities,
	                               const GenericTopology* topo,
	                               const Vector3D& centerOfMass)
	{
		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		// Loop through all the atoms and compute angular momentum
		// using an advanced method -- Kahan's magic addition algorithm to get
		// rid of round-off errors: Scientific Computing pp34.

		Vector3D momentum(((*positions)[0] - centerOfMass).cross((*velocities)[0]) * topo->atoms[0].scaledMass);
		Vector3D tempC = Vector3D(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < topo->atoms.size(); i++)
		{
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				Vector3D tempX(((*positions)[i] - centerOfMass).cross((*velocities)[i]) * topo->atoms[i].scaledMass);
				Vector3D tempY(tempX - tempC);
				Vector3D tempT(momentum + tempY);
				tempC = (tempT - momentum) - tempY;
				momentum = tempT;
			}
		}
		return momentum;
	}

	//___________________________________________________________inertiaMomentum
	Matrix3by3 inertiaMomentum(const Vector3DBlock* positions,
	                           const GenericTopology* topo,
	                           const Vector3D& centerOfMass)
	{
		Matrix3by3 inertia;
		inertia.zeroMatrix();
		for (unsigned int i = 0; i < topo->atoms.size(); i++)
		{
			Vector3D r((*positions)[i] - centerOfMass);
			Real x = r.x;
			Real y = r.y;
			Real z = r.z;
			inertia += Matrix3by3(y * y + z * z, -x * y, -x * z,
			                      -x * y, x * x + z * z, -y * z,
			                      -x * z, -y * z, x * x + y * y) * topo->atoms[i].scaledMass;
		}
		return inertia;
	}

	//___________________________________________________________inertiaMomentumSolute
	Matrix3by3 inertiaMomentumSolute(const Vector3DBlock* positions,
	                                 const GenericTopology* topo,
	                                 const Vector3D& centerOfMass)
	{
		Matrix3by3 inertia;
		inertia.zeroMatrix();
		for (unsigned int i = 0; i < topo->atoms.size(); i++)
		{
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				Vector3D r((*positions)[i] - centerOfMass);
				Real x = r.x;
				Real y = r.y;
				Real z = r.z;
				inertia += Matrix3by3(y * y + z * z, -x * y, -x * z,
				                      -x * y, x * x + z * z, -y * z,
				                      -x * z, -y * z, x * x + y * y) * topo->atoms[i].scaledMass;
			}
		}
		return inertia;
	}

	//___________________________________________________________removeAngularMomentum
	Vector3D removeAngularMomentum(const Vector3DBlock* positions,
	                               Vector3DBlock* velocities,
	                               const GenericTopology* topo)
	{
		// UPDATE FRI APR 13th 2007 by PRB
		// ONLY REMOVES ANGULAR MOMENTUM OF THE SOLUTE

		if (velocities->empty())
			return Vector3D(0.0, 0.0, 0.0);

		// Center of mass
		Vector3D center = centerOfMass(positions, topo);

		// Angular momentum
		//Vector3D momentum = angularMomentum(positions,velocities,topo,center);
		Vector3D momentum = angularMomentumSolute(positions, velocities, topo, center);
		if (momentum.normSquared() == 0)
		{
			return Vector3D(0.0, 0.0, 0.0);
		}

		// Inertia moment
		// Matrix3by3 inertia = inertiaMomentum(positions,topo,center);
		Matrix3by3 inertia = inertiaMomentumSolute(positions, topo, center);

		if (!inertia.invert())
		{
			return Vector3D(0.0, 0.0, 0.0);
		}
		// Remove angular momentum
		Vector3D w(inertia * momentum);
		for (unsigned int i = 0; i < topo->atoms.size(); i++)
		{
			if (topo->molecules[topo->atoms[i].molecule].water == false)
			{
				Vector3D r((*positions)[i] - center);
				(*velocities)[i] -= w.cross(r);
			}
		}

		return momentum;
	}


	//___________________________________________________________computePressure
	Real computePressure(const GenericTopology* topology,
	                     const Vector3DBlock* positions,
	                     const Vector3DBlock* velocities,
	                     const ScalarStructure* energies)
	{
		return computePressure(energies, topology->getVolume(*positions), kineticEnergy(topology, velocities));
	}

	//___________________________________________________________computePressure
	Real computePressure(const ScalarStructure* energies,
	                     Real volume,
	                     Real kineticEnergy)
	{
		return (energies->pressure(volume) + kineticEnergy * 2.0 / 3.0 / volume * Constant::PRESSUREFACTOR);
	}

	//___________________________________________________________computeMolecularPressure
	Real computeMolecularPressure(const ScalarStructure* energies,
	                              Real volume,
	                              Real kineticEnergy)
	{
		return (energies->molecularPressure(volume) + kineticEnergy * 2.0 / 3.0 / volume * Constant::PRESSUREFACTOR);
	}


	//___________________________________________________________molecularMomentum
	Vector3D molecularMomentum(const vector<int>& atomList,
	                           const Vector3DBlock* velocities,
	                           const GenericTopology* topo)
	{
		//  Loop through all the atoms and remove the motion of center of mass
		//  using an advanced method -- Kahan's magic addition algorithm to get
		//  rid of round-off errors: Scientific Computing pp34.

		Vector3D sumMomentum((*velocities)[atomList[0]] * topo->atoms[atomList[0]].scaledMass);
		Vector3D tempC(0.0, 0.0, 0.0);
		for (unsigned int i = 1; i < atomList.size(); i++)
		{
			Vector3D tempX((*velocities)[atomList[i]] * topo->atoms[atomList[i]].scaledMass);
			Vector3D tempY(tempX - tempC);
			Vector3D tempT(sumMomentum + tempY);
			tempC = (tempT - sumMomentum) - tempY;
			sumMomentum = tempT;
		}
		return sumMomentum;
	}

	//___________________________________________________________molecularCenterOfMass
	Vector3D molecularCenterOfMass(const vector<int>& atomList,
	                               const Vector3DBlock* positions,
	                               const GenericTopology* topo)
	{
		if (atomList.empty())
			return Vector3D(0.0, 0.0, 0.0);

		// Loop through all the atoms  and compute center of mass using an advanced
		// method -- Kahan's magic addition algorithm to get rid of round-off
		// errors: Heath, "Scientific Computing", pp48.

		const int numAtoms = atomList.size();

		Real mass = topo->atoms[atomList[0]].scaledMass;
		Real sumMass = mass;
		Vector3D center((*positions)[atomList[0]] * mass);
		Vector3D tempC(0.0, 0.0, 0.0);

		for (int i = 1; i < numAtoms; i++)
		{
			mass = topo->atoms[atomList[i]].scaledMass;
			sumMass += mass;
			Vector3D tempX((*positions)[atomList[i]] * mass);
			Vector3D tempY(tempX - tempC);
			Vector3D tempT(center + tempY);
			tempC = (tempT - center) - tempY;
			center = tempT;
		}

		return (center / sumMass);
	}

	//___________________________________________________________buildMolecularCenterOfMass
	void buildMolecularCenterOfMass(const Vector3DBlock* positions,
	                                GenericTopology* topo)
	{
		for (unsigned int i = 0; i < topo->molecules.size(); ++i)
		{
			topo->molecules[i].position = molecularCenterOfMass(topo->molecules[i].atoms, positions, topo);
		}
	}

	//___________________________________________________________buildMolecularMomentum
	void buildMolecularMomentum(const Vector3DBlock* velocities,
	                            GenericTopology* topo)
	{
		for (unsigned int i = 0; i < topo->molecules.size(); ++i)
		{
			topo->molecules[i].momentum = molecularMomentum(topo->molecules[i].atoms, velocities, topo);
		}
	}


	void buildRattleShakeBondConstraintList(GenericTopology* topology, vector<Bond::Constraint>& bondConstraints)
	{
		// here we go through the bond list first, then the angle list to extract the possible
		// third pair if they are both hydrogen. Thus, all bond lengths plus the third pair in 
		// the waters will be constrained. The angles rather than H-(heavy)-H are left free 
		// of vibrations 

		bondConstraints.clear();

		for (unsigned int i = 0; i < topology->bonds.size(); i++)
		{
			bondConstraints.push_back(Bond::Constraint(topology->bonds[i].atom1, topology->bonds[i].atom2, topology->bonds[i].restLength));
		}

		// now we have used all the spaces allocated for bondConstraints.
		// we will use push_back method for more constraints.    
		for (unsigned int i = 0; i < topology->angles.size(); i++)
		{
			int a1 = topology->angles[i].atom1;
			int a2 = topology->angles[i].atom2;
			int a3 = topology->angles[i].atom3;

			PairIntSorted p1(a1, a2);
			PairIntSorted p2(a2, a3);

			Real M1 = topology->atoms[a1].scaledMass;
			Real M3 = topology->atoms[a3].scaledMass;
			// For regular water M3 = M1 = 1.008, and for heavy waters, 
			// M1 or M3 should be 2 ~ 3 times heavier
			if ((M1 < 5) && (M3 < 5))
			{
				// this exclude (heavy atom)-H pairs and (heavy)-(heavy) pairs
				//           sqrt( 2 * 0.957 * 0.957 * ( 1 - cos( 104.52 * 0.017453 ) ) );
				//bondConstraints.push_back(Bond::Constraint(a1,a3,1.51356642665));

				// ... properly retrieve the right length!
				int b1 = -1;
				int b2 = -1;
				const std::vector<int>& mybonds1 = topology->atoms[a1].mybonds;
				const std::vector<int>& mybonds3 = topology->atoms[a3].mybonds;

				for (unsigned int j = 0; j < mybonds1.size(); j++)
				{
					PairIntSorted s1(topology->bonds[mybonds1[j]].atom1, topology->bonds[mybonds1[j]].atom2);
					if (p1 == s1)
					{
						for (unsigned int k = 0; k < mybonds3.size(); k++)
						{
							PairIntSorted s2(topology->bonds[mybonds3[k]].atom1, topology->bonds[mybonds3[k]].atom2);
							if (p2 == s2)
							{
								if (b1 > -1 || b2 > -1)
								{
									report << warning << "Found already two bonds (" << b1 << "," << b2 << ") for angle " << i << "." << endr;
								}
								else
								{
									b1 = mybonds1[j];
									b2 = mybonds3[k];
								}
							}
						}
					}
				}
				if (b1 > -1 && b2 > -1)
				{
					Real b = topology->bonds[b1].restLength;
					Real c = topology->bonds[b2].restLength;
					Real alpha = topology->angles[i].restAngle;
					bondConstraints.push_back(Bond::Constraint(a1, a3, sqrt(b * b + c * c - 2 * b * c * cos(alpha))));
					report << debug(10) << i << " \t :b=" << b << ", c=" << c << ", alpha=" << alpha << ", H-H r=" << sqrt(b * b + c * c - 2 * b * c * cos(alpha)) << endr;
				}
				else
				{
					report << debug(10) << "Could not find two matching bonds for angle " << i << "." << endr;
				}
			}
		}
		if (bondConstraints.size() == topology->bonds.size())
			report << hint << "No additional H-X-H SHAKE/RATTLE constraint contributions." << endr;
	}

	//___________________________________________________________getAtomsBondedtoDihedral
	// this function gets all the atoms bonded to ONE side of the dihedral
	void getAtomsBondedtoDihedral(const GenericTopology* topology,
	                              set<int, std::less<int>>* atomSet,
	                              const int atomID,
	                              const int inAtomID,
	                              const int outAtomID,
	                              const int exclAtomID)
	{
		atomSet->insert(atomID);

		for (unsigned int j = 0; j < topology->atoms[atomID].mybonds.size(); j++)
		{
			int bondindex = topology->atoms[atomID].mybonds[j];
			int atom1 = topology->bonds[bondindex].atom1;
			int atom2 = topology->bonds[bondindex].atom2;

			//MUST IMPLEMENT ERROR CASE FOR A BOND LOOP THROUGH THE DIHEDRAL

			if (!((atom1 == inAtomID && atom2 == outAtomID) || (atom1 == outAtomID && atom2 == inAtomID)))
			{
				if (atom1 == atomID)
				{
					if (atomSet->find(atom2) == atomSet->end())
					{
						getAtomsBondedtoDihedral(topology, atomSet, atom2, inAtomID, outAtomID, exclAtomID);
					}
				}
				else
				{
					if (atomSet->find(atom1) == atomSet->end())
					{
						getAtomsBondedtoDihedral(topology, atomSet, atom1, inAtomID, outAtomID, exclAtomID);
					}
				}
			}
		}
	}

	void rotateDihedral(const GenericTopology* topology,
	                    Vector3DBlock* positions,
	                    Vector3DBlock* velocities,
	                    const int dihedralID, Real angle)
	{
		int exclusionAtomID = topology->dihedrals[dihedralID].atom1;
		int innerAtomID1 = topology->dihedrals[dihedralID].atom2;
		int innerAtomID2 = topology->dihedrals[dihedralID].atom3;
		int atomID = topology->dihedrals[dihedralID].atom4;
		vector<AngleInfo> angles;

		build_angle_list(topology, atomID, innerAtomID1, innerAtomID2, exclusionAtomID, angle, &angles);
		general_rotation(innerAtomID1, innerAtomID2, positions, velocities, &angles);
	}

	void rotateDihedral(const GenericTopology* topology,
	                    Vector3DBlock* positions,
	                    const int dihedralID, Real angle)
	{
		int exclusionAtomID = topology->dihedrals[dihedralID].atom1;
		int innerAtomID1 = topology->dihedrals[dihedralID].atom2;
		int innerAtomID2 = topology->dihedrals[dihedralID].atom3;
		int atomID = topology->dihedrals[dihedralID].atom4;
		vector<AngleInfo> angles;

		build_angle_list(topology, atomID, innerAtomID1, innerAtomID2, exclusionAtomID, angle, &angles);
		general_rotation(innerAtomID1, innerAtomID2, positions, &angles);
	}

	void set_angles(Stack<unsigned int>* nodeStack, vector<AngleInfo>* angles, bool lastIsInnerAtom, Real wholeAngle)
	{
		Real startAngle = 0.0;
		Real finishAngle = 0.0;
		bool done = false;
		unsigned int atomID;
		unsigned int pos;
		unsigned int count;

		count = nodeStack->getNumElements();
		finishAngle = (*angles)[nodeStack->getElement(count - 1)].getAngle();
		if (lastIsInnerAtom)
		{
			finishAngle = wholeAngle;
		}
		pos = nodeStack->getNumElements() - 2;
		while (!done)
		{
			atomID = nodeStack->getElement(pos);
			if ((*angles)[atomID].getAngleType() == AngleInfo::ANGLE_POINTER)
			{
				pos = pos - 1;
			}
			else if ((*angles)[atomID].getAngleType() == AngleInfo::ANGLE_VALUE)
			{
				startAngle = (*angles)[atomID].getAngle();
				done = true;
			}
			else if (pos > nodeStack->getNumElements() - 2)
			{
				report << "Error: pos out of range!!!" << endr;
			}
			else
			{
				report << "Error: No specified Angle Type!!!" << endr;
			}
		}

		count = count - pos - 1;
		finishAngle = startAngle - finishAngle;

		finishAngle = finishAngle / count;

		for (unsigned int x = pos; x < nodeStack->getNumElements(); x++)
		{
			atomID = nodeStack->getElement(x);
			(*angles)[atomID].setAngle(startAngle);
			startAngle = startAngle - finishAngle;
		}
	}


	void build_angle_list(const GenericTopology* topo, const unsigned int atomID, const unsigned int inAtomID, const unsigned int outAtomID, const unsigned int exclAtomID, Real rotAngle, vector<AngleInfo>* angles)
	{
		Stack<unsigned int> nodeStack, indexStack;
		unsigned int curElement;
		unsigned int index;
		unsigned int tempBond1, tempBond2;
		bool done;
		AngleInfo temp;

		for (unsigned int x = 0; x < topo->atoms.size(); x++)
		{
			angles->push_back(AngleInfo());
		}

		for (unsigned int x = 0; x < topo->bonds.size(); x++)
		{
			tempBond1 = topo->bonds[x].atom1;
			tempBond2 = topo->bonds[x].atom2;
			if (tempBond2 != atomID)
			{
				(*angles)[tempBond1].addBond(tempBond2);
			}
			if (tempBond1 != atomID)
			{
				(*angles)[tempBond2].addBond(tempBond1);
			}
		}
		(*angles)[exclAtomID].setExclusionAtom();
		(*angles)[exclAtomID].setAngle(0);
		(*angles)[inAtomID].setInnerAtom();
		(*angles)[inAtomID].setVisited();
		(*angles)[outAtomID].setInnerAtom();
		(*angles)[outAtomID].setVisited();

		nodeStack.addElement(atomID);
		(*angles)[atomID].setAngle(rotAngle);
		curElement = atomID;
		index = 0;
		while (nodeStack.getNumElements() > 0)
		{
			done = false;
			(*angles)[curElement].setVisited();
			while (!done)
			{
				if (index < ((*angles)[curElement].numBonds()))
				{
					if ((*angles)[(*angles)[curElement].getBond(index)].isInnerAtom())
					{
						if (curElement != atomID)
						{
							nodeStack.addElement((*angles)[curElement].getBond(index));
							set_angles(&nodeStack, angles, true, rotAngle);
							nodeStack.popElement();
						}
						index++;
					}
					else if ((*angles)[(*angles)[curElement].getBond(index)].isExclusionAtom())
					{
						nodeStack.addElement((*angles)[curElement].getBond(index));
						set_angles(&nodeStack, angles, false, rotAngle);
						nodeStack.popElement();
						index++;
					}
					else if ((*angles)[(*angles)[curElement].getBond(index)].isVisited())
					{
						if ((*angles)[curElement].getBond(index) != nodeStack.getElement(nodeStack.getNumElements() - 2))
						{
							nodeStack.addElement((*angles)[curElement].getBond(index));
							set_angles(&nodeStack, angles, false, rotAngle);
							nodeStack.popElement();
						}
						index++;
					}
					else
					{
						indexStack.addElement(index);
						(*angles)[(*angles)[curElement].getBond(index)].setPointer(curElement);
						curElement = (*angles)[curElement].getBond(index);
						nodeStack.addElement(curElement);
						done = true;
						index = 0;
					}
				}
				if ((index > ((*angles)[curElement].numBonds() - 1)) || ((*angles)[curElement].numBonds() == 0))
				{
					index = 0;
					nodeStack.popElement();
					if (indexStack.getNumElements() > 0)
					{
						index = indexStack.popElement() + 1;
					}
					curElement = nodeStack.getElement(nodeStack.getNumElements() - 1);
					done = true;
				}
			}
		}

		done = false;
		nodeStack.reset(false);
		index = 0;
		rotAngle = -1.0;
		while (!done)
		{
			if ((*angles)[index].getAngleType() == AngleInfo::ANGLE_POINTER)
			{
				nodeStack.addElement(index);
				curElement = (*angles)[index].getPointer();
				while (nodeStack.getNumElements() > 0)
				{
					if ((*angles)[curElement].getAngleType() == AngleInfo::ANGLE_VALUE)
					{
						rotAngle = (*angles)[curElement].getAngle();
						curElement = nodeStack.popElement();
					}
					else if (((*angles)[curElement].getAngleType() == AngleInfo::ANGLE_POINTER) && (rotAngle != -1.0))
					{
						(*angles)[curElement].setAngle(rotAngle);
						curElement = nodeStack.popElement();
					}
					else
					{
						nodeStack.addElement(curElement);
						curElement = (*angles)[curElement].getPointer();
					}
				}
				(*angles)[curElement].setAngle(rotAngle);
			}
			index++;
			if (index >= angles->size())
			{
				done = true;
			}
		}
	}

	//___________________________________________________________
	Real computePhiDihedral(const GenericTopology* topo, const Vector3DBlock* positions, int index)
	{
		int a1 = topo->dihedrals[index].atom1;
		int a2 = topo->dihedrals[index].atom2;
		int a3 = topo->dihedrals[index].atom3;
		int a4 = topo->dihedrals[index].atom4;

		//Vector3D r12 = (*positions)[a1] - (*positions)[a2];  // Vector from atom 1 to atom 2
		Vector3D r12 = topo->minimalDifference((*positions)[a2], (*positions)[a1]);
		//Vector3D r23 = (*positions)[a2] - (*positions)[a3];  // Vector from atom 2 to atom 3
		Vector3D r23 = topo->minimalDifference((*positions)[a3], (*positions)[a2]);
		//Vector3D r34 = (*positions)[a3] - (*positions)[a4];  // Vector from atom 3 to atom 4
		Vector3D r34 = topo->minimalDifference((*positions)[a4], (*positions)[a3]);

		// Cross product of r12 and r23, represents the plane shared by these two vectors
		Vector3D a = r12.cross(r23);
		// Cross product of r12 and r23, represents the plane shared by these two vectors
		Vector3D b = r23.cross(r34);

		Vector3D c = r23.cross(a);

		// Calculate phi.
		Real cosPhi = a.dot(b) / (a.norm() * b.norm());
		Real sinPhi = c.dot(b) / (c.norm() * b.norm());

		return -atan2(sinPhi, cosPhi);
	}

	//___________________________________________________________
	Real computePhiDihedralEnergy(const GenericTopology* topo, int index, Real phi)
	{
		Torsion currTorsion = topo->dihedrals[index];
		Real energy = 0.0;
		for (int i = 0; i < currTorsion.multiplicity; i++)
		{
			if (currTorsion.periodicity[i] > 0)
			{
				// Add energy
				energy += currTorsion.forceConstant[i] * (1.0 + cos(
					currTorsion.periodicity[i] * phi
					+ currTorsion.phaseShift[i]));
				//report << "force constant = " << currTorsion.forceConstant[i] << endr;
				//report << "periodicity    = " << currTorsion.periodicity[i] << endr;
				//report << "trial angle    = " << phi << endr;
				//report << "phase shift    = " << currTorsion.phaseShift[i] << endr;
			}

			else
			{
				Real diff = phi - currTorsion.phaseShift[i];
				if (diff < -Constant::M_PI)
					diff += 2 * Constant::M_PI;
				else if (diff > Constant::M_PI)
					diff -= 2 * Constant::M_PI;

				// Add energy
				energy += currTorsion.forceConstant[i] * diff * diff;
			}
		}

		return energy;
	}

	void general_rotation(unsigned int innerAtom1, unsigned int innerAtom2, Vector3DBlock* positions, vector<AngleInfo>* angles)
	{
		unsigned int numAtoms;
		Vector3D a;
		Real norm_a;
		Real d;
		Real cosTheta;
		Real sinTheta;
		Real temp[9];
		Real angle;

		numAtoms = angles->size();

		a.x = (*positions)[innerAtom2].x - (*positions)[innerAtom1].x;
		a.y = (*positions)[innerAtom2].y - (*positions)[innerAtom1].y;
		a.z = (*positions)[innerAtom2].z - (*positions)[innerAtom1].z;
		norm_a = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
		a.x = a.x / norm_a;
		a.y = a.y / norm_a;
		a.z = a.z / norm_a;

		d = sqrt(a.y * a.y + a.z * a.z);
		if (d == 0.0)
		{
			for (unsigned int c = 0; c < numAtoms; c++)
			{
				angle = (*angles)[c].getAngle();
				cosTheta = cos(angle);
				sinTheta = sin(angle);

				if (0.0 < angle)
				{
					temp[1] = (*positions)[c].y - (*positions)[innerAtom1].y;
					temp[2] = (*positions)[c].z - (*positions)[innerAtom1].z;

					(*positions)[c].y = temp[1] * cosTheta - temp[2] * sinTheta + (*positions)[innerAtom1].y;
					(*positions)[c].z = temp[1] * sinTheta + temp[2] * cosTheta + (*positions)[innerAtom1].z;
				}
			}
		}
		else
		{
			for (unsigned int c = 0; c < numAtoms; c++)
			{
				angle = (*angles)[c].getAngle();
				cosTheta = cos(angle);
				sinTheta = sin(angle);

				if (0.0 < angle)
				{
					cosTheta = cos(angle);
					sinTheta = sin(angle);
					temp[0] = (*positions)[c].x - (*positions)[innerAtom1].x;
					temp[1] = (*positions)[c].y - (*positions)[innerAtom1].y;
					temp[8] = (*positions)[c].z - (*positions)[innerAtom1].z;;
					temp[2] = d * temp[0] + -a.x * ((temp[1] * a.y + temp[8] * a.z) / d);
					temp[3] = (temp[1] * a.z - temp[8] * a.y) / d;
					temp[4] = (cosTheta * temp[2] - sinTheta * temp[3]);
					temp[5] = (a.x * temp[0] + (temp[1] * a.y + temp[8] * a.z));
					temp[6] = sinTheta * temp[2] + cosTheta * temp[3];
					temp[7] = -a.x * temp[4] + d * temp[5];

					(*positions)[c].x = d * temp[4] + a.x * temp[5] + (*positions)[innerAtom1].x;
					(*positions)[c].y = (a.z * (temp[6]) + a.y * (temp[7])) / d + (*positions)[innerAtom1].y;
					(*positions)[c].z = (-a.y * (temp[6]) + a.z * (temp[7])) / d + (*positions)[innerAtom1].z;
				}
			}
		}

		return;
	}


	void general_rotation(unsigned int innerAtom1, unsigned int innerAtom2, Vector3DBlock* positions, Vector3DBlock* velocities, vector<AngleInfo>* angles)
	{
		unsigned int numAtoms;
		Vector3D a;
		Real norm_a;
		Real d;
		Real cosTheta;
		Real sinTheta;
		Real temp[9];
		Real tempVelo[9];
		Real angle;

		numAtoms = angles->size();

		a.x = (*positions)[innerAtom2].x - (*positions)[innerAtom1].x;
		a.y = (*positions)[innerAtom2].y - (*positions)[innerAtom1].y;
		a.z = (*positions)[innerAtom2].z - (*positions)[innerAtom1].z;
		norm_a = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
		a.x = a.x / norm_a;
		a.y = a.y / norm_a;
		a.z = a.z / norm_a;

		d = sqrt(a.y * a.y + a.z * a.z);

		if (d == 0.0)
		{
			for (unsigned int c = 0; c < numAtoms; c++)
			{
				angle = (*angles)[c].getAngle();
				cosTheta = cos(angle);
				sinTheta = sin(angle);

				if (0.0 < angle)
				{
					temp[1] = (*positions)[c].y - (*positions)[innerAtom1].y;
					temp[2] = (*positions)[c].z - (*positions)[innerAtom1].z;

					tempVelo[1] = (*velocities)[c].y + temp[1];
					tempVelo[2] = (*velocities)[c].z + temp[2];

					(*positions)[c].y = temp[1] * cosTheta - temp[2] * sinTheta;
					(*positions)[c].z = temp[1] * sinTheta + temp[2] * cosTheta;

					(*velocities)[c].y = tempVelo[1] * cosTheta - tempVelo[2] * sinTheta - (*positions)[c].y;
					(*velocities)[c].z = tempVelo[1] * sinTheta + tempVelo[2] * cosTheta - (*positions)[c].z;

					(*positions)[c].y += (*positions)[innerAtom1].y;
					(*positions)[c].z += (*positions)[innerAtom1].z;
				}
			}
		}
		else
		{
			for (unsigned int c = 0; c < numAtoms; c++)
			{
				angle = (*angles)[c].getAngle();
				cosTheta = cos(angle);
				sinTheta = sin(angle);

				if (0.0 < angle)
				{
					cosTheta = cos(angle);
					sinTheta = sin(angle);
					temp[0] = (*positions)[c].x - (*positions)[innerAtom1].x;
					temp[1] = (*positions)[c].y - (*positions)[innerAtom1].y;
					temp[8] = (*positions)[c].z - (*positions)[innerAtom1].z;
					temp[2] = d * temp[0] + -a.x * ((temp[1] * a.y + temp[8] * a.z) / d);
					temp[3] = (temp[1] * a.z - temp[8] * a.y) / d;
					temp[4] = (cosTheta * temp[2] - sinTheta * temp[3]);
					temp[5] = (a.x * temp[0] + (temp[1] * a.y + temp[8] * a.z));
					temp[6] = sinTheta * temp[2] + cosTheta * temp[3];
					temp[7] = -a.x * temp[4] + d * temp[5];

					tempVelo[0] = (*velocities)[c].x + temp[0];
					tempVelo[1] = (*velocities)[c].y + temp[1];
					tempVelo[8] = (*velocities)[c].z + temp[8];
					tempVelo[2] = d * tempVelo[0] + -a.x * ((tempVelo[1] * a.y + tempVelo[8] * a.z) / d);
					tempVelo[3] = (tempVelo[1] * a.z - tempVelo[8] * a.y) / d;
					tempVelo[4] = (cosTheta * tempVelo[2] - sinTheta * tempVelo[3]);
					tempVelo[5] = (a.x * tempVelo[0] + (tempVelo[1] * a.y + tempVelo[8] * a.z));
					tempVelo[6] = sinTheta * tempVelo[2] + cosTheta * tempVelo[3];
					tempVelo[7] = -a.x * tempVelo[4] + d * tempVelo[5];

					(*positions)[c].x = d * temp[4] + a.x * temp[5];
					(*positions)[c].y = (a.z * (temp[6]) + a.y * (temp[7])) / d;
					(*positions)[c].z = (-a.y * (temp[6]) + a.z * (temp[7])) / d;

					(*velocities)[c].x = d * tempVelo[4] + a.x * tempVelo[5] - (*positions)[c].x;
					(*velocities)[c].y = (a.z * (tempVelo[6]) + a.y * (tempVelo[7])) / d - (*positions)[c].y;
					(*velocities)[c].z = (-a.y * (tempVelo[6]) + a.z * (tempVelo[7])) / d - (*positions)[c].z;

					(*positions)[c].x += (*positions)[innerAtom1].x;
					(*positions)[c].y += (*positions)[innerAtom1].y;
					(*positions)[c].z += (*positions)[innerAtom1].z;
				}
			}
		}

		return;
	}
}
