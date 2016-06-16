/* -*- c++ -*- */
#ifndef TOPOLOGYUTILITIES_H
#define TOPOLOGYUTILITIES_H

#include "Vector3DBlock.h"
#include "Vector3D.h"
#include "Matrix3by3.h"
#include "Real.h"
#include "Torsion.h"
#include "Stack.h"
#include "AngleInfo.h"
#include "Bond.h"

#include <string>
#include <vector>
#include <set>


namespace ProtoMol
{
	enum
	{
		X_AXIS,
		Y_AXIS,
		Z_AXIS
	};

	class GenericTopology;
	class ScalarStructure;

	//___________________________________________________________randomVelocity
	/// Gaussian random distributed velocities
	void randomVelocity(Real temperature,
	                    const GenericTopology* topology,
	                    Vector3DBlock* velocities,
	                    unsigned int seed = 1234);

	//___________________________________________________________randomVelocity
	/// Gaussian random distributed velocities rescaled to the given interval with optional re-movements of linear and/or angular momentum
	void randomVelocity(Real temperatureFrom,
	                    Real temperatureTo,
	                    const GenericTopology* topology,
	                    const Vector3DBlock* positions,
	                    Vector3DBlock* velocities,
	                    bool removeLinear = false,
	                    bool removeAngular = false,
	                    unsigned int seed = 1234);

	//___________________________________________________________getAtomsBondedtoDihedral
	/// this function gets all the atoms bonded to ONE side of the dihedral

	void getAtomsBondedtoDihedral(const GenericTopology* topology,
	                              std::set<int>* atomSet,
	                              const int atomID,
	                              const int inAtomID,
	                              const int outAtomID,
	                              const int exclAtomID);

	//___________________________________________________________rotateDihedral
	/// this function rotates all the atoms bonded to ONE side of the dihedral

	void rotateDihedral(const GenericTopology* topology,
	                    Vector3DBlock* positions,
	                    const int dihedralID, Real angle);

	void rotateDihedral(const GenericTopology* topology,
	                    Vector3DBlock* positions,
	                    Vector3DBlock* velocities,
	                    const int dihedralID, Real angle);

	//___________________________________________________________kineticEnergy
	Real kineticEnergy(const GenericTopology* topology,
	                   const Vector3DBlock* velocities);

	//___________________________________________________________molecularKineticEnergy
	Real molecularKineticEnergy(const GenericTopology* topology,
	                            const Vector3DBlock* velocitiies);

	//___________________________________________________________kineticEnergyForAtomType
	enum waterOption
	{
		IGNORE_WATER,
		ONLY_WATER,
		ALL
	};

	Real kineticEnergyForAtomType(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int atomType,
	                              waterOption option);
	Real kineticEnergyForAtomType(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int atomType,
	                              waterOption option,
	                              int& atomCount);

	//___________________________________________________________kineticEnergyForWater
	Real kineticEnergyForWater(const GenericTopology* topology,
	                           const Vector3DBlock* velocities);
	Real kineticEnergyForWater(const GenericTopology* topology,
	                           const Vector3DBlock* velocities,
	                           int& waterCount);

	//___________________________________________________________kineticEnergyForNonWater
	Real kineticEnergyForNonWater(const GenericTopology* topology,
	                              const Vector3DBlock* velocities);
	Real kineticEnergyForNonWater(const GenericTopology* topology,
	                              const Vector3DBlock* velocities,
	                              int& nonWaterCount);

	//___________________________________________________________temperature
	Real temperature(const GenericTopology* topology,
	                 const Vector3DBlock* velocities);

	//___________________________________________________________temperature
	Real temperature(Real kineticEnergy,
	                 unsigned int degreesOfFreedom);

	//___________________________________________________________temperatureForAtomType
	Real temperatureForAtomType(const GenericTopology* topology,
	                            const Vector3DBlock* velocities,
	                            int atomType,
	                            waterOption option);

	//___________________________________________________________temperatureForWater
	Real temperatureForWater(const GenericTopology* topology,
	                         const Vector3DBlock* velocities);

	//___________________________________________________________temperatureForNonWater
	Real temperatureForNonWater(const GenericTopology* topology,
	                            const Vector3DBlock* velocities);

	//___________________________________________________________getNonWaterAtoms
	int getNonWaterAtoms(const GenericTopology* topology);

	//___________________________________________________________atomTypeToSymbolName
	std::string atomTypeToSymbolName(const std::string& type);

	//___________________________________________________________linearMomentum
	Vector3D linearMomentum(const Vector3DBlock* velocities,
	                        const GenericTopology* topo);

	//___________________________________________________________linearMomentumSolute
	Vector3D linearMomentumSolute(const Vector3DBlock* velocities,
	                              const GenericTopology* topo);

	//___________________________________________________________centerOfMass
	Vector3D centerOfMass(const Vector3DBlock* positions,
	                      const GenericTopology* topo);

	//___________________________________________________________angularMomentum
	Vector3D angularMomentum(const Vector3DBlock* positions,
	                         const Vector3DBlock* velocities,
	                         const GenericTopology* topo);

	//___________________________________________________________angularMomentum
	Vector3D angularMomentum(const Vector3DBlock* positions,
	                         const Vector3DBlock* velocities,
	                         const GenericTopology* topo,
	                         const Vector3D& centerOfMass);

	//___________________________________________________________angularMomentumSolute
	Vector3D angularMomentumSolute(const Vector3DBlock* positions,
	                               const Vector3DBlock* velocities,
	                               const GenericTopology* topo,
	                               const Vector3D& centerOfMass);

	//___________________________________________________________inertiaMomentum
	Matrix3by3 inertiaMomentum(const Vector3DBlock* positions,
	                           const GenericTopology* topo,
	                           const Vector3D& centerOfMass);

	//___________________________________________________________inertiaMomentumSolute
	Matrix3by3 inertiaMomentumSolute(const Vector3DBlock* positions,
	                                 const GenericTopology* topo,
	                                 const Vector3D& centerOfMass);

	//___________________________________________________________removeLinearMomentum
	Vector3D removeLinearMomentum(Vector3DBlock* velocities,
	                              const GenericTopology* topo);

	//___________________________________________________________removeAngularMomentum
	Vector3D removeAngularMomentum(const Vector3DBlock* positions,
	                               Vector3DBlock* velocities,
	                               const GenericTopology* topo);

	//___________________________________________________________velocityVirial
	ScalarStructure velocityVirial(const GenericTopology* topology,
	                               const Vector3DBlock* velocities);

	//___________________________________________________________addVelocityVirial
	void addVelocityVirial(ScalarStructure* energies,
	                       const GenericTopology* topology,
	                       const Vector3DBlock* velocities);

	//___________________________________________________________computePressure
	Real computePressure(const GenericTopology* topology,
	                     const Vector3DBlock* positions,
	                     const Vector3DBlock* velocities,
	                     const ScalarStructure* energies);

	//___________________________________________________________computePressure
	Real computePressure(const ScalarStructure* energies,
	                     Real volume,
	                     Real kineticEnergy);

	//___________________________________________________________computeMolecularPressure
	Real computeMolecularPressure(const ScalarStructure* energies,
	                              Real volume,
	                              Real kineticEnergy);

	//___________________________________________________________molecularMomentum
	Vector3D molecularMomentum(const std::vector<int>&,
	                           const Vector3DBlock*,
	                           const GenericTopology*);

	//___________________________________________________________molecularCenterOfMass
	Vector3D molecularCenterOfMass(const std::vector<int>&,
	                               const Vector3DBlock*,
	                               const GenericTopology*);

	//___________________________________________________________buildMolecularCenterOfMass
	void buildMolecularCenterOfMass(const Vector3DBlock* positions, GenericTopology* topo);

	//___________________________________________________________buildMolecularMomentum
	void buildMolecularMomentum(const Vector3DBlock* velocities, GenericTopology* topo);


	//___________________________________________________________buildRattleShakeBondConstraintList
	void buildRattleShakeBondConstraintList(GenericTopology* topo, std::vector<Bond::Constraint>& bondConstraints);

	void build_angle_list(const GenericTopology* topo, const unsigned int atomID, const unsigned int inAtomID, const unsigned int outAtomID, const unsigned int exclAtomID, Real rotAngle, std::vector<AngleInfo>* angles);
	void set_angles(Stack<unsigned int>* nodeStack, std::vector<AngleInfo>* angles, bool lastIsInnerAtom, Real wholeAngle);
	void general_rotation(unsigned int innerAtom1, unsigned int innerAtom2, Vector3DBlock* positions, std::vector<AngleInfo>* angles);
	void general_rotation(unsigned int innerAtom1, unsigned int innerAtom2, Vector3DBlock* positions, Vector3DBlock* velocities, std::vector<AngleInfo>* angles);
	Real computePhiDihedral(const GenericTopology* topo, const Vector3DBlock* positions, int index);
	Real computePhiDihedralEnergy(const GenericTopology* topo, int index, Real phi);
}

#endif // TOPOLOGYUTILITIES_H
