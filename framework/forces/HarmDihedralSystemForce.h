/* -*- c++ -*- */
#ifndef HARMDIHEDRALSYSTEMFORCE_H
#define HARMDIHEDRALSYSTEMFORCE_H

#include "MTorsionSystemForce.h"
#include "HarmDihedralSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"
#include "topologyutilities.h"

//for debugging
#include "Report.h"

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ HarmDihedralSystemForce

	template <class TBoundaryConditions>
	class HarmDihedralSystemForce: public MTorsionSystemForce<TBoundaryConditions>, private HarmDihedralSystemForceBase
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		// Default constructor
		HarmDihedralSystemForce(): k(1), myDihedral(0), myDihedralReference(0.0),
		                           computeOthers(false)
		{
		}


		// Constructor with parameters
		HarmDihedralSystemForce(Real kbias, int dihedral, Real dihedralReference, bool other) : k(kbias), myDihedral(dihedral),
		                                                                                        myDihedralReference(dtor(dihedralReference)),
		                                                                                        computeOthers(other)
		{
		}


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class HarmDihedralSystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void harmCalcTorsion(const GenericTopology* topo,
		                     const TBoundaryConditions& boundary,
		                     const Torsion& currentTorsion,
		                     const Vector3DBlock* positions,
		                     Vector3DBlock* forces,
		                     Real& energy,
		                     ScalarStructure* energies);


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class SystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* positions,
		                      Vector3DBlock* forces,
		                      ScalarStructure* energies);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              Vector3DBlock* forces,
		                              ScalarStructure* energies);

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
	private:
		virtual Force* doMake(std::string&, std::vector<Value> values) const
		{
			return (new HarmDihedralSystemForce(values[0], values[1], values[2], values[3]));
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual unsigned int getParameterSize() const
		{
			return 4;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const
		{
			parameters.push_back(Parameter("-kbias", Value(k, ConstraintValueType::NotNegative()), Text("potential bias constant")));
			parameters.push_back(Parameter("-dihedral", Value(myDihedral, ConstraintValueType::NotNegative())));
			parameters.push_back(Parameter("-angle", Value(rtod(myDihedralReference)), Text("reference angle -180 to 180")));
			parameters.push_back(Parameter("-others", Value(computeOthers, ConstraintValueType::NoConstraints())));
		}

	private:
		virtual void doSetParameters(std::string&, std::vector<Value> values)
		{
			k = values[0];
			myDihedral = values[1];
			myDihedralReference = dtor((Real)values[2]);
			computeOthers = values[3];
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		// The harmonic potential biasing constant
		Real k;
		// The dihedral reference angle ID
		int myDihedral;
		// The dihedral reference angle value
		Real myDihedralReference;
		bool computeOthers;
	};

	//______________________________________________________________________ INLINES


	template <class TBoundaryConditions>
	inline void HarmDihedralSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
	                                                                   const Vector3DBlock* positions,
	                                                                   Vector3DBlock* forces,
	                                                                   ScalarStructure* energies)
	{
		const TBoundaryConditions& boundary =
			(dynamic_cast<const SemiGenericTopology<TBoundaryConditions>&>(*topo)).boundaryConditions;

		// Examine each dihedral
		for (unsigned int i = 0; i < topo->dihedrals.size(); i++)
		{
			if (computeOthers)
			{
				calcTorsion(boundary, topo->dihedrals[i], positions, forces, (*energies)[ScalarStructure::DIHEDRAL], energies);
			}
			if (static_cast<int>(i) == myDihedral)
			{
				harmCalcTorsion(topo, boundary, topo->dihedrals[i], positions, forces, (*energies)[ScalarStructure::DIHEDRAL], energies);
			}
		}
	}

	template <class TBoundaryConditions>
	inline void HarmDihedralSystemForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
	                                                                           const Vector3DBlock* positions,
	                                                                           Vector3DBlock* forces,
	                                                                           ScalarStructure* energies)
	{
		const TBoundaryConditions& boundary =
			(dynamic_cast<const SemiGenericTopology<TBoundaryConditions>&>(*topo)).boundaryConditions;

		unsigned int n = topo->dihedrals.size();
		unsigned int count = numberOfBlocks(topo, positions);

		for (unsigned int i = 0; i < count; i++)
		{
			if (Parallel::next())
			{
				int to = (n * (i + 1)) / count;
				if (to > static_cast<int>(n))
					to = n;
				int from = (n * i) / count;
				for (int j = from; j < to; j++)
				{
					if (j != myDihedral)
					{
						calcTorsion(boundary, topo->dihedrals[j], positions, forces, (*energies)[ScalarStructure::DIHEDRAL], energies);
					}
					else
					{
						harmCalcTorsion(topo, boundary, topo->dihedrals[j], positions, forces, (*energies)[ScalarStructure::DIHEDRAL], energies);
					}
				}
			}
		}
	}


	template <class TBoundaryConditions>
	inline void HarmDihedralSystemForce<TBoundaryConditions>::harmCalcTorsion(const GenericTopology* topo,
	                                                                          const TBoundaryConditions& boundary,
	                                                                          const Torsion& /*currTorsion*/,
	                                                                          const Vector3DBlock* positions,
	                                                                          Vector3DBlock* forces,
	                                                                          Real& energy,
	                                                                          ScalarStructure* energies)
	{
		// Calculate the energy and forces due to the harmonic potential
		// -------------------------------------------------------------
		// The dihedral angle
		Real dihedralAngle = computePhiDihedral(topo, positions, myDihedral);

		// -------------------------------------------------------------
		// The dihedral potential
		//Real k = 3.0; // kcal / (mol*K)

		Real diff = (dihedralAngle - myDihedralReference);
		if (diff < -Constant::M_PI)
			diff += 2 * Constant::M_PI;
		else if (diff > Constant::M_PI)
			diff -= 2 * Constant::M_PI;
		Real V = k * diff * diff;

		/* Replaced by PRB 5/24/05
		Real V = k * (dihedralAngle - myDihedralReference) * (dihedralAngle - myDihedralReference);
		report << debug(1) << "actual dihedral angle: " << dihedralAngle
		       << "   target dihedral angle: " << myDihedralReference <<  endr;
		*/
		//report << debug(1) << "diff: " << diff << "   potential bias V: " << V << endr;

		// -------------------------------------------------------------
		// The atoms in the dihedral
		int ai = topo->dihedrals[myDihedral].atom1;
		int aj = topo->dihedrals[myDihedral].atom2;
		int ak = topo->dihedrals[myDihedral].atom3;
		int al = topo->dihedrals[myDihedral].atom4;

		// -------------------------------------------------------------
		// The vector coordinates of the atoms
		Vector3D ri((*positions)[ai]);
		Vector3D rj((*positions)[aj]);
		Vector3D rk((*positions)[ak]);
		Vector3D rl((*positions)[al]);
		// -------------------------------------------------------------
		// Vectors between atoms (rij = ri - rj)
		Vector3D rij = boundary.minimalDifference(rj, ri);
		Vector3D rkj = boundary.minimalDifference(rj, rk);
		Vector3D rkl = boundary.minimalDifference(rl, rk);
		// -------------------------------------------------------------
		// Normals
		Vector3D m = rij.cross(rkj);
		Vector3D n = rkj.cross(rkl);
		// -------------------------------------------------------------
		// The derivate of V with respect to the dihedral angle
		Real dVdPhi = 2 * (k) * diff;
		/* Replaced by PRB 5/24/05
		Real dVdPhi = (k) * (dihedralAngle - myDihedralReference);
		*/
		// -------------------------------------------------------------
		// Miscellaneous quantities needed to compute the forces
		Real rkj_norm = rkj.norm();
		Real rkj_normsq = rkj.normSquared();
		Real m_normsq = m.normSquared();
		Real n_normsq = n.normSquared();
		Real rij_dot_rkj = rij.dot(rkj);
		Real rkl_dot_rkj = rkl.dot(rkj);
		// -------------------------------------------------------------
		// Forces on the atoms
		// Atom i
		Vector3D fi = m * (-dVdPhi * rkj_norm / m_normsq);
		// Atom l
		Vector3D fl = n * (dVdPhi * rkj_norm / n_normsq);
		// Atom j
		Vector3D fj = fi * (-1 + rij_dot_rkj / rkj_normsq) - fl * (rkl_dot_rkj / rkj_normsq);
		// Atom k
		Vector3D fk = - (fi + fj + fl);


		// *************************************************************
		// -------------------------------------------------------------
		// Update the dihedral system energy
		energy += V;
		// Update the atom forces
		(*forces)[ai] += fi;
		(*forces)[aj] += fj;
		(*forces)[ak] += fk;
		(*forces)[al] += fl;
		// -------------------------------------------------------------
		// Compute the virial energy
		if (energies->virial())
		{
			Vector3D f1 = fi;
			Vector3D f2 = fi + fj;
			Vector3D f3 = -fl;
			Vector3D r12 = rij;
			Vector3D r23 = -rkj;
			Vector3D r34 = rkl;

			Real xy = f1.x * r12.y + f2.x * r23.y + f3.x * r34.y;
			Real xz = f1.x * r12.z + f2.x * r23.z + f3.x * r34.z;
			Real yz = f1.y * r12.z + f2.y * r23.z + f3.y * r34.z;

			(*energies)[ScalarStructure::VIRIALXX] += f1.x * r12.x + f2.x * r23.x + f3.x * r34.x;
			(*energies)[ScalarStructure::VIRIALXY] += xy;
			(*energies)[ScalarStructure::VIRIALXZ] += xz;
			(*energies)[ScalarStructure::VIRIALYX] += xy;
			(*energies)[ScalarStructure::VIRIALYY] += f1.y * r12.y + f2.y * r23.y + f3.y * r34.y;
			(*energies)[ScalarStructure::VIRIALYZ] += yz;
			(*energies)[ScalarStructure::VIRIALZX] += xz;
			(*energies)[ScalarStructure::VIRIALZY] += yz;
			(*energies)[ScalarStructure::VIRIALZZ] += f1.z * r12.z + f2.z * r23.z + f3.z * r34.z;
		}
	}

	template <class TBoundaryConditions>
	inline unsigned int HarmDihedralSystemForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology* topo,
	                                                                                 const Vector3DBlock*)
	{
		return std::min(Parallel::getAvailableNum(), static_cast<int>(topo->dihedrals.size()));
	}
}
#endif /* HARMDIHEDRALSYSTEMFORCE_H */
