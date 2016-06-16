/* -*- c++ -*- */
#ifndef PAULTRAPEXTENDEDFORCE_H
#define PAULTRAPEXTENDEDFORCE_H

#include "ExtendedForce.h"
#include "PaulTrapExtendedForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"
#include "pmconstants.h"

namespace ProtoMol
{
	//_________________________________________________________________ PaulTrapExtendedForce
	template <class TBoundaryConditions>
	class PaulTrapExtendedForce : public ExtendedForce, private PaulTrapExtendedForceBase
	{
		/* @latexonly Coulomb crystal systems are defined by ions with an additional
		 magnetic trap~\cite{HaAv91,Horn00,Schi93,TKTT02}.
		 A common used trap is the Paul Trap attraction, which is given by
		
		 \begin{equation}
		 U = \sum^N_{i=1} \left( \frac {1}{2} m_i\omega_r^{2}(x_{i}^{2}+y_{i}^{2}) + \frac {1}{2} m_i\omega_z^{2}z_{i}^{2} \right)
		 \label{eq:PT}
		 \end{equation}
		
		 \begin{equation}
		 \omega_r(j) = \left(1+\frac{\alpha}{r_0^4}x^2y^2\right)\sqrt{(\omega_r(1)^2+0.5\omega_z(1)^2)\left(\frac{q_j}{q_1} \frac{m_1}{m_j}\right) ^2 - 0.5\omega_z(1)^2  \frac{q_j}{q_1}\frac{m_1}{m_j}}
		 \end{equation}
		 \begin{equation}
		 \omega_z(j) = \omega_z(1) \sqrt{\frac{q_j}{q_1}  \frac{m_1}{m_j}}
		 \end{equation}
		 $\omega_r(i)$ and $\omega_z(i)$ denote the real $\omega$'s for atom type $i$ in the Paul Trap equation 
		 (\ref{eq:PT}), where $\omega_r(1)$ and $\omega_z(1)$ are  the \ProtoMol\ input parameter defining the 
		 spherical isotropic oscillation frequencies for the ion with $q_1$ and $m_1$. $\alpha$ and $r_0$ define 
		 the non-spherical field behavior.
		
		 NOTE: In the equations above  $\omega_r(1)$ and $\omega_z(1)$ do not be identical for the ions with $q_1$ and $m_1$. @endlatexonly
		
		*/


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		PaulTrapExtendedForce();
		PaulTrapExtendedForce(Real omegaR, Real omegaZ, Real alpha, Real r0, Real u0, Vector3D k);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class PaulTrapExtendedForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Real cohesive(Real potentialEnergy)
		{
			return (potentialEnergy - myPaulUHom) * myPaulF;
		}

		Real epsilonEnergy(Real potentialEnergy)
		{
			return potentialEnergy * myPaulF;
		}

	private:
		void doEvaluate(const GenericTopology* topo,
		                const Vector3DBlock* positions,
		                const Vector3DBlock* velocities,
		                Vector3DBlock* forces,
		                ScalarStructure* energies, int from, int to);
		void initialize(const GenericTopology* Topo);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class SystemForce
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void evaluate(const GenericTopology* topo,
		                      const Vector3DBlock* positions,
		                      const Vector3DBlock* velocities,
		                      Vector3DBlock* forces,
		                      ScalarStructure* energies);

		virtual void parallelEvaluate(const GenericTopology* topo,
		                              const Vector3DBlock* positions,
		                              const Vector3DBlock* velocities,
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

		virtual void uncache()
		{
			myCached = false;
		};

	private:
		virtual Force* doMake(std::string&, std::vector<Value>) const;

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
			return 8;
		}

		virtual void getParameters(std::vector<Parameter>&) const;


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Real myOmegaR1;
		Real myOmegaZ1;
		Real myAlpha;
		Real myR0;
		Real myU0;
		Vector3D myK;

		// Short cuts/internal data
		Real myAlphaR04;
		Real myOmega;
		Real myA;
		Real myPaulF;
		Real myPaulUHom;
		bool myCached;
		std::vector<Real> myOmegaR;
		std::vector<Real> myOmegaZ;
	};

	//______________________________________________________________________ INLINES
	template <class TBoundaryConditions>
	PaulTrapExtendedForce<TBoundaryConditions>::PaulTrapExtendedForce(): ExtendedForce(), myOmegaR1(0.0), myOmegaZ1(0.0), myAlpha(0.0), myR0(0.0), myU0(0.0), myK(Vector3D(0.0, 0.0, 0.0)), myAlphaR04(0.0), myOmega(0.0), myA(0.0), myCached(false)
	{
	}

	template <class TBoundaryConditions>
	PaulTrapExtendedForce<TBoundaryConditions>::PaulTrapExtendedForce(Real omegaR, Real omegaZ, Real alpha, Real r0, Real u0, Vector3D k): ExtendedForce(), myOmegaR1(omegaR), myOmegaZ1(omegaZ), myU0(u0), myK(k), myOmega(0.0), myA(0.0), myCached(false)
	{
		if (alpha != 0.0 && r0 > 0.0)
		{
			myAlphaR04 = alpha / (r0 * r0 * r0 * r0);
			myAlpha = alpha;
			myR0 = r0;
		}
		else
		{
			myAlphaR04 = 0.0;
			myAlpha = 0.0;
			myR0 = 0;
		}
	}

	template <class TBoundaryConditions>
	inline void PaulTrapExtendedForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-omegaR", Value(myOmegaR1, ConstraintValueType::NotNegative())));
		parameters.push_back(Parameter("-omegaZ", Value(myOmegaZ1, ConstraintValueType::NotNegative())));
		parameters.push_back(Parameter("-alpha", Value(myAlpha), 0.0));
		parameters.push_back(Parameter("-r0", Value(myR0), 0.0));
		parameters.push_back(Parameter("-u0", Value(myU0), 0.0));
		parameters.push_back(Parameter("-k", Value(myK.x), 0.0));
		parameters.push_back(Parameter("", Value(myK.y), 0.0));
		parameters.push_back(Parameter("", Value(myK.z), 0.0));
	}

	template <class TBoundaryConditions>
	inline Force* PaulTrapExtendedForce<TBoundaryConditions>::doMake(std::string&, std::vector<Value> values) const
	{
		return new PaulTrapExtendedForce(values[0], values[1], values[2], values[3], values[4], Vector3D(values[5], values[6], values[7]));
	}

	template <class TBoundaryConditions>
	inline void PaulTrapExtendedForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
	                                                                 const Vector3DBlock* positions,
	                                                                 const Vector3DBlock* velocities,
	                                                                 Vector3DBlock* forces,
	                                                                 ScalarStructure* energies)
	{
		doEvaluate(topo, positions, velocities, forces, energies, 0, (int)positions->size());
	}

	template <class TBoundaryConditions>
	inline void PaulTrapExtendedForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
	                                                                         const Vector3DBlock* positions,
	                                                                         const Vector3DBlock* velocities,
	                                                                         Vector3DBlock* forces,
	                                                                         ScalarStructure* energies)
	{
		unsigned int n = positions->size();
		unsigned int count = numberOfBlocks(topo, positions);

		for (unsigned int i = 0; i < count; i++)
		{
			if (Parallel::next())
			{
				int to = (n * (i + 1)) / count;
				if (to > static_cast<int>(n))
					to = n;
				int from = (n * i) / count;
				doEvaluate(topo, positions, velocities, forces, energies, from, to);
			}
		}
	}

	template <class TBoundaryConditions>
	inline void PaulTrapExtendedForce<TBoundaryConditions>::doEvaluate(const GenericTopology* topo,
	                                                                   const Vector3DBlock* positions,
	                                                                   const Vector3DBlock* /*velocities*/,
	                                                                   Vector3DBlock* forces,
	                                                                   ScalarStructure* energies, int from, int to)
	{
		if (!myCached)
			initialize(topo);

		const TBoundaryConditions& boundary =
			(dynamic_cast<const SemiGenericTopology<TBoundaryConditions>&>(*topo)).boundaryConditions;
		const Real f1 = myU0 * Constant::BOLTZMANN * 2.0 * Constant::M_PI;
		const Real u1 = -myU0 * Constant::BOLTZMANN;

		Real e = 0.0;

		for (int i = from; i < to; i++)
		{
			const unsigned int type = topo->atoms[i].type;
			Vector3D pos(boundary.minimalPosition((*positions)[i]));
			Real x2 = pos.x * pos.x;
			Real y2 = pos.y * pos.y;
			Real z2 = pos.z * pos.z;
			Real x = myK.dot(pos) * 2.0 * Constant::M_PI;
			Real f2 = f1 * sin(x);
			Real omegaR2 = myOmegaR[type] * power<2>(1 + myAlphaR04 * x2 * y2);
			Real omegaZ2 = myOmegaZ[type];
			Vector3D f(omegaR2 * pos.x,
			           omegaR2 * pos.y,
			           omegaZ2 * pos.z);
			Vector3D a(f2 * myK.x,
			           f2 * myK.y,
			           f2 * myK.z);

			(*forces)[i] -= f + a;
			Real b = u1 * cos(x);

			e += 0.5 * (omegaR2 * (x2 + y2) + omegaZ2 * z2) + b;


			//    Report::report << Report::plain
			//           << (0.5*(omegaR2*(x2+y2)+omegaZ2*z2)+b)*ENERGY_TO_SI<<","
			//           <<b*ENERGY_TO_SI<<","<<x*180.0/Constant::M_PI<<Report::endr;
			//
			//    Report::report << Report::plain
			//           << f*FORCE_TO_SI<<","
			//           << a*FORCE_TO_SI<<Report::endr;
		}
		(*energies)[ScalarStructure::OTHER] += e;
		//Report::report << e*ENERGY_TO_SI<<Report::endr;
	}

	template <class TBoundaryConditions>
	inline void PaulTrapExtendedForce<TBoundaryConditions>::initialize(const GenericTopology* topo)
	{
		const unsigned int numberOfAtoms = topo->atoms.size();
		const unsigned int numberOfTypes = topo->atomTypes.size();

		unsigned int pos = 0;
		unsigned int neg = 0;
		for (unsigned int i = 0; i < numberOfTypes; i++)
		{
			if (topo->atomTypes[i].charge < 0.0)
				neg++;
			else if (topo->atomTypes[i].charge > 0.0)
				pos++;
		}
		if (pos != numberOfTypes && neg != numberOfTypes)
			Report::report << Report::error << "PaulTrap force expects only positive or negative charges!" << Report::endr;

		// Derived the real omegas from the first omega definition
		myOmegaR.resize(numberOfTypes);
		myOmegaZ.resize(numberOfTypes);
		for (unsigned int i = 0; i < numberOfTypes; i++)
		{
			Real f = topo->atomTypes[i].charge / topo->atomTypes[0].charge * topo->atomTypes[0].mass / topo->atomTypes[i].mass;
			myOmegaR[i] = (myOmegaR1 * myOmegaR1 + 0.5 * myOmegaZ1 * myOmegaZ1) * f * f - 0.5 * myOmegaZ1 * myOmegaZ1 * f;
			myOmegaZ[i] = myOmegaZ1 * myOmegaZ1 * f;
		}

		// Average omega
		myOmega = 0.0;
		for (unsigned int i = 0; i < numberOfAtoms; i++)
		{
			myOmega += (2.0 * myOmegaR[topo->atoms[i].type] + myOmegaZ[topo->atoms[i].type]) / 3.0;
		}
		myOmega = sqrt(myOmega / static_cast<Real>(numberOfAtoms));

		// Number of each ion type
		std::vector<int> n(numberOfTypes, 0);
		for (unsigned int i = 0; i < numberOfAtoms; i++)
		{
			n[topo->atoms[i].type]++;
		}

		// Average Q
		Real q = 0.0;
		Real m = 0.0;
		for (unsigned int i = 0; i < numberOfTypes; i++)
		{
			if (numberOfAtoms > 1)
				q += topo->atomTypes[i].charge * topo->atomTypes[i].charge * n[i] * (n[i] - 1.0) / 2.0;
			else
				q += topo->atomTypes[i].charge * topo->atomTypes[i].charge * n[i];
			m += topo->atomTypes[i].mass * n[i];
			for (unsigned int j = i + 1; j < numberOfTypes; j++)
			{
				q += topo->atomTypes[i].charge * topo->atomTypes[j].charge * n[i] * n[j];
			}
		}
		if (numberOfAtoms > 1)
			q = Constant::SQRTCOULOMBCONSTANT * sqrt(q / ((Real)numberOfAtoms * ((Real)numberOfAtoms - 1.0) / 2.0));
		else
			q = Constant::SQRTCOULOMBCONSTANT * sqrt(q);
		if (neg > 0)
			q = -q;

		// Wigner Seitz radius
		Real k = m / (Real)numberOfAtoms * myOmega * myOmega * 1.0e7 / 4184.0;
		myA = pow(q * q / k, 1.0 / 3.0);
		Real r = pow((Real)numberOfAtoms * q * q / k, 1.0 / 3.0);
		myPaulUHom = 9.0 / 10.0 * pow((Real)numberOfAtoms, 5.0 / 3.0) * q * q / myA;
		myPaulF = 1.0 / ((Real)numberOfAtoms * q * q / myA);

		// Output
		Report::report.precision(13);
		Report::report << Report::plain << "Paul Trap F : w           = " << myOmega << "[fs-1]" << Report::endr;
		for (unsigned int i = 0; i < numberOfTypes; i++)
		{
			Report::report << Report::plain << "            : w_" << i << "         = ("
				<< sqrt((2.0 * myOmegaR[topo->atoms[i].type] + myOmegaZ[topo->atoms[i].type]) / 3.0) << ","
				<< sqrt(myOmegaR[i]) << ","
				<< sqrt(myOmegaZ[i]) << ")[fs-1]" << Report::endr;
		}
		Report::report << Report::plain
			<< "            : m           = " << m / (Real)numberOfAtoms << "[amu]" << std::endl
			<< "            : k           = " << k << "[]" << std::endl
			<< "            : q           = " << q / Constant::SQRTCOULOMBCONSTANT << "[e]" << std::endl
			<< "            : a           = " << myA * 1e-4 << "[um]" << std::endl
			<< "            : R           = " << r * 1e-4 << "[um]" << std::endl
			<< "            : Uhom        = " << myPaulUHom << "[kcal/mol]" << std::endl
			<< "            :             = " << myPaulUHom * myPaulF << "[Nq^2/a]" << std::endl
			<< "            : 1[K]        = " << (1.0 * Constant::BOLTZMANN * 3.0 * myA) / (2.0 * q * q) << std::endl;


		// Scaling factor
		for (unsigned int i = 0; i < numberOfTypes; i++)
		{
			Real c = topo->atomTypes[i].mass * 1e-3 * Constant::SI::KCAL * Constant::SI::TIME_FS * Constant::SI::TIME_FS / Constant::SI::LENGTH_AA / Constant::SI::LENGTH_AA;
			myOmegaR[i] *= c;
			myOmegaZ[i] *= c;
		}


		myCached = true;
	}

	template <class TBoundaryConditions>
	inline unsigned int PaulTrapExtendedForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology*,
	                                                                               const Vector3DBlock* pos)
	{
		return std::min(Parallel::getAvailableNum(), static_cast<int>(pos->size()));
	}
}
#endif /* PAULTRAPEXTENDEDFORCE_H */
