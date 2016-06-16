#include "OutputCache.h"
#include "Output.h"
#include "PDB.h"
#include "Configuration.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "Integrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "systemutilities.h"
#include "topologyutilities.h"
#include "pmconstants.h"
#include <limits>
using std::vector;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//________________________________________________________ OutputCache
	OutputCache::OutputCache(): myConfig(NULL),
	                            myTopology(NULL),
	                            myIntegrator(NULL),
	                            myEnergies(NULL),
	                            myPositions(NULL),
	                            myVelocities(NULL),
	                            myInitialPositions(new Vector3DBlock()),
	                            myMinimalPositions(new Vector3DBlock()),
	                            myPositionsNoWater(new Vector3DBlock()),
	                            myDihedralPhi(Constant::REAL_NAN),
	                            myDihedralPhis(new vector<Real>()),
	                            myBrentMaxima(new vector<vector<Real>>())
	{
		uncache();
	}

	OutputCache::~OutputCache()
	{
		if (myInitialPositions != NULL)
			delete myInitialPositions;
		if (myMinimalPositions != NULL)
			delete myMinimalPositions;
		if (myDihedralPhis != NULL)
			delete myDihedralPhis;
		if (myBrentMaxima != NULL)
			delete myBrentMaxima;
	}

	void OutputCache::initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
	                             const Vector3DBlock* pos, const Vector3DBlock* vel, const ScalarStructure* energies)
	{
		myConfig = config;
		myTopology = topo;
		myIntegrator = integrator;
		myEnergies = energies;
		myPositions = pos;
		myVelocities = vel;
		(*myInitialPositions) = (*pos);
	}

	Real OutputCache::totalEnergy() const
	{
		return (potentialEnergy() + kineticEnergy());
	}

	Real OutputCache::potentialEnergy() const
	{
		if (!myCachedPE)
		{
			myPE = myEnergies->potentialEnergy();
			myCachedPE = true;
		}
		return myPE;
	}

	Real OutputCache::kineticEnergy() const
	{
		if (!myCachedKE)
		{
			myKE = ProtoMol::kineticEnergy(myTopology, myVelocities);
			myT = ProtoMol::temperature(myKE, myTopology->degreesOfFreedom);
			myCachedKE = true;
		}
		return myKE;
	}

	Real OutputCache::temperature() const
	{
		if (!myCachedKE)
		{
			myKE = ProtoMol::kineticEnergy(myTopology, myVelocities);
			myT = ProtoMol::temperature(myKE, myTopology->degreesOfFreedom);
			myCachedKE = true;
		}
		return myT;
	}

	Real OutputCache::temperatureForWater() const
	{
		if (!myCachedWaterT)
		{
			myWaterT = ProtoMol::temperatureForWater(myTopology, myVelocities);
			myCachedWaterT = true;
		}
		return myWaterT;
	}

	Real OutputCache::temperatureForNonWater() const
	{
		if (!myCachedNonWaterT)
		{
			myNonWaterT = ProtoMol::temperatureForNonWater(myTopology, myVelocities);
			myCachedNonWaterT = true;
		}
		return myNonWaterT;
	}

	Real OutputCache::molecularTemperature() const
	{
		if (!myCachedMolT)
		{
			myMolT = ProtoMol::temperature(molecularKineticEnergy(), 3 * myTopology->molecules.size());
			myCachedMolT = true;
		}
		return myMolT;
	}

	Real OutputCache::molecularKineticEnergy() const
	{
		if (!myCachedMolKE)
		{
			myMolKE = ProtoMol::molecularKineticEnergy(myTopology, myVelocities);
			myCachedMolKE = true;
		}
		return myMolKE;
	}


	Real OutputCache::pressure() const
	{
		if (!myCachedP)
		{
			if (!myEnergies->virial())
				myP = 0.0;
			else if (volume() > 0.0)
				myP = ProtoMol::computePressure(myEnergies, volume(), kineticEnergy());
			else
				myP = Constant::REAL_INFINITY;
			myCachedP = true;
		}
		return myP;
	}

	Real OutputCache::molecularPressure() const
	{
		if (!myCachedMolP)
		{
			if (!myEnergies->molecularVirial())
				myMolP = 0.0;
			else if (volume() > 0.0)
				myMolP = ProtoMol::computeMolecularPressure(myEnergies, volume(), molecularKineticEnergy());
			else
				myMolP = Constant::REAL_INFINITY;
			myCachedMolP = true;
		}
		return myMolP;
	}

	Real OutputCache::volume() const
	{
		if (!myCachedV)
		{
			myV = myTopology->getVolume(*myPositions);
			myCachedV = true;
		}
		return myV;
	}

	Vector3D OutputCache::linearMomentum() const
	{
		if (!myCachedLinearMomentum)
		{
			myLinearMomentum = ProtoMol::linearMomentum(myVelocities, myTopology);
			myCachedLinearMomentum = true;
		}
		return myLinearMomentum;
	}

	Vector3D OutputCache::angularMomentum() const
	{
		if (!myCachedAngularMomentum)
		{
			myAngularMomentum = ProtoMol::angularMomentum(myPositions, myVelocities, myTopology, OutputCache::centerOfMass());
			myCachedAngularMomentum = true;
		}
		return myAngularMomentum;
	}

	Vector3D OutputCache::centerOfMass() const
	{
		if (!myCachedCenterOfMass)
		{
			myCenterOfMass = ProtoMol::centerOfMass(myPositions, myTopology);
			myCachedCenterOfMass = true;
		}
		return myCenterOfMass;
	}

	Real OutputCache::diffusion() const
	{
		if (!myCachedDiffusion)
		{
			myDiffusion = 0.0;
			unsigned int numberOfAtoms = myPositions->size();
			for (unsigned int i = 0; i < numberOfAtoms; i++)
				myDiffusion += ((*myPositions)[i] - (*myInitialPositions)[i]).normSquared();
			myDiffusion /= (6.0 * numberOfAtoms);
			myCachedDiffusion = true;
		}
		return myDiffusion;
	}

	Real OutputCache::density() const
	{
		if (!myCachedDensity)
		{
			myDensity = (volume() > 0.0 ? (mass() / volume() * Constant::SI::AMU * power<3>(Constant::SI::LENGTH_AA) * 1e-3) : Constant::REAL_NAN);
			myCachedDensity = true;
		}
		return myDensity;
	}

	Real OutputCache::mass() const
	{
		if (!myCachedMass)
		{
			myMass = 0.0;
			unsigned int numberOfAtoms = myPositions->size();
			for (unsigned int i = 0; i < numberOfAtoms; i++)
				myMass += myTopology->atoms[i].scaledMass;
			myCachedMass = true;
		}
		return myMass;
	}

	Real OutputCache::time() const
	{
		return myTopology->time;
	}

	const Vector3DBlock* OutputCache::minimalPositions() const
	{
		if (!myCachedMinimalPositions)
		{
			*myMinimalPositions = *myPositions;
			(const_cast<GenericTopology*>(myTopology))->minimalImage(*myMinimalPositions);
		}
		myCachedMinimalPositions = true;
		return myMinimalPositions;
	}

	const Vector3DBlock* OutputCache::PositionsNoWater() const
	{
		//cout << "Before clear myPositionsNoWater.size() = " << (*myPositionsNoWater).size() << endl;
		(*myPositionsNoWater).clear();
		//cout << "After clear myPositionsNoWater.size() = " << (*myPositionsNoWater).size() << endl;
		int atomsNoWater = ProtoMol::getNonWaterAtoms(myTopology);
		//cout << "Number of atomsNoWater = " << atomsNoWater << endl;
		for (unsigned int i = 0; i < atomsNoWater; ++i)
		{
			(*myPositionsNoWater).push_back((*myPositions)[i]);
		}
		//cout << "After push_back loop myPositionsNoWater.size() = " << (*myPositionsNoWater).size() << endl;
		return myPositionsNoWater;
	}


	Real OutputCache::dihedralPhi(int index) const
	{
		if (index < 0 || index >= static_cast<int>(myTopology->dihedrals.size()))
			index = -1;

		if (index < 0)
		{
			myDihedralPhi = Constant::REAL_NAN;
			myCachedDihedralPhi = index;
			return myDihedralPhi;
		}

		/*
		if(myCachedDihedralPhi < 0){
		  myDihedralPhi  = computePhiDihedral(myTopology,myPositions,index);
		}
		myCachedDihedralPhi = index;
		*/
		myDihedralPhi = computePhiDihedral(myTopology, myPositions, index);
		return myDihedralPhi;
	}

	vector<Real> OutputCache::dihedralPhis(vector<int> dihedralset) const
	{
		if (!myCachedDihedralPhis)
		{
			myDihedralPhis->resize(dihedralset.size());
			for (unsigned int i = 0; i < dihedralset.size(); ++i)
			{
				(*myDihedralPhis)[i] = computePhiDihedral(myTopology, myPositions, dihedralset[i]);
			}
			myCachedDihedralPhis = false; // different functions require different dihedralset
		}
		return (*myDihedralPhis);
	}

	// Brent's Maxima function goes here to use topology for dihedral well calculation
	vector<vector<Real>> OutputCache::brentMaxima(vector<int> dihedralset, bool max) const
	{
		if (!myCachedBrentMaxima)
		{
			//The Brent algorithm gets the maxima if maxmin = -1
			int maxmin = -1;
			if (!max)
				maxmin = 1;

			myBrentMaxima->clear();
			myBrentMaxima->resize(dihedralset.size());

			for (unsigned int i = 0; i < dihedralset.size(); ++i)
			{
				//((*myBrentMaxima)[i]).resize(1);

				for (unsigned int j = 0; j <= 99; j++)
				{ //note the function evaluates one step past 2 pi
					Real lradangle = (Constant::M_PI * 2 / 100 * j);
					Real radangle = (Constant::M_PI * 2 / 100 * (j + 1));
					Real rradangle = (Constant::M_PI * 2 / 100 * (j + 2));

					Real valLangle = maxmin * computePhiDihedralEnergy(myTopology, dihedralset[i], lradangle);
					//report << hint << "Val left angle = " << valLangle << endr;
					Real valRangle = maxmin * computePhiDihedralEnergy(myTopology, dihedralset[i], rradangle);
					//report << "Val right angle = " << valRangle << endr;
					Real valAngle = maxmin * computePhiDihedralEnergy(myTopology, dihedralset[i], radangle);
					//report << "Val midpoint angle = " << valAngle << endr;

					Real xmax = 0.0;
					//Real ymax = 0.0;
					Real tol = 0.01;
					if ((valLangle > valAngle) && (valRangle > valAngle))
					{
						//report << "Suitable function segment" << endr;
						brent(lradangle, radangle, rradangle, tol, xmax, dihedralset[i], max);
						//report << std::endl << std::endl << endr;
						((*myBrentMaxima)[i]).push_back(xmax);
					}
				}
				// Throws Warning if no maxima were found
				if (((*myBrentMaxima)[i]).size() == 0)
				{
					report << warning << "No dihedral maxima found for dihedral index: " << dihedralset[i]
						<< " Check dihedral energy equation" << endr;
					((*myBrentMaxima)[i]).push_back(0.0);
				}
			}
			//temp hack that allows multiple calls but defeats the purpose of the cache... please fix me!
			//myCachedBrentMaxima= true;
		}
		return (*myBrentMaxima);
	}


	//BRENT FUNCTION
	Real OutputCache::brent(Real ax, Real bx, Real cx, Real tol, Real& xmin, int dihindex, bool max) const
	{
		const int ITMAX = 100;
		const Real CGOLD = 0.3819660;
		const Real ZEPS = std::numeric_limits<Real>::epsilon() * 1.0e-3;
		Real a, b, d = 0.0, etemp, fu, fv, fw, fx;
		Real p, q, r, tol1, tol2, u, v, w, x, xm;
		Real e = 0.0;

		//The Brent algorithm gets the maxima if maxmin = -1
		int maxmin = -1;
		if (!max)
			maxmin = 1;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = maxmin * computePhiDihedralEnergy(myTopology, dihindex, x);
		for (int iter = 0; iter < ITMAX; iter++)
		{
			xm = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
			if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
			{
				xmin = x;
				return fx;
			}
			if (fabs(e) > tol1)
			{
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q < 0.0) p = -p;
				q = fabs(q);
				etemp = e;
				e = d;
				if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				else
				{
					d = p / q;
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = sign(tol1, xm - x);
				}
			}
			else
			{
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (fabs(d) >= tol1 ? x + d : x + sign(tol1, d));
			fu = maxmin * computePhiDihedralEnergy(myTopology, dihindex, u);
			if (fu <= fx)
			{
				if (u >= x) a = x;
				else b = x;
				shift(v, w, x, u);
				shift(fv, fw, fx, fu);
			}
			else
			{
				if (u < x) a = u;
				else b = u;
				if (fu <= fw || w == x)
				{
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				}
				else if (fu <= fv || v == x || v == w)
				{
					v = u;
					fv = fu;
				}
			}
		}
		//too many iterations in brent
		xmin = x;
		return fx;
	}


	void OutputCache::uncache() const
	{
		myCachedKE = false;
		myCachedPE = false;
		myCachedV = false;
		myCachedP = false;
		myCachedMolP = false;
		myCachedLinearMomentum = false;
		myCachedAngularMomentum = false;
		myCachedCenterOfMass = false;
		myCachedDiffusion = false;
		myCachedDensity = false;
		myCachedMass = false;
		myCachedDihedralPhis = false;
		myCachedDihedralPhi = -1;
		myCachedBrentMaxima = false;
		myCachedMolT = false;
		myCachedMolKE = false;
		myCachedWaterT = false;
		myCachedNonWaterT = false;
		myCachedMinimalPositions = false;
		myCachedPositionsNoWater = false;
	}
}
