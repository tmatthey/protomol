#include "EquilibriumMOLLYIntegrator.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "Topology.h"
#include <algorithm>
#include "Report.h"
using namespace ProtoMol::Report;
using std::find;
using std::vector;
using std::string;

namespace ProtoMol
{
	//_____________________________________________________EquilibriumMOLLYIntegrator

	const string EquilibriumMOLLYIntegrator::keyword("EquilibriumMOLLY");

	//const int EquilibriumMOLLYIntegrator::maxGDim;

	EquilibriumMOLLYIntegrator::EquilibriumMOLLYIntegrator(): MOLLYIntegrator(), myAveragedPositions(NULL)
	{
	}

	EquilibriumMOLLYIntegrator::EquilibriumMOLLYIntegrator(int cycles,
	                                                       ForceGroup* overloadedForces,
	                                                       StandardIntegrator* nextIntegrator)
		: MOLLYIntegrator(cycles, overloadedForces, nextIntegrator), myAveragedPositions(new Vector3DBlock())
	{
	}


	EquilibriumMOLLYIntegrator::~EquilibriumMOLLYIntegrator()
	{
		delete myAveragedPositions;
	}

	void EquilibriumMOLLYIntegrator::initialize(GenericTopology* topo,
	                                            Vector3DBlock* positions,
	                                            Vector3DBlock* velocities,
	                                            ScalarStructure* energies)
	{
		MOLLYIntegrator::initialize(topo, positions, velocities, energies);

		myMOLLYStepsize = getBottomTimestep() / Constant::TIMEFACTOR;
		// make myMOLLYStepsize same as deepest level (STS level) stepsize

		myAveragedPositions->resize(myPositions->size());

		vector<int> hAtoms;
		for (unsigned int i = 0; i < myTopo->atomTypes.size(); i++)
		{
			if (myTopo->atomTypes[i].name[0] == 'H')
			{
				hAtoms.push_back(i);
				//report << debug << "Using this as hydrogen : "<<myTopo->atomTypes[i].name <<endr; 
			}
		}

		// Find all Hydrogen Groups
		// Collect all bonds for each hydrogen group

		for (unsigned int i = 0; i < myTopo->bonds.size(); i++)
		{
			int a1 = myTopo->bonds[i].atom1;
			int a2 = myTopo->bonds[i].atom2;
			if (hAtoms.end() != find(hAtoms.begin(), hAtoms.end(), myTopo->atoms[a1].type) ||
				hAtoms.end() != find(hAtoms.begin(), hAtoms.end(), myTopo->atoms[a2].type))
			{
				bool append = true;
				HydrogenBond bond;
				bond.a1 = a1;
				bond.a2 = a2;
				bond.lambda = 0.0;
				bond.l12 = power<2>(myTopo->bonds[i].restLength);
				for (unsigned int k = 0; k < myHydrogenAtomGroups.size(); k++)
				{
					bool found1 = false;
					bool found2 = false;
					short i1 = myHydrogenAtomGroups[k].size();
					short i2 = myHydrogenAtomGroups[k].size();
					for (unsigned int l = 0; l < myHydrogenAtomGroups[k].size(); l++)
					{
						if (myHydrogenAtomGroups[k][l] == a1)
						{
							found1 = true;
							i1 = l;
						}
						else if (myHydrogenAtomGroups[k][l] == a2)
						{
							found2 = true;
							i2 = l;
						}
					}

					// Add the atom if we do not already have it in the list ...
					// ... and keep the bond.
					if (found1 && !found2)
					{
						myHydrogenAtomGroups[k].push_back(a2);
						bond.i1 = i1;
						bond.i2 = i2;
						myHydrogenConstraintGroups[k].push_back(bond);
						append = false;
						break;
					}
					else if (!found1 && found2)
					{
						myHydrogenAtomGroups[k].push_back(a1);
						bond.i1 = i1;
						bond.i2 = i2;
						myHydrogenConstraintGroups[k].push_back(bond);
						append = false;
						break;
					}
					else if (found1 && found2)
					{
						bond.i1 = i1;
						bond.i2 = i2;
						myHydrogenConstraintGroups[k].push_back(bond);
						append = false;
						break;
					}
				}

				// We found a new pair 
				if (append)
				{
					vector<int> v;
					v.push_back(a1);
					v.push_back(a2);
					myHydrogenAtomGroups.push_back(v);

					vector<HydrogenBond> bonds;
					bond.i1 = 0;
					bond.i2 = 1;
					bonds.push_back(bond);
					myHydrogenConstraintGroups.push_back(bonds);
				}
			}
		}

		// Check if the dimension is ok
		unsigned int max = 1;
		for (unsigned int k = 0; k < myHydrogenAtomGroups.size(); k++)
		{
			if (myHydrogenAtomGroups[k].size() > max)
				max = myHydrogenAtomGroups[k].size();
			if (myHydrogenConstraintGroups[k].size() > max)
				max = myHydrogenConstraintGroups[k].size();
		}
		if (max > (static_cast<unsigned int>(maxGDim)))
			report << error << "[EquilibriumMOLLYIntegrator::initialize] Dimension of G (" << maxGDim << ")"
				<< " smaller than the requested size " << max << "." << endr;

		initializeForces();
	}


	Vector3DBlock* EquilibriumMOLLYIntegrator::doAveragingPositions()
	{
		// Some user defined values
		const unsigned int ntrial = 100;
		const Real tolf = 0.00001;
		const Real tolx = tolf;

		// Set the myAveragedPositions to actual positions
		for (unsigned int i = 0; i < myPositions->size(); i++)
			(*myAveragedPositions)[i] = (*myPositions)[i];

		// Loop over all hydrogen groups
		for (unsigned int i = 0; i < myHydrogenConstraintGroups.size(); i++)
		{
			// Short cuts
			vector<HydrogenBond>& constraints = myHydrogenConstraintGroups[i];
			const unsigned int numConstraints = constraints.size();
			const vector<int>& atoms = myHydrogenAtomGroups[i];
			const unsigned int numAtoms = atoms.size();

			// Clear the lambda's and
			for (unsigned int j = 0; j < numConstraints; j++)
			{
				constraints[j].lambda = 0.0;
			}

			// ... and precompute the inverse mass for the actual atoms of the hydrogen group
			Real imass[maxGDim];
			for (unsigned int j = 0; j < numAtoms; j++)
			{
				imass[j] = 1.0 / myTopo->atoms[atoms[j]].scaledMass;
			}

			// Initial grhs
			Vector3D grhs[maxGDim][maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
				for (unsigned int l = 0; l < numAtoms; l++)
					grhs[k][l] = Vector3D(0.0, 0.0, 0.0);
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				const HydrogenBond& bond = constraints[k];
				grhs[k][bond.i1] = ((*myPositions)[bond.a1] - (*myPositions)[bond.a2]) * 2;
				grhs[k][bond.i2] = -grhs[k][bond.i1];
			}

			// Iterations
			for (unsigned int j = 0; j < ntrial; j++)
			{
				//report << debug << "Iteration :" <<j<<endr;

				Vector3D auxrhs[maxGDim][maxGDim];
				for (unsigned int k = 0; k < numConstraints; k++)
					for (unsigned int l = 0; l < numAtoms; l++)
						auxrhs[k][l] = grhs[k][l] * constraints[k].lambda * imass[l];

				Vector3D tmp[maxGDim];
				for (unsigned int l = 0; l < numAtoms; l++)
				{
					tmp[l] = Vector3D(0.0, 0.0, 0.0);
					for (unsigned int k = 0; k < numConstraints; k++)
						tmp[l] += auxrhs[k][l];
				}

				for (unsigned int l = 0; l < numAtoms; l++)
					(*myAveragedPositions)[atoms[l]] = (*myPositions)[atoms[l]] + tmp[l];

				Vector3D avgab[maxGDim];
				for (unsigned int k = 0; k < numConstraints; k++)
					avgab[k] = (*myAveragedPositions)[constraints[k].a1] -
						(*myAveragedPositions)[constraints[k].a2];

				Vector3D glhs[maxGDim][maxGDim];
				for (unsigned int k = 0; k < numConstraints; k++)
					for (unsigned int l = 0; l < numAtoms; l++)
						glhs[k][l] = Vector3D(0.0, 0.0, 0.0);
				for (unsigned int k = 0; k < numConstraints; k++)
				{
					const HydrogenBond& bond = constraints[k];
					glhs[k][bond.i1] = avgab[k] * 2;
					glhs[k][bond.i2] = -glhs[k][bond.i1];
				}

				// Update with the masse
				Real fjac[maxGDim][maxGDim];
				for (unsigned int k = 0; k < numConstraints; k++)
				{
					for (unsigned int l = 0; l < numConstraints; l++)
					{
						fjac[k][l] = 0.0;
						for (unsigned int m = 0; m < numAtoms; m++)
						{
							fjac[k][l] += glhs[k][m].dot(grhs[l][m]) * imass[m];
						}
					}
				}

				Real err = 0.0;
				//Real gij[maxGDim];
				//Real fvec[maxGDim];
				Real p[maxGDim];
				for (unsigned int k = 0; k < numConstraints; k++)
				{
					Real x = avgab[k].normSquared() - constraints[k].l12;
					//gij[k]  = x;
					//fvec[k] = x;
					p[k] = -x;
					err += fabs(x);
				}

				if (err < tolf)
					break;

				int index[maxGDim];
				Real d;
				luDcmp(fjac, numConstraints, index, d);
				luBksb(fjac, numConstraints, index, p);

				Real errx = 0.0;
				for (unsigned int k = 0; k < numConstraints; k++)
				{
					errx += fabs(p[k]);
					constraints[k].lambda += p[k];
				}

				if (errx < tolx)
					break;
			}
		}

		return myAveragedPositions;
	}

	void EquilibriumMOLLYIntegrator::doMollification(Vector3DBlock*)
	{
		// Loop over all hydrogen groups
		for (unsigned int i = 0; i < myHydrogenConstraintGroups.size(); i++)
		{
			// Short cuts
			vector<HydrogenBond>& constraints = myHydrogenConstraintGroups[i];
			const unsigned int numConstraints = constraints.size();
			const vector<int>& atoms = myHydrogenAtomGroups[i];
			const unsigned int numAtoms = atoms.size();

			// ... and precompute the inverse mass for the actual atoms of the hydrogen group

			Real imass[maxGDim];
			Vector3D tmpforce[maxGDim];
			for (unsigned int j = 0; j < numAtoms; j++)
			{
				imass[j] = 1.0 / myTopo->atoms[atoms[j]].scaledMass;
				tmpforce[j] = (*myForces)[atoms[j]] * imass[j];
			}

			Vector3D avgab[maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
				avgab[k] = (*myAveragedPositions)[constraints[k].a1] -
					(*myAveragedPositions)[constraints[k].a2];

			Vector3D grhs[maxGDim][maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
				for (unsigned int l = 0; l < numAtoms; l++)
					grhs[k][l] = Vector3D(0.0, 0.0, 0.0);
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				const HydrogenBond& bond = constraints[k];
				grhs[k][bond.i1] = ((*myPositions)[bond.a1] - (*myPositions)[bond.a2]) * 2;
				grhs[k][bond.i2] = -grhs[k][bond.i1];
			}

			Vector3D glhs[maxGDim][maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
				for (unsigned int l = 0; l < numAtoms; l++)
					glhs[k][l] = Vector3D(0.0, 0.0, 0.0);
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				const HydrogenBond& bond = constraints[k];
				glhs[k][bond.i1] = avgab[k] * 2;
				glhs[k][bond.i2] = -glhs[k][bond.i1];
			}

			// Update with the masse
			Real fjac[maxGDim][maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				for (unsigned int l = 0; l < numConstraints; l++)
				{
					fjac[k][l] = 0.0;
					for (unsigned int m = 0; m < numAtoms; m++)
					{
						fjac[k][l] += glhs[k][m].dot(grhs[l][m]) * imass[m];
					}
				}
			}

			Real aux[maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				aux[k] = 0.0;
				for (unsigned int l = 0; l < numAtoms; l++)
				{
					aux[k] += grhs[k][l].dot(tmpforce[l]);
				}
			}

			int index[maxGDim];
			Real d;
			luDcmp(fjac, numConstraints, index, d);

			Real fjacinv[maxGDim][maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
				for (unsigned int l = 0; l < numConstraints; l++)
					fjacinv[k][l] = (k == l) ? 1.0 : 0.0;
			for (unsigned int k = 0; k < numConstraints; k++)
				luBksb(fjac, numConstraints, index, fjacinv[k]);
			Real aux2[maxGDim];
			for (unsigned int k = 0; k < numConstraints; k++)
			{
				aux2[k] = 0.0;
				for (unsigned int l = 0; l < numConstraints; l++)
				{
					aux2[k] += fjacinv[k][l] * aux[l];
				}
			}
			for (unsigned int k = 0; k < numConstraints; k++)
				aux[k] = aux2[k];

			//luBksb(fjac,numConstraints,index,aux);


			Vector3D y[maxGDim];
			for (unsigned int l = 0; l < numAtoms; l++)
			{
				y[l] = Vector3D(0.0, 0.0, 0.0);
				for (unsigned int k = 0; k < numConstraints; k++)
				{
					y[l] += glhs[k][l] * aux[k];
				}
			}

			for (unsigned int l = 0; l < numAtoms; l++)
				y[l] = (*myForces)[atoms[l]] - y[l];

			Vector3D tmpforce2[maxGDim];
			for (unsigned int l = 0; l < numAtoms; l++)
				tmpforce2[l] = y[l] * imass[l];

			for (unsigned int l = 0; l < numAtoms; l++)
				tmpforce[l] = Vector3D(0.0, 0.0, 0.0);

			for (unsigned int k = 0; k < numConstraints; k++)
			{
				Vector3D tmpf = (tmpforce2[constraints[k].i1] - tmpforce2[constraints[k].i2]) * 2.0 * constraints[k].lambda;
				tmpforce[constraints[k].i1] += tmpf;
				tmpforce[constraints[k].i2] -= tmpf;
			}

			for (unsigned int l = 0; l < numAtoms; l++)
				(*myForces)[atoms[l]] = tmpforce[l] + y[l];
		}
	}

	void EquilibriumMOLLYIntegrator::luDcmp(Real (&m)[maxGDim][maxGDim],
	                                        int dim, int (&index)[maxGDim],
	                                        Real& d) const
	{
		int imax = 0;
		Real big, dum, sum, temp;
		Real vv[maxGDim];
		d = 1.0;
		for (int i = 0; i < dim; i++)
		{
			big = 0.0;
			for (int j = 0; j < dim; j++)
				if ((temp = fabs(m[i][j])) > big)
					big = temp;
			if (big == 0.0)
				report << error << "[EquilibriumMOLLYIntegrator::luDcmp] "
					<< "Singular matrix ... bye bye ... " << endr;
			vv[i] = 1.0 / big;
		}
		for (int j = 0; j < dim; j++)
		{
			for (int i = 0; i < j; i++)
			{
				sum = m[i][j];
				for (int k = 0; k < i; k++) sum -= m[i][k] * m[k][j];
				m[i][j] = sum;
			}
			big = 0.0;
			for (int i = j; i < dim; i++)
			{
				sum = m[i][j];
				for (int k = 0; k < j; k++)
					sum -= m[i][k] * m[k][j];
				m[i][j] = sum;
				if ((dum = vv[i] * fabs(sum)) >= big)
				{
					big = dum;
					imax = i;
				}
			}
			if (j != imax)
			{
				for (int k = 0; k < dim; k++)
				{
					dum = m[imax][k];
					m[imax][k] = m[j][k];
					m[j][k] = dum;
				}
				d = -d;
				vv[imax] = vv[j];
			}
			index[j] = imax;
			if (m[j][j] == 0.0) m[j][j] = Constant::TINY;
			if (j != dim - 1)
			{
				dum = 1.0 / (m[j][j]);
				for (int i = j + 1; i < dim; i++) m[i][j] *= dum;
			}
		}
	}

	MTSIntegrator* EquilibriumMOLLYIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const
	{
		return new EquilibriumMOLLYIntegrator(values[0], fg, nextIntegrator);
	}


	void EquilibriumMOLLYIntegrator::luBksb(Real (&m)[maxGDim][maxGDim],
	                                        int dim, const int (&index)[maxGDim],
	                                        Real (&b)[maxGDim]) const
	{
		int ii = -1;

		for (int i = 0; i < dim; i++)
		{
			int ip = index[i];
			Real sum = b[ip];
			b[ip] = b[i];
			if (ii >= 0)
				for (int j = ii; j < i; j++)
					sum -= m[i][j] * b[j];
			else if (sum)
				ii = i;
			b[i] = sum;
		}
		for (int i = dim - 1; i >= 0; i--)
		{
			Real sum = b[i];
			for (int j = i + 1; j < dim; j++)
				sum -= m[i][j] * b[j];
			b[i] = sum / m[i][i];
		}
	}
}
