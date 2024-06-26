#include "LeapfrogTruncatedShadow.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"
#include "ReducedHessAngle.h"
#include "reducedHessBond.h"
#include "ReducedHessLennardJones.h"
#include "LennardJonesForce.h"
#include "ReducedHessCoulomb.h"
#include "ReducedHessCoulombDiElec.h"
#include "ReducedHessCoulombSCPISM.h"
#include "ReducedHessCoulombBornRadii.h"
#include "CoulombForce.h"
#include "CoulombForceDiElec.h"
#include "CoulombSCPISMForce.h"
#include "CoulombBornRadiiForce.h"
#include "CnSwitchingFunction.h"
#include "C2SwitchingFunction.h"
#include "C1SwitchingFunction.h"

using namespace ProtoMol::Report;
using std::string;
using std::vector;
using std::map;

namespace ProtoMol
{
	//__________________________________________________ LeapfrogTruncatedShadow

	const string LeapfrogTruncatedShadow::keyword("LeapfrogTruncatedShadow");

	LeapfrogTruncatedShadow::LeapfrogTruncatedShadow() : STSIntegrator()
	{
	}

	LeapfrogTruncatedShadow::LeapfrogTruncatedShadow(Real timestep,
	                                                 ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces)
	{
		vector<Force*> ListForces = overloadedForces->getForces();
		//std::vector<Parameter> param;
		lCutoff = lSwitchon = lSwitch = cCutoff = cSwitchon = cSwitch = 0.0;
		lOrder = cOrder = lSwitchoff = cSwitchoff = 0.0;
		D = 78.0;
		S = 0.3;
		epsi = 1.0;
		myBond = myAngle = myCoulomb = myCoulombDielec = myCoulombSCPISM = myCoulombBornRadii = myLennardJones = myDihedral = myImproper = false;
		for (unsigned int i = 0; i < ListForces.size(); i++)
		{
			if (equalNocase(ListForces[i]->getId(), "Bond"))
			{
				myBond = true;
			}
			else if (equalNocase(ListForces[i]->getId(), "Angle"))
			{
				myAngle = true;
			}
			else if (equalStartNocase("CoulombBornRadii", ListForces[i]->getId()))
			{
				myCoulombBornRadii = true;
				swt = 1;
				swt = myTopo->doSCPISM;
				/*std::vector<Parameter> Fparam = ListForces[i]->getParameters();
                for(unsigned int j=0;j<Fparam.size();j++){

                   if(equalNocase(Fparam[j].keyword,"-bornswitch")){
                      swt=Fparam[j].value;
                   }
                }*/
			}
			else if (equalStartNocase("CoulombDiElec", ListForces[i]->getId()))
			{
				myCoulombDielec = true;
				std::vector<Parameter> Fparam = ListForces[i]->getParameters();
				for (unsigned int j = 0; j < Fparam.size(); j++)
				{
					if (equalNocase(Fparam[j].keyword, "-cutoff"))
					{
						cCutoff = Fparam[j].value;
						if (cSwitch == 0) cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchon"))
					{
						cSwitchon = Fparam[j].value;
						if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
						else if (equalStartNocase("Cn", Fparam[j].text))
						{
							cSwitch = 3;
						}
						else cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-n"))
					{
						cOrder = Fparam[j].value;
						cSwitch = 3;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchoff"))
					{
						cSwitchoff = Fparam[j].value;
						cSwitch = 3;
					}
					else if (equalNocase(Fparam[j].keyword, "-D"))
					{
						D = Fparam[j].value;
					}
					else if (equalNocase(Fparam[j].keyword, "-S"))
					{
						S = Fparam[j].value;
					}
					else if (equalNocase(Fparam[j].keyword, "-EPS"))
					{
						epsi = Fparam[j].value;
					}
				}
			}
			else if (equalStartNocase("CoulombSCPISM", ListForces[i]->getId()))
			{
				myCoulombSCPISM = true;
				std::vector<Parameter> Fparam = ListForces[i]->getParameters();
				for (unsigned int j = 0; j < Fparam.size(); j++)
				{
					if (equalNocase(Fparam[j].keyword, "-cutoff"))
					{
						cCutoff = Fparam[j].value;
						if (cSwitch == 0) cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchon"))
					{
						cSwitchon = Fparam[j].value;
						if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
						else if (equalStartNocase("Cn", Fparam[j].text))
						{
							cSwitch = 3;
						}
						else cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-n"))
					{
						cOrder = Fparam[j].value;
						cSwitch = 3;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchoff"))
					{
						cSwitchoff = Fparam[j].value;
						cSwitch = 3;
					}
				}
			}
			else if (equalStartNocase("Coulomb", ListForces[i]->getId()))
			{
				myCoulomb = true;
				std::vector<Parameter> Fparam = ListForces[i]->getParameters();
				for (unsigned int j = 0; j < Fparam.size(); j++)
				{
					if (equalNocase(Fparam[j].keyword, "-cutoff"))
					{
						cCutoff = Fparam[j].value;
						if (cSwitch == 0) cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchon"))
					{
						cSwitchon = Fparam[j].value;
						if (equalStartNocase("C2", Fparam[j].text)) cSwitch = 2;
						else if (equalStartNocase("Cn", Fparam[j].text))
						{
							cSwitch = 3;
						}
						else cSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-n"))
					{
						cOrder = Fparam[j].value;
						cSwitch = 3;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchoff"))
					{
						cSwitchoff = Fparam[j].value;
						cSwitch = 3;
					}
				}
			}
			else if (equalStartNocase("LennardJones", ListForces[i]->getId()))
			{
				myLennardJones = true;
				std::vector<Parameter> Fparam = ListForces[i]->getParameters();
				for (unsigned int j = 0; j < Fparam.size(); j++)
				{
					if (equalNocase(Fparam[j].keyword, "-cutoff"))
					{
						lCutoff = Fparam[j].value;
						if (lSwitch == 0) lSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchon"))
					{
						lSwitchon = Fparam[j].value;
						if (equalStartNocase("C2", Fparam[j].text)) lSwitch = 2;
						else if (equalStartNocase("Cn", Fparam[j].text))
						{
							lSwitch = 3;
						}
						else lSwitch = 1;
					}
					else if (equalNocase(Fparam[j].keyword, "-n"))
					{
						lOrder = Fparam[j].value;
						lSwitch = 3;
					}
					else if (equalNocase(Fparam[j].keyword, "-switchoff"))
					{
						lSwitchoff = Fparam[j].value;
						lSwitch = 3;
					}
				}
			}
			else if (equalStartNocase("Dihedral", ListForces[i]->getId()))
			{
				myDihedral = true;
			}
			else if (equalStartNocase("Improper", ListForces[i]->getId()))
			{
				myImproper = true;
			}
		}
		report << hint << "[LeapfrogTruncatedShadow::initialize] Coulumb cutoff=" << cCutoff << ", switchon=" << cSwitchon << ", switch=" << cSwitch << ", order=" << cOrder << ", sOff=" << cSwitchoff << endr;
		report << hint << "[LeapfrogTruncatedShadow::initialize] LJ cutoff=" << lCutoff << ", switchon=" << lSwitchon << ", switch=" << lSwitch << ", order=" << lOrder << ", sOff=" << lSwitchoff << endr;
		report << hint << "[LeapfrogTruncatedShadow::initialize] myCoulumb=" << myCoulomb << ", myCoulombDielec=" << myCoulombDielec << ", myCoulombSCPISM=" << myCoulombSCPISM << ", myCoulombBornRadii=" << myCoulombBornRadii << ", myLJ=" << myLennardJones << ", myBond=" << myBond << ", myAngle=" << myAngle << ", nyDihedral=" << myDihedral << ", myImproper=" << myImproper << endr;
		report << hint << "[LeapfrogIntCrs::initialize] myCoulombDielec=" << myCoulombDielec << ", D=" << D << ", S=" << S << ", epsi=" << epsi << endr;
	}


	LeapfrogTruncatedShadow::~LeapfrogTruncatedShadow()
	{
	}


	void LeapfrogTruncatedShadow::initialize(GenericTopology* topo,
	                                         Vector3DBlock* positions,
	                                         Vector3DBlock* velocities,
	                                         ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
	}

	void LeapfrogTruncatedShadow::doHalfKickdoDrift()
	{
		if (anyPreDriftOrNextModify())
		{
			doHalfKick();
			doDriftOrNextIntegrator();
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			for (unsigned int i = 0; i < count; ++i)
			{
				(*myVelocities)[i] += (*myForces)[i] * h * 0.5 / myTopo->atoms[i].scaledMass;
				(*myPositions)[i] += (*myVelocities)[i] * h;
			}
			buildMolecularCenterOfMass(myPositions, myTopo);
			buildMolecularMomentum(myVelocities, myTopo);
			postDriftOrNextModify();
		}
	}

	void LeapfrogTruncatedShadow::doKickdoDrift()
	{
		if (anyPreDriftOrNextModify() || anyPreStepModify() || anyPostStepModify())
		{
			if (anyPreStepModify() || anyPostStepModify())
			{
				doHalfKick();
				postStepModify();
				preStepModify();
				doHalfKick();
			}
			else
			{
				doKick();
			}
			doDriftOrNextIntegrator();
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			for (unsigned int i = 0; i < count; ++i)
			{
				(*myVelocities)[i] += (*myForces)[i] * h / myTopo->atoms[i].scaledMass;
				(*myPositions)[i] += (*myVelocities)[i] * h;
			}


			buildMolecularCenterOfMass(myPositions, myTopo);
			buildMolecularMomentum(myVelocities, myTopo);
			postDriftOrNextModify();
		}
	}

	void LeapfrogTruncatedShadow::run(int numTimesteps)
	{
		preStepModify();
		doHalfKickdoDrift();
		calculateForces();
		for (int i = 1; i < numTimesteps; i++)
		{
			doKickdoDrift();
			calculateForces();
		}
		doHalfKick();
		postStepModify();
		calculateForces();
		(*myEnergies)[ScalarStructure::SHADOW] = calcNearHamiltonian();
		//
	}

	STSIntegrator* LeapfrogTruncatedShadow::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new LeapfrogTruncatedShadow(values[0], fg);
	}

	Real LeapfrogTruncatedShadow::outputHessian()
	{
		return 0;
	}

	// calculate the Verlet 'nearby' Hamiltonian based on the 
	// Hessian matrices for bond and/or angle. 
	Real LeapfrogTruncatedShadow::calcNearHamiltonian()
	{
		Real accum, gl2;
		ReducedHessAngle rh;
		Matrix3by3 rha;
		Vector3D v[4], sm, zerov(0.0, 0.0, 0.0);
		int a1, a2, a3, a4;

		accum = 0.0;
		//Impropers
		if (myImproper)
		{
			double hessD[144];
			for (unsigned int i = 0; i < myTopo->impropers.size(); i++)
			{
				bool nonZForce = false; //test for force constants
				for (int j = 0; j < myTopo->impropers[i].multiplicity; j++)
				{
					if (myTopo->impropers[i].forceConstant[j]) nonZForce = true;
				}
				if (nonZForce)
				{
					//report <<hint<<"[LeapfrogTruncatedShadow::initialize] impropers= "<<i<<" multi= "<<myTopo->impropers[i].multiplicity<<" per="<<myTopo->impropers[i].periodicity[0]<<endr;
					TorsionHess(myTopo->impropers[i], hessD);
					a1 = myTopo->impropers[i].atom1;
					a2 = myTopo->impropers[i].atom2;
					a3 = myTopo->impropers[i].atom3;
					a4 = myTopo->impropers[i].atom4;
					//calc p^TM^{-1}(p^TM^{-1}U'')^T
					v[0] = (*myVelocities)[a1]; // p/m
					v[1] = (*myVelocities)[a2];
					v[2] = (*myVelocities)[a3];
					v[3] = (*myVelocities)[a4];
					for (int ii = 0; ii < 4; ii++)
					{
						sm = zerov; //zero sm
						for (int kk = 0; kk < 4; kk++)
						{
							int mi = kk * 36 + ii * 3;
							Matrix3by3 rhd(hessD[mi], hessD[mi + 1], hessD[mi + 2], hessD[mi + 12], hessD[mi + 13], hessD[mi + 14], hessD[mi + 24], hessD[mi + 25], hessD[mi + 26]);
							rhd.transpose();
							sm += rhd * v[kk];
						}
						accum += v[ii].dot(sm);
					}
				}
			}
		}
		//Dihedrals
		if (myDihedral)
		{
			double hessD[144];
			for (unsigned int i = 0; i < myTopo->dihedrals.size(); i++)
			{
				bool nonZForce = false; //test for force constants
				for (int j = 0; j < myTopo->dihedrals[i].multiplicity; j++)
				{
					if (myTopo->dihedrals[i].forceConstant[j]) nonZForce = true;
				}
				if (nonZForce)
				{
					TorsionHess(myTopo->dihedrals[i], hessD);//Dihedral(i, hessD);
					a1 = myTopo->dihedrals[i].atom1;
					a2 = myTopo->dihedrals[i].atom2;
					a3 = myTopo->dihedrals[i].atom3;
					a4 = myTopo->dihedrals[i].atom4;
					//calc p^TM^{-1}(p^TM^{-1}U'')^T
					v[0] = (*myVelocities)[a1]; // p/m
					v[1] = (*myVelocities)[a2];
					v[2] = (*myVelocities)[a3];
					v[3] = (*myVelocities)[a4];
					for (int ii = 0; ii < 4; ii++)
					{
						sm = zerov; //zero sm
						for (int kk = 0; kk < 4; kk++)
						{
							int mi = kk * 36 + ii * 3;
							Matrix3by3 rhd(hessD[mi], hessD[mi + 1], hessD[mi + 2], hessD[mi + 12], hessD[mi + 13], hessD[mi + 14], hessD[mi + 24], hessD[mi + 25], hessD[mi + 26]);
							rhd.transpose();
							sm += rhd * v[kk];
						}
						accum += v[ii].dot(sm);
					}
				}
			}
		}
		//Bonds
		if (myBond)
		{
			for (unsigned int i = 0; i < myTopo->bonds.size(); i++)
			{
				a1 = myTopo->bonds[i].atom1;
				a2 = myTopo->bonds[i].atom2;
				Real r_0 = myTopo->bonds[i].restLength;
				Real k = myTopo->bonds[i].springConstant;
				Matrix3by3 rha = reducedHessBond((*myPositions)[a1], (*myPositions)[a2], k, r_0);
				//calc p^TM^{-1}(p^TM^{-1}U'')^T
				v[0] = (*myVelocities)[a1]; // p/m
				v[1] = (*myVelocities)[a2];
				sm = zerov; //zero sm
				rha.transpose();
				for (int kk = 0; kk < 2; kk++)
				{
					if (kk) rha *= -1;
					sm += rha * v[kk];
				}
				accum += v[0].dot(sm);
				accum += v[1].dot(-sm);
			}
		}
		//Angles
		if (myAngle)
		{
			for (unsigned int i = 0; i < myTopo->angles.size(); i++)
			{
				//report <<hint<<"[LeapfrogTruncatedShadow::initialize] ang= "<<i<<" ubc="<<myTopo->angles[i].ureyBradleyConstant<<" urrl= "<<myTopo->angles[i].ureyBradleyRestLength<<endr;
				a1 = myTopo->angles[i].atom1;
				a2 = myTopo->angles[i].atom2;
				a3 = myTopo->angles[i].atom3;
				Real theta0 = myTopo->angles[i].restAngle;
				Real k_t = myTopo->angles[i].forceConstant;
				Real ubConst = myTopo->angles[i].ureyBradleyConstant;
				Real ubRestL = myTopo->angles[i].ureyBradleyRestLength;
				// ReducedHessAngle for atoms a1, a2 and a3
				rh.evaluate((*myPositions)[a1], (*myPositions)[a2], (*myPositions)[a3], k_t, theta0);
				//ureyBradley
				if (ubConst)
				{
					//Cheat using bond hessian as same as UB!!!!
					Matrix3by3 ubm = reducedHessBond((*myPositions)[a1], (*myPositions)[a3], ubConst, ubRestL);
					rh.accumulateTo(0, 0, ubm);
					rh.accumulateTo(2, 2, ubm);
					rh.accumulateNegTo(2, 0, ubm);
					rh.accumulateNegTo(0, 2, ubm);
				}
				//calc p^TM^{-1}(p^TM^{-1}U'')^T
				v[0] = (*myVelocities)[a1]; // p/m
				v[1] = (*myVelocities)[a2];
				v[2] = (*myVelocities)[a3];
				for (int ii = 0; ii < 3; ii++)
				{
					//for (int d1 = 0; d1 < 3; d1++) sm[d1] = 0.0;
					sm = zerov; //zero sm
					for (int kk = 0; kk < 3; kk++)
					{
						rha = rh(kk, ii);
						rha.transpose();
						sm += rha * v[kk];
					}
					accum += v[ii].dot(sm);
				}
			}
		}
		//Lennard Jones forces
		if (myLennardJones)
		{
			accum += calcPairInteractionHess(lSwitch, lSwitchon, lCutoff, lOrder, lSwitchoff, 0);
		}
		//Coulombic forces
		if (myCoulomb)
		{
			Real coulaccum;
			coulaccum = calcPairInteractionHess(cSwitch, cSwitchon, cCutoff, cOrder, cSwitchoff, 1);
			accum += coulaccum;
		}
		//Coulombic forces for implicit solvent
		if (myCoulombDielec)
		{
			Real couldielecaccum;
			couldielecaccum = calcPairInteractionHess(cSwitch, cSwitchon, cCutoff, cOrder, cSwitchoff, 2);
			accum += couldielecaccum;
		}
		if (myCoulombSCPISM)
		{
#if 0
		Real scpismAccum;
		scpismAccum = calcPairInteractionHess(cSwitch,cSwitchon,cCutoff,cOrder,cSwitchoff,3);
		accum += scpismAccum;
			//accum += calcPairInteractionHess(cSwitch,cSwitchon,cCutoff,cOrder,cSwitchoff,3);
#endif
		}
		if (myCoulombBornRadii)
		{
#if 1
			Real bornAccum;
			bornAccum = calcPairInteractionHess(cSwitch, cSwitchon, cCutoff, cOrder, cSwitchoff, 4);
			accum += bornAccum;
#endif
		}
		//
		gl2 = 2.0 * accum;
		//calc (U')^TM^{-1}U'
		accum = 0.0;
		//


		for (unsigned int i = 0; i < (*myForces).size(); i++)
		{
			accum += (*myForces)[i].dot((*myForces)[i]) / myTopo->atoms[i].scaledMass;
			//accum += (*myVelocities)[i].dot((*myVelocities)[i])/(myTopo->atoms[i].scaledMass*myTopo->atoms[i].scaledMass);
		}
		gl2 -= accum;
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		gl2 *= h * h / 24.0;
		//add in rest of energy
		gl2 += kineticEnergy(myTopo, myVelocities) + myEnergies->potentialEnergy();
		return gl2;
	}

	Real LeapfrogTruncatedShadow::calcPairInteractionHess(Real doSwitch, Real switchon, Real cutoff, Real order, Real switchoff, int lj)
	{
		Real accum;
		Matrix3by3 rha;
		Vector3D v[3], sm, zerov(0.0, 0.0, 0.0);
		ExclusionClass ec;

		accum = 0.0;
		for (unsigned int i = 0; i < myTopo->atoms.size(); i++)
		{
			for (unsigned int j = i + 1; j < myTopo->atoms.size(); j++)
			{
				//if not bonded/dihedral
				ec = myTopo->exclusions.check(i, j);
				if (ec != EXCLUSION_FULL)
				{
					Vector3D rij = myTopo->minimalDifference((*myPositions)[i], (*myPositions)[j]);
					//Vector3D rij = (*myPositions)[j]-(*myPositions)[i];
					Matrix3by3 mz(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
					Real a = rij.normSquared();
					if (lj == 0)
					{ // Lennard Jones pairs
						ReducedHessLennardJones rHess;
						LennardJonesForce hForce;
						Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
						if (doSwitch)
						{
							if (doSwitch == 3)
							{
								CnSwitchingFunction cnsf(switchon, cutoff, order, switchoff);
								cnsf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = cnsf.hessian(rij, a);
								}
							}
							else if (doSwitch == 2)
							{
								C2SwitchingFunction c2sf(switchon, cutoff);
								c2sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c2sf.hessian(rij, a);
								}
							}
							else
							{
								C1SwitchingFunction c1sf(cutoff);
								c1sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c1sf.hessian(rij, a);
								}
							}
						}
						rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec);
					}
					else if (lj == 2)
					{ // Implicit Colombic
						ReducedHessCoulombDielec rHess;
						CoulombForceDiElec hForce;
						Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
						if (doSwitch)
						{
							if (doSwitch == 3)
							{
								CnSwitchingFunction cnsf(switchon, cutoff, order, switchoff);
								cnsf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = cnsf.hessian(rij, a);
								}
							}
							else if (doSwitch == 2)
							{
								C2SwitchingFunction c2sf(switchon, cutoff);
								c2sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c2sf.hessian(rij, a);
								}
							}
							else
							{
								C1SwitchingFunction c1sf(cutoff);
								c1sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c1sf.hessian(rij, a);
								}
							}
						}
						rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec, D, S, epsi);
					}
					else if (lj == 3)
					{
						ReducedHessCoulombSCPISM rHess;
						CoulombSCPISMForce hForce;
						Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
						if (!myTopo->atoms[i].mySCPISM ||
							!myTopo->atoms[j].mySCPISM)
							report << error << "[LeapfrogTruncatedShadow::calcPairInteractionHess] SCPISM data not set." << endr;
						Real alpha_ij = myTopo->atoms[i].mySCPISM->sqrtalphaSCPISM * myTopo->atoms[j].mySCPISM->sqrtalphaSCPISM;
						if (doSwitch)
						{
							if (doSwitch == 3)
							{
								CnSwitchingFunction cnsf(switchon, cutoff, order, switchoff);
								cnsf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = cnsf.hessian(rij, a);
								}
							}
							else if (doSwitch == 2)
							{
								C2SwitchingFunction c2sf(switchon, cutoff);
								c2sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c2sf.hessian(rij, a);
								}
							}
							else
							{
								C1SwitchingFunction c1sf(cutoff);
								c1sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c1sf.hessian(rij, a);
								}
							}
						}
						rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec, alpha_ij, 80.0);
					}
					else if (lj == 4)
					{
						ReducedHessCoulombBornRadii rHess;
						CoulombBornRadiiForce hForce(swt);
						Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
						rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec, myPositions, hForce);
					}
					else
					{ // Colombic
						ReducedHessCoulomb rHess;
						CoulombForce hForce;
						Real rawE = 0.0, rawF = 0.0, swtchV = 1.0, swtchD = 0.0;
						if (doSwitch)
						{
							if (doSwitch == 3)
							{
								CnSwitchingFunction cnsf(switchon, cutoff, order, switchoff);
								cnsf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = cnsf.hessian(rij, a);
								}
							}
							else if (doSwitch == 2)
							{
								C2SwitchingFunction c2sf(switchon, cutoff);
								c2sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c2sf.hessian(rij, a);
								}
							}
							else
							{
								C1SwitchingFunction c1sf(cutoff);
								c1sf(swtchV, swtchD, a);
								if (swtchV > 0.0 && swtchV < 1.0)
								{
									hForce(rawE, rawF, a, 1.0 / a, rij, myTopo, i, j, ec);
									mz = c1sf.hessian(rij, a);
								}
							}
						}
						rha = rHess(rawE, rawF, a, a, rij, myTopo, i, j, swtchV, swtchD, mz, ec);
					}
					v[0] = (*myVelocities)[i]; // p/m
					v[1] = (*myVelocities)[j];
					sm = zerov; //zero sm
					rha.transpose();
					for (int kk = 0; kk < 2; kk++)
					{
						if (kk) rha *= -1;
						sm += rha * v[kk];
					}
					accum += v[0].dot(sm);
					accum += v[1].dot(-sm);
				}
			}
		}
		return accum;
	}

	void LeapfrogTruncatedShadow::TorsionHess(const Torsion& currTorsion, double* hessD)
	{
		//
		//ofstream myFile;
		int a1 = currTorsion.atom1;
		int a2 = currTorsion.atom2;
		int a3 = currTorsion.atom3;
		int a4 = currTorsion.atom4;

		//actual positions
		double a[9] = {(*myPositions)[a1][0],(*myPositions)[a2][0],(*myPositions)[a3][0],
			(*myPositions)[a1][1],(*myPositions)[a2][1],(*myPositions)[a3][1],
			(*myPositions)[a1][2],(*myPositions)[a2][2],(*myPositions)[a3][2]};
		//Vector3D rxy // Vector from atom a to atom b
		Vector3D r12(myTopo->minimalDifference((*myPositions)[a2], (*myPositions)[a1]));
		Vector3D r23v(myTopo->minimalDifference((*myPositions)[a3], (*myPositions)[a2]));
		Vector3D r34(myTopo->minimalDifference((*myPositions)[a4], (*myPositions)[a3]));
		// Cross product of r12 and r23v, represents the plane shared by these two vectors
		double cosPhi = -r12.dot(r23v) / (r12.norm() * r23v.norm());
		double sinPhi = (r12.cross(r23v)).norm() / (r12.norm() * r23v.norm());//  sin(acos(cosPhi));
		//double sinPhi= sin(acos(cosPhi));
		double nsa1 = (*myPositions)[a1].normSquared();
		double nsa2 = (*myPositions)[a2].normSquared();
		double nsa3 = (*myPositions)[a3].normSquared();
		//
		double x, y, z;
		z = (nsa3 - nsa2 - r23v.norm() * r23v.norm()) / (-2.0 * r23v.norm());
		x = (nsa1 - nsa2 - r12.norm() * r12.norm() + 2.0 * z * r12.norm() * cosPhi) / (-2.0 * r12.norm() * sinPhi);
		y = sqrt(nsa2 - x * x - z * z);
		//target positions
		double opC[9] = {x - r12.norm() * sinPhi, x, x,
			y, y, y,
			z - r12.norm() * cosPhi, z, z - r23v.norm()};
		//DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
		double dta = a[0] * (a[8] * a[4] - a[7] * a[5]) - a[3] * (a[8] * a[1] - a[7] * a[2]) + a[6] * (a[5] * a[1] - a[4] * a[2]);
		//| a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
		//| a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
		//| a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
		//inverse of a
		double ia[9] = {a[8] * a[4] - a[7] * a[5], -(a[8] * a[1] - a[7] * a[2]), a[5] * a[1] - a[4] * a[2],
			-(a[8] * a[3] - a[6] * a[5]), a[8] * a[0] - a[6] * a[2], -(a[5] * a[0] - a[3] * a[2]),
			a[7] * a[3] - a[6] * a[4], -(a[7] * a[0] - a[6] * a[1]), a[4] * a[0] - a[3] * a[1]};
		for (int i = 0; i < 9; i++) ia[i] /= dta;
		//rotation matrix
		double aRot[9];
		for (int i = 0; i < 9; i++) aRot[i] = 0.0;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 3; j++) aRot[i] += opC[(i / 3) * 3 + j] * ia[i % 3 + j * 3];
		}
		//fourth body position
		double in4[3] = {(*myPositions)[a4][0],(*myPositions)[a4][1],(*myPositions)[a4][2]};
		//rotate it
		double out4[3] = {0.0, 0.0, 0.0};
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++) out4[i] += aRot[i * 3 + j] * in4[j];
		}
		double dgd[12] = {0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0};
		// Calculate phi
		//#######################
		double g;
		g = -atan((out4[1] - opC[5]) / (out4[0] - opC[2]));
		//dfdxdy := -k*cos(n*g(x, y, z)+ps)*n^2*(diff(g(x, y, z), y))*(diff(g(x, y, z), x))
		//				-k*sin(n*g(x, y, z)+ps)*n*(diff(g(x, y, z), y, x))
		//multiplicity
		double fact1 = 0.0;
		double fact2 = 0.0;
		for (int dm = 0; dm < currTorsion.multiplicity; dm++)
		{ //Dihedrals
			if (currTorsion.periodicity[dm] > 0)
			{
				fact1 += -currTorsion.forceConstant[dm] *
					sin(currTorsion.periodicity[dm] * g + currTorsion.phaseShift[dm]) *
					currTorsion.periodicity[dm];
				fact2 += -currTorsion.forceConstant[dm] *
					cos(currTorsion.periodicity[dm] * g + currTorsion.phaseShift[dm]) *
					currTorsion.periodicity[dm] * currTorsion.periodicity[dm];
			}
			else
			{ //Impropers
				Real diff = g - currTorsion.phaseShift[dm];
				if (diff < -Constant::M_PI)
					diff += 2 * Constant::M_PI;
				else if (diff > Constant::M_PI)
					diff -= 2 * Constant::M_PI;
				fact1 += 2.0 * currTorsion.forceConstant[dm] * diff;
				fact2 += 2.0 * currTorsion.forceConstant[dm];
			}
		}
		//Setup variables
		Real x21, x21_2, x43, x43_2, x43_3, x43_4, x43_5, y43, y43_2, y43_3, z21, z21_2, z43, z43_2, r23, r23_2;
		x21 = opC[1] - opC[0];
		x43 = out4[0] - opC[2];
		y43 = out4[1] - opC[5];
		z21 = opC[7] - opC[6];
		z43 = out4[2] - opC[8];
		r23 = r23v.norm();
		x21_2 = x21 * x21;
		x43_2 = x43 * x43;
		y43_2 = y43 * y43;
		r23_2 = r23 * r23;
		z21_2 = z21 * z21;
		z43_2 = z43 * z43;
		x43_3 = x43_2 * x43;
		y43_3 = y43_2 * y43;
		x43_4 = x43_2 * x43_2;
		x43_5 = x43_2 * x43_3;
		//a1 -force_x=0
		//   -force_y=-1/x21
		dgd[1] = -1 / x21;
		//a4 -force_x=y43/(x43_2*(1+y43_2/x43_2))
		//   -force_y=-1/(x43*(1+y43_2/x43_2))
		dgd[9] = y43 / (x43_2 * (1 + y43_2 / x43_2));
		dgd[10] = -1 / (x43 * (1 + y43_2 / x43_2));
		//a2 -force_x=-y43*z43/(x43_2*r23*(1+y43_2/x43_2))
		//   -force_y=(1-z21/r23)/x21+z43/(r23*x43*(1+y43_2/x43_2))
		dgd[3] = -y43 * z43 / (x43_2 * r23 * (1 + y43_2 / x43_2));
		dgd[4] = (1 - z21 / r23) / x21 + z43 / (r23 * x43 * (1 + y43_2 / x43_2));
		//a3 -force_x=y43*(-1+z43/r23)/(x43_2*(1+y43_2/x43_2))
		//   -force_y=z21/(r23*x21)-(-1+z43/r23)/(x43*(1+y43_2/x43_2))
		dgd[6] = y43 * (-1 + z43 / r23) / (x43_2 * (1 + y43_2 / x43_2));
		dgd[7] = z21 / (r23 * x21) - (-1 + z43 / r23) / (x43 * (1 + y43_2 / x43_2));
		//Hessian 
		//dfdxdy := -k*cos(n*g(x, y, z)+ps)*n^2*(diff(g(x, y, z), y))*(diff(g(x, y, z), x))
		//				-k*sin(n*g(x, y, z)+ps)*n*(diff(g(x, y, z), y, x))
		//create output hessian in x_1 y_1 z_1
		//                         y_1
		//                         z_1
		//format
		for (int i = 0; i < 144; i++) hessD[i] = 0.0;
		//
		//common factors
		double sqTerm1 = (1 + y43_2 / x43_2) * (1 + y43_2 / x43_2);
		double sqTerm2 = (-1 + z43 / r23) * (-1 + z43 / r23);
		//#d2x1y1
		hessD[0 * 12 + 1] = hessD[1 * 12 + 0] = -1 / x21_2;
		//d2x1y2
		hessD[0 * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 0] = (1 - z21 / r23) / x21_2;
		//d2x1y3
		hessD[0 * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 0] = z21 / (r23 * x21_2);
		//**d2x2
		hessD[3 * 12 + 3] = -2 * y43 * z43_2 / (x43_3 * r23_2 * (1 + y43_2 / x43_2)) -
			y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * y43_3 * z43_2 / (x43_5 * r23_2 * sqTerm1);//changed 2*y43/
		//**d2x3
		hessD[6 * 12 + 6] = -2 * y43 * sqTerm2 / (x43_3 * (1 + y43_2 / x43_2)) -
			y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * y43_3 * sqTerm2 / (x43_5 * sqTerm1); //changed from 2*y43/(
		//#d2x4
		hessD[9 * 12 + 9] = -2 * y43 / (x43_3 * (1 + y43_2 / x43_2)) + 2 * y43_3 / (x43_5 * sqTerm1);
		//**d2y2
		hessD[(3 + 1) * 12 + 3 + 1] = y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * z43_2 * y43 / (r23_2 * x43_3 * sqTerm1);//2*
		//**d2y3
		hessD[(6 + 1) * 12 + 6 + 1] = y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * sqTerm2 * y43 / (x43_3 * sqTerm1);//2*
		//d2y4
		hessD[(9 + 1) * 12 + 9 + 1] = 2 * y43 / (x43_3 * sqTerm1);
		//
		//**d2x2x3
		hessD[(3) * 12 + 6] = hessD[(6) * 12 + 3] = 2 * y43 * z43 * (-1 + z43 / r23) / (x43_3 * r23 * (1 + y43_2 / x43_2)) +
			y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) - 2 * y43_3 * z43 * (-1 + z43 / r23) / (x43_5 * r23 * sqTerm1);//2*
		//#d2x2x4
		hessD[(3) * 12 + 9] = hessD[(9) * 12 + 3] = 2 * y43 * z43 / (x43_3 * r23 * (1 + y43_2 / x43_2)) - 2 * y43_3 * z43 / (x43_5 * r23 * sqTerm1);
		//%%%%d2x2y1
		hessD[(3) * 12 + 1] = hessD[(1) * 12 + 3] = (1 - z21 / r23) / x21_2; //changed from 1/x21_2
		//%%%%d2x2y2
		hessD[(3) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 3] =
			-(1 - z21 / r23) * (1 - z21 / r23) / x21_2 + z43_2 / (r23_2 * x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 * z43_2 / (x43_4 * r23_2 * sqTerm1); //changed from -(1-z21/r23)/x21_2+
		//**d2x2y3
		hessD[(3) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 3] = -z21 * (1 - z21 / r23) / (r23 * x21_2) - (-1 + z43 / r23) * z43 / (x43_2 * r23 * (1 + y43_2 / x43_2)) +
			2 * y43_2 * z43 * (-1 + z43 / r23) / (x43_4 * r23 * sqTerm1);
		//#d2x2y4
		hessD[(3) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 3] = -z43 / (x43_2 * (1 + y43_2 / x43_2) * r23) + 2 * y43_2 * z43 / (x43_4 * sqTerm1 * r23); //=x4y2
		//
		//#d2x3x4
		hessD[(6) * 12 + 9] = hessD[(9) * 12 + 6] = -2 * y43 * (-1 + z43 / r23) / (x43_3 * (1 + y43_2 / x43_2)) + 2 * y43_3 * (-1 + z43 / r23) / (x43_5 * sqTerm1);
		//d2x3y2!!!!!!!!!!!
		hessD[(6) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 6] =
			-z43 * (-1 + z43 / r23) / (r23 * x43_2 * (1 + y43_2 / x43_2)) + 2 * y43_2 * (-1 + z43 / r23) * z43 / (x43_4 * sqTerm1 * r23)
			- z21 * (1 - z21 / r23) / (x21_2 * r23);
		//d2x3y3!!!!!!!!!!!
		hessD[(6) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 6] =
			sqTerm2 / (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 * sqTerm2 / (x43_4 * sqTerm1)
			- z21_2 / (r23_2 * x21_2);
		//#d2x3y4
		hessD[(6) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 6] = (-1 + z43 / r23) / (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 * (-1 + z43 / r23) / (x43_4 * sqTerm1); //=x4y3
		//##
		//#d2x4y2
		//hessD[(9)*12 + 3+1] = hessD[(3)*12 + 9+1]; //d2x2y4; //-z43/(r23*x43^2*(1+y43^2/x43^2))+2*y43^2*z43/(x43^4*(1+y43^2/x43^2)^2*r23) //=x2y4
		hessD[(9) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 9] = -z43 / (r23 * x43_2 * (1 + y43_2 / x43_2)) + 2 * y43_2 * z43 / (x43_4 * sqTerm1 * r23); //=x2y4
		//#d2x4y3
		hessD[(9) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 9] = hessD[(6) * 12 + 9 + 1]; //d2x3y4; //(-1+z43/r23)/(x43^2*(1+y43^2/x43^2))-2*y43^2*(-1+z43/r23)/(x43^4*(1+y43^2/x43^2)^2) //=x3y4
		//#d2x4y4
		hessD[(9) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 9] = 1 / (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 / (x43_4 * sqTerm1);
		//
		//**d2y2y3
		hessD[(3 + 1) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 3 + 1] =
			-y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) - 2 * z43 * y43 * (-1 + z43 / r23) / (r23 * x43_3 * sqTerm1);//2*
		//#d2y2y4
		hessD[(3 + 1) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 3 + 1] = -2 * z43 * y43 / (r23 * x43_3 * sqTerm1);
		//
		//#d2y3y4
		hessD[(6 + 1) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 6 + 1] = 2 * (-1 + z43 / r23) * y43 / (x43_3 * sqTerm1);
		//
		//**d2x2z3
		double temphd = y43 / (x43_2 * r23 * (1 + y43_2 / x43_2)); //same but +/-
		double temphdn = y43 * (1 - z43 / r23) / (x43_2 * r23 * (1 + y43_2 / x43_2));
		hessD[(3) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 3] = temphdn;//temphd;
		//#d2x2z4
		hessD[(3) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 3] = -temphd; //-d2x2z3; //-y43/(x43^2*r23*(1+y43^2/x43^2))
		//**d2x3z3
		hessD[(6) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 6] = -temphdn;//-temphd; //-d2x2z3; //-y43/(x43^2*r23*(1+y43^2/x43^2))
		//#d2x3z4
		hessD[(6) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 6] = temphd; //d2x2z3;	   //y43/(x43^2*r23*(1+y43^2/x43^2))
		//
		//d2y2z1
		temphd = 1 / (r23 * x21);
		hessD[(3 + 1) * 12 + 2] = hessD[(2) * 12 + 3 + 1] = temphd; //1/(r23*x21);
		//**Temp for y2z2 and -y3z2
		double temphd3 = -(1 / r23 - z21 / r23_2) / x21 - z43 / (r23_2 * x43 * (1 + y43_2 / x43_2));
		//** Tempf for y3z3 and -y2z3
		double temphd4 = z21 / (r23_2 * x21) - (-1 / r23 + z43 / r23_2) / (x43 * (1 + y43_2 / x43_2));
		//**d2y2z2
		hessD[(3 + 1) * 12 + 3 + 2] = hessD[(3 + 2) * 12 + 3 + 1] = temphd3;//-temphd; //-d2y2z1; //-1/(r23*x21)
		//#d2y2z4
		double temphd2 = 1 / (r23 * x43 * (1 + y43_2 / x43_2));
		hessD[(3 + 1) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 3 + 1] = temphd2;
		//**d2y2z3
		hessD[(3 + 1) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 3 + 1] = -temphd4;//-temphd2; //-d2y2z4; //-1/(r23*x43*(1+y43^2/x43^2))
		//
		//d2y3z1
		hessD[(6 + 1) * 12 + 2] = hessD[(2) * 12 + 6 + 1] = -temphd; //-d2y2z1; //-1/(r23*x21)
		//**d2y3z2
		hessD[(6 + 1) * 12 + 3 + 2] = hessD[(3 + 2) * 12 + 6 + 1] = -temphd3;//temphd; //d2y2z1; // 1/(r23*x21)
		//**d2y3z3
		hessD[(6 + 1) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 6 + 1] = temphd4;//temphd2; //d2y2z4; //1/(r23*x43*(1+y43^2/x43^2))
		//#d2y3z4
		hessD[(6 + 1) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 6 + 1] = -temphd2; //-d2y2z4; //-1/(r23*x43*(1+y43^2/x43^2))
		//d2x3y1!!!!!!!!!!
		hessD[(6) * 12 + 1] = hessD[(1) * 12 + 6] = z21 / (x21_2 * r23);
		///////
		int eye, eyeh, jay, kay;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{ //####for test! should be j=i; so 1/2 matrix
				//create the i,j hessian matrix
				eye = i * 3;
				eyeh = eye * 12;
				jay = j * 3;
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						kay = eyeh + jay + k * 12 + l;// (k/3)*12 + jay + k%3;
						//if(i==1 && j==0 && k==0 && l==1) 
						//	report <<hint<<"[HessianInt::initialize] hess= "<<hessD[kay]<<" dgd x2 = "<<dgd[eye+k]<<" dgd y1 = "<<dgd[jay+l]<<endr;
						hessD[kay] = hessD[kay] * fact1 + dgd[eye + k] * dgd[jay + l] * fact2;//
					}
				}
				for (int k = 0; k < 3; k++)
				{
					double v[3] = {hessD[eyeh + jay + k], hessD[eyeh + jay + k + 12], hessD[eyeh + jay + k + 24]};
					rotateV3D(aRot, v);
					hessD[eyeh + jay + k] = v[0];
					hessD[eyeh + jay + k + 12] = v[1];
					hessD[eyeh + jay + k + 24] = v[2];
				}
				for (int k = 0; k < 3; k++)
				{
					rotateV3D(aRot, &hessD[eyeh + jay + 12 * k]);
				}
			}
		}
	}

	//Use aRot to rotate the vector back into real space
	double* LeapfrogTruncatedShadow::rotateV3D(double* aRot, double* mf)
	{
		double out[3] = {0.0,0.0,0.0};
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++) out[i] += aRot[i + j * 3] * mf[j];
		}
		for (int i = 0; i < 3; i++) mf[i] = out[i];
		return mf;
	}
}
