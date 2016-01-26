/* -*- c++ -*- */
#ifndef ISGBONDSYSTEMFORCE_H
#define ISGBONDSYSTEMFORCE_H

#include "SystemForce.h"
#include "iSGBondSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"

namespace ProtoMol {
  //_________________________________________________________________ iSGBondSystemForce

  template<class TBoundaryConditions>
  class iSGBondSystemForce : public SystemForce, private iSGBondSystemForceBase{
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGBondSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void calcBond(const TBoundaryConditions &boundary,
		  const Bond& currentBond, 
		  const Vector3DBlock* positions, 
		  Vector3DBlock* forces,
		  ScalarStructure* energies);

    Real calcBondEnergy(const TBoundaryConditions &boundary,
			const Bond& currentBond, 
			const Vector3DBlock* positions);


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
    virtual std::string getKeyword() const{return keyword;}
    virtual unsigned int numberOfBlocks(const GenericTopology* topo,
					const Vector3DBlock* pos);
  private:
    virtual Force* doMake(std::string&, std::vector<Value>) const {
      return (new iSGBondSystemForce());
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 0;}
    virtual void getParameters(std::vector<Parameter>&) const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  }; 

  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions> 
  inline void iSGBondSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
							        const Vector3DBlock* positions, 
							        Vector3DBlock* forces,
							        ScalarStructure* energies)
  {
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;

    for (unsigned int i = 0; i < topo->bonds.size(); i++)
      calcBond(boundary, topo->bonds[i], positions, forces, energies);
  }

  template<class TBoundaryConditions> 
  inline void iSGBondSystemForce<TBoundaryConditions>::calcBond(const TBoundaryConditions &boundary,
							        const Bond& currentBond,
							        const Vector3DBlock* positions, 
							        Vector3DBlock* forces,
							        ScalarStructure* energies)
  {
    int a1 = currentBond.atom1;
    int a2 = currentBond.atom2;

    Real restLength     = currentBond.restLength;
    Real springConstant = currentBond.springConstant;
    Real DeltaK         = currentBond.DeltaK;
    Real DeltaR0        = currentBond.DeltaR0;

    Vector3D atom1 = (*positions)[a1];
    Vector3D atom2 = (*positions)[a2];
  
    //Vector3D r12 = atom1 - atom2;                       // Vector from atom 1 to atom 2.
    Vector3D r12 = boundary.minimalDifference(atom2,atom1); // Vector from atom 1 to atom 2.
    Real r       = r12.norm();                            // Distance between atom 1 and 2.

    Real dpotdr = 2.0 * springConstant * (r - restLength);   // Calculate dpot/dr
  
    // Calculate force on atom1 due to atom2.
    Vector3D force1(r12 * (-dpotdr/r));
  
    // Add to the total force.
    (*forces)[a1] += force1;
    (*forces)[a2] -= force1;
  
    // Add energy
    Real Vbond = springConstant * (r - restLength)*(r - restLength);
    (*energies)[ScalarStructure::BOND] += Vbond;

    // Add chemical potential difference
    if (DeltaK != 0.0 && DeltaR0 != 0.0)
      (*energies)[ScalarStructure::BOND_DELTAMU] += DeltaK / springConstant * Vbond - DeltaR0 * dpotdr;

//#define ISG_DEBUG_BOND
#ifdef ISG_DEBUG_BOND
    Report::report.precision(8);
    Report::report << "bond atoms: " << a1+1 << ", " << a2+1 << Report::endr;
    Report::report << "k = " << springConstant <<  ", r0 = " << restLength << Report::endr;
    Report::report << "DeltaK = " << DeltaK <<  ", DeltaR0 = " << DeltaR0 << Report::endr;
    Report::report << "bond length = " << r << Report::endr;
    Report::report << "current energy (bond) = "<< (*energies)::[ScalarStructure::BOND] << Report::endr;
    Report::report << "current DeltaMu (bond) = "<< (*energies)::[ScalarStructure::BOND_DELTAMU] << Report::endr;
#endif
    // Add virial
    energies->addVirial(force1,r12);
  }


  template<class TBoundaryConditions> 
  inline Real iSGBondSystemForce<TBoundaryConditions>::calcBondEnergy(const TBoundaryConditions &boundary,
								      const Bond& currentBond, 
								      const Vector3DBlock* positions)
  {

    int a1 = currentBond.atom1;
    int a2 = currentBond.atom2;
    Real restLength     = currentBond.restLength;
    Real springConstant = currentBond.springConstant;

    Vector3D atom1 = (*positions)[a1];
    Vector3D atom2 = (*positions)[a2];

    //Vector3D r12 = atom1 - atom2;                       // Vector from atom 1 to atom 2.
    Vector3D r12 = boundary.minimalDifference(atom2,atom1); // Vector from atom 1 to atom 2.
    Real r       = r12.norm();                            // Distance between atom 1 and 2.
 
    //Calculate energy.
    return (springConstant * ( r - restLength)*(r - restLength));
  }

  template<class TBoundaryConditions> 
  inline void iSGBondSystemForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
								        const Vector3DBlock* positions, 
								        Vector3DBlock* forces,
								        ScalarStructure* energies)
  {
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;

    unsigned int n = topo->bonds.size();
    unsigned int count = numberOfBlocks(topo,positions);
  
    for(unsigned int i = 0;i<count;i++){
      if(Parallel::next()){
	int to = (n*(i+1))/count;
	if(to > static_cast<int>(n))
	  to = n;
	int from = (n*i)/count;
	for (int j = from; j < to; j++)
	  calcBond(boundary, topo->bonds[j], positions, forces, energies);
      }
    }
  }

  template<class TBoundaryConditions> 
  inline unsigned int iSGBondSystemForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology* topo,
									      const Vector3DBlock*){    
    return std::min(Parallel::getAvailableNum(),static_cast<int>(topo->bonds.size()));
  }
}

#endif /* ISGBONDSYSTEMFORCE_H */
