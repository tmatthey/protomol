/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFBORNFORCE_H
#define NONBONDEDCUTOFFBORNFORCE_H

#include "SystemForce.h"
#include "NonbondedCutoffForce.h"
#include "evaluateBorn.h"
#include "Report.h"
using namespace ProtoMol::Report;
#include <iostream>
#include <iomanip>
using namespace std;


namespace ProtoMol {
  //_________________________________________________________________ NonbondedCutoffSystemForce


  template<class TCellManager, class TOneAtomPair>
  class NonbondedCutoffBornForce: public NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffBornForce<TCellManager,TOneAtomPair> > {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    typedef typename TOneAtomPair::BoundaryConditions BoundaryConditions;
    typedef Topology<BoundaryConditions, TCellManager> RealTopologyType;
    typedef typename RealTopologyType::Enumerator EnumeratorType;
    typedef typename RealTopologyType::Enumerator::CellPair CellPairType;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedCutoffBornForce() : 
      NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffBornForce>(), keyword("NonbondedCutoffBorn") {}
    NonbondedCutoffBornForce(Real cutoff, TOneAtomPair oneAtomPair) : 
      NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffBornForce>(cutoff,oneAtomPair), keyword("NonbondedCutoffBorn") {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual void evaluate( GenericTopology*, 
			  const Vector3DBlock*, 
			  Vector3DBlock*, 
			  ScalarStructure*);
    virtual void parallelEvaluate( GenericTopology* topo, 
				  const Vector3DBlock* pos, 
				  Vector3DBlock* f, 
				  ScalarStructure* e);
    virtual std::string getKeyword() const{return keyword;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  public:
    const std::string keyword;
	Vector3DBlock bforces;

  };

  //const string NonbondedCutoffBornForce<TCellManager,TOneAtomPair>::keyword("NonbondedCutoffBorn");

  template<class TCellManager, class TOneAtomPair>
  void NonbondedCutoffBornForce<TCellManager,TOneAtomPair>::evaluate( GenericTopology* topo, 
								       const Vector3DBlock* positions, 
								       Vector3DBlock* forces, 
								       ScalarStructure* energies) {
    // Standard for forces
    if (!topo->doSCPISM)
      report << error << "To use SCPISM forces, please set doscpism in the configuration file" << endr;
    RealTopologyType* realTopo = dynamic_cast< RealTopologyType*>(topo);
    this->myOneAtomPair.getNonbondedForceFunction()->resize(topo->atoms.size());
    for (unsigned int i = 0; i < topo->atoms.size(); i++)
      this->myOneAtomPair.getNonbondedForceFunction()->resizeDR(i, topo->atoms.size());
    this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
    realTopo->updateCellLists(positions);
    this->enumerator.initialize(realTopo, this->myCutoff);
    this->doEvaluate(topo,realTopo->cellLists.size());
    evaluateBorn(this->myOneAtomPair, topo, forces, energies);
    
  }

  template<class TCellManager, class TOneAtomPair>
  void NonbondedCutoffBornForce<TCellManager,TOneAtomPair>::parallelEvaluate( GenericTopology* topo, 
									       const Vector3DBlock* positions, 
									       Vector3DBlock* forces, 
									       ScalarStructure* energies) {
  
     RealTopologyType* realTopo = dynamic_cast< RealTopologyType*>(topo);
     
    this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
    realTopo->updateCellLists(positions);
    this->enumerator.initialize(realTopo, this->myCutoff);

    unsigned int n = realTopo->cellLists.size();
    unsigned int count = numberOfBlocks(realTopo,positions);

    for(unsigned int i = 0;i<count;i++){
      unsigned int l = (n*(i+1))/count - (n*i)/count;
      if(Parallel::next())
	this->doEvaluate(topo,l);
      else
	this->enumerator.nextNewPair(l);	
    }
  }

}
#endif /* NONBONDEDCUTOFFSYSTEMFORCE_H */
