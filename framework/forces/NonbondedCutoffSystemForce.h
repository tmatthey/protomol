/* -*- c++ -*- */
#ifndef NONBONDEDCUTOFFSYSTEMFORCE_H
#define NONBONDEDCUTOFFSYSTEMFORCE_H

#include "SystemForce.h"
#include "NonbondedCutoffForce.h"
//#include "evaluateBorn.h"
namespace ProtoMol {
  //_________________________________________________________________ NonbondedCutoffSystemForce


  template<class TCellManager, class TOneAtomPair>
  class NonbondedCutoffSystemForce: public NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffSystemForce<TCellManager,TOneAtomPair> > {
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
    NonbondedCutoffSystemForce() : 
      NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffSystemForce>(){}
    NonbondedCutoffSystemForce(Real cutoff, TOneAtomPair oneAtomPair) : 
      NonbondedCutoffForce<TCellManager,TOneAtomPair,SystemForce,NonbondedCutoffSystemForce>(cutoff,oneAtomPair) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual void evaluate(const GenericTopology*, 
			  const Vector3DBlock*, 
			  Vector3DBlock*, 
			  ScalarStructure*);
//     virtual void evaluate(GenericTopology*, 
// 			  const Vector3DBlock*, 
// 			  Vector3DBlock*, 
// 			  ScalarStructure*);
    virtual void parallelEvaluate(const GenericTopology* topo, 
				  const Vector3DBlock* pos, 
				  Vector3DBlock* f, 
				  ScalarStructure* e);
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  private:
  };



  template<class TCellManager, class TOneAtomPair>
  void NonbondedCutoffSystemForce<TCellManager,TOneAtomPair>::evaluate(const GenericTopology* topo, 
								       const Vector3DBlock* positions, 
								       Vector3DBlock* forces, 
								       ScalarStructure* energies) {
  
    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);
    this->myOneAtomPair.initialize(realTopo, positions, forces, energies);
    realTopo->updateCellLists(positions);
    this->enumerator.initialize(realTopo, this->myCutoff);
    this->doEvaluate(topo,realTopo->cellLists.size());
  }

//   template<class TCellManager, class TOneAtomPair>
//   void NonbondedCutoffSystemForce<TCellManager,TOneAtomPair>::evaluate(GenericTopology* topo, 
// 								       const Vector3DBlock* positions, 
// 								       Vector3DBlock* forces, 
// 								       ScalarStructure* energies) {
  
//     evaluate(const_cast<GenericTopology*>(static_cast<const GenericTopology*>(topo)), positions, forces, energies);
//     if (born)
//       evaluateBorn(myOneAtomPair, topo, positions, forces, energies);
//   }

  template<class TCellManager, class TOneAtomPair>
  void NonbondedCutoffSystemForce<TCellManager,TOneAtomPair>::parallelEvaluate(const GenericTopology* topo, 
									       const Vector3DBlock* positions, 
									       Vector3DBlock* forces, 
									       ScalarStructure* energies) {
  
    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);

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
