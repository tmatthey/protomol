/* -*- c++ -*- */
#ifndef ISGIMPROPERSYSTEMFORCE_H
#define ISGIMPROPERSYSTEMFORCE_H

#include "iSGMTorsionSystemForce.h"
#include "iSGImproperSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"

namespace ProtoMol {
  //_________________________________________________________________ iSGImproperSystemForce

  template<class TBoundaryConditions>
  class iSGImproperSystemForce: public iSGMTorsionSystemForce<TBoundaryConditions>, private iSGImproperSystemForceBase {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

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
    virtual Force* doMake(std::string&, std::vector<Value>) const { return (new iSGImproperSystemForce()); }

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
  inline void iSGImproperSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
								    const Vector3DBlock* positions, 
								    Vector3DBlock* forces,
								    ScalarStructure* energies)
  {    
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;
    for (unsigned int i = 0; i < topo->impropers.size(); i++)
      calcTorsion(boundary,topo->impropers[i], positions, forces,
      (*energies)[ScalarStructure::IMPROPER],
      (*energies)[ScalarStructure::IMPROPER_DELTAMU],
      energies);
  }
  template<class TBoundaryConditions> 
  inline void iSGImproperSystemForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
									    const Vector3DBlock* positions, 
                                                                            Vector3DBlock* forces,
                                                                            ScalarStructure* energies)
  {
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;

    unsigned int n = topo->impropers.size();
    unsigned int count = numberOfBlocks(topo,positions);
  
    for(unsigned int i = 0;i<count;i++){
      if(Parallel::next()){
	int to = (n*(i+1))/count;
	if(to > static_cast<int>(n))
	  to = n;
	int from = (n*i)/count;
	for (int j = from; j < to; j++)
	  calcTorsion(boundary, topo->impropers[j], positions, forces,
          (*energies)[ScalarStructure::IMPROPER],
          (*energies)[ScalarStructure::IMPROPER_DELTAMU],
          energies);
      }
    }
  }

  template<class TBoundaryConditions> 
  inline unsigned int iSGImproperSystemForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology* topo,
									          const Vector3DBlock*){    
    return std::min(Parallel::getAvailableNum(),static_cast<int>(topo->impropers.size()));
  }
}
#endif /* ISGIMPROPERSYSTEMFORCE_H */
