/* -*- c++ -*- */
#ifndef IMPROPERSYSTEMFORCE_H
#define IMPROPERSYSTEMFORCE_H

#include "MTorsionSystemForce.h"
#include "ImproperSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"

namespace ProtoMol {
  //_________________________________________________________________ ImproperSystemForce

  template<class TBoundaryConditions>
  class ImproperSystemForce: public MTorsionSystemForce<TBoundaryConditions>, private ImproperSystemForceBase {

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
    virtual Force* doMake(std::string&, std::vector<Value>) const { return (new ImproperSystemForce()); }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 0;}
    virtual void getParameters(std::vector<Parameter>&) const {}
  private:
    virtual void doSetParameters(std::string&, std::vector<Value>){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  }; 

  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions> 
  inline void ImproperSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
								 const Vector3DBlock* positions, 
								 Vector3DBlock* forces,
								 ScalarStructure* energies)
  {    
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;
    for (unsigned int i = 0; i < topo->impropers.size(); i++)
      calcTorsion(boundary,topo->impropers[i], positions, forces, (*energies)[ScalarStructure::IMPROPER],energies);
  }
  template<class TBoundaryConditions> 
  inline void ImproperSystemForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
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
	  calcTorsion(boundary, topo->impropers[j], positions, forces, (*energies)[ScalarStructure::IMPROPER],energies);
      }
    }
  }

  template<class TBoundaryConditions> 
  inline unsigned int ImproperSystemForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology* topo, 
									       const Vector3DBlock*){    
    return std::min(Parallel::getAvailableNum(),static_cast<int>(topo->impropers.size()));
  }
}
#endif /* IMPROPERSYSTEMFORCE_H */
