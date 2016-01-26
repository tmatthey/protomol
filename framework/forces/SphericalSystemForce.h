/* -*- c++ -*- */
#ifndef SPHERICALSYSTEMFORCE_H
#define SPHERICALSYSTEMFORCE_H

#include "SystemForce.h"
#include "Vector3DBlock.h"

namespace ProtoMol {
  //_________________________________________________________________ SphericalSystemForce
  class SphericalSystemForce : public SystemForce {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SphericalSystemForce();
    SphericalSystemForce(Vector3D center, Real radius, Real k, int j/*, Real u*/);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class SphericalSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doEvaluate(const GenericTopology* topo,
		    const Vector3DBlock* positions, 
		    Vector3DBlock* forces,
		    ScalarStructure* energies, int from, int to);

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
    virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 4;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3D myCenter;
    Real myRadius;
    Real myK;
    int  myJ;
    //    Real myU;
    Real myRadius2;
  }; 

  //______________________________________________________________________ INLINES
  inline void SphericalSystemForce::evaluate(const GenericTopology* topo,
					     const Vector3DBlock* positions, 
					     Vector3DBlock* forces,
					     ScalarStructure* energies){
    doEvaluate(topo,positions,forces,energies,0,(int)positions->size());
  }
}
#endif /* SPHERICALSYSTEMFORCE_H */
