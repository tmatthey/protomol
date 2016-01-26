/* -*- c++ -*- */
#ifndef EXTERNALMAGNETICFIELDEXTENDEDFORCE_H
#define EXTERNALMAGNETICFIELDEXTENDEDFORCE_H

#include "ExtendedForce.h"

namespace ProtoMol {
  class Vector3DBlock;
  class ScalarStructure;
  class GenericTopology;
  //_________________________________________________________________ ExternalMagneticFieldExtendedForce
  class ExternalMagneticFieldExtendedForce : public ExtendedForce {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ExternalMagneticFieldExtendedForce();
    ExternalMagneticFieldExtendedForce(const Vector3D& b);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ExternalMagneticFieldExtendedForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doEvaluate(const GenericTopology* topo,
		    const Vector3DBlock* positions, 
		    const Vector3DBlock* velocities,
		    Vector3DBlock* forces,
		    ScalarStructure* energies, int from, int to); 

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
    virtual std::string getKeyword() const{return keyword;}
    virtual unsigned int numberOfBlocks(const GenericTopology* topo, 
					const Vector3DBlock* pos);
  private:
    virtual Force* doMake(std::string&, std::vector<Value> values) const;
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 1;}
    virtual void getParameters(std::vector<Parameter>&) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3D myB;

  }; 

  //______________________________________________________________________ INLINES
}
#endif /* EXTERNALMAGNETICFIELDEXTENDEDFORCE_H */
