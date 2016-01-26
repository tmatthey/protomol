/*  -*- c++ -*-  */
#ifndef STSINTEGRATOR_H
#define STSINTEGRATOR_H

#include "StandardIntegrator.h"
#include "GenericTopology.h"

namespace ProtoMol {

  class ScalarStructure;
  class ForceGroup;

  //_________________________________________________________________ STSIntegrator

  class STSIntegrator: public StandardIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    STSIntegrator();
    STSIntegrator(Real timestep, ForceGroup *overloadedForces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    STSIntegrator* make(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg)const;
    virtual Real setTimestep( Real );
  protected:
    virtual void doDrift();
  private:
    virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg)const=0;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void getParameters(std::vector<Parameter> &parameter) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(GenericTopology *topo, 
			    Vector3DBlock   *positions, 
			    Vector3DBlock   *velocities, 
			    ScalarStructure *energies);
    virtual Integrator* next(){return NULL;}
    virtual const Integrator* next() const {return NULL;}
    virtual Real getTimestep() const;
  protected:
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doDriftOrNextIntegrator();
    virtual void calculateForces();
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myTimestep;

  };
  //______________________________________________________________________ INLINES

  inline Real STSIntegrator::getTimestep() const {
    return ( isForward() ? myTimestep : -myTimestep );
  }

  inline Real STSIntegrator::setTimestep( Real newTimestep ) {

    Real oldTimestep = myTimestep;

    myTimestep = newTimestep;

    return ( oldTimestep );

  }

}
#endif
