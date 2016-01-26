/*  -*- c++ -*-  */
#ifndef MODIFIERFORCEPROJECTION_H
#define MODIFIERFORCEPROJECTION_H

#include "Modifier.h"
#include "NormalModeUtilities.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierForceProjection
  class ModifierForceProjection : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierForceProjection(NormalModeUtilities* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->forceProjection();
    }
    virtual std::string doPrint()const{return std::string("ForceProjection");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    NormalModeUtilities* myTheIntegrator;
  };

}
#endif /* MODIFIER_H */
