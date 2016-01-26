/*  -*- c++ -*-  */
#ifndef MODIFIERBERENDSEN_H
#define MODIFIERBERENDSON_H

#include "Modifier.h"
#include "BerendsenIntegrator.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierBerendsen
  class ModifierBerendsen : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierBerendsen(BerendsenIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->doBerendsen();
    }
    virtual std::string doPrint()const{return std::string("BerendsenMod");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    BerendsenIntegrator* myTheIntegrator;
  };

}
#endif /* MODIFIER_H */
