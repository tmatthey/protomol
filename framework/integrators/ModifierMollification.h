/*  -*- c++ -*-  */
#ifndef MODIFIERMOLLIFICATION_H
#define MODIFIERMOLLIFICATION_H

#include "Modifier.h"
#include "MOLLYIntegrator.h"
#include "Vector3DBlock.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierMollification
  class ModifierMollification : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierMollification(MOLLYIntegrator* i):Modifier(),myTheIntegrator(i){}
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->mollification();
    }
    virtual std::string doPrint()const{return std::string("Mollification");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    MOLLYIntegrator* myTheIntegrator;
  };


}
#endif /* MODIFIERMOLLIFICATION_H */
