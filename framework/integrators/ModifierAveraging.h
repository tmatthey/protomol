/*  -*- c++ -*-  */
#ifndef MODIFIERAVERAGING_H
#define MODIFIERAVERAGING_H

#include "Modifier.h"
#include "MOLLYIntegrator.h"
#include "Vector3DBlock.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierAveraging
  class ModifierAveraging : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierAveraging(MOLLYIntegrator* i):Modifier(),myTheIntegrator(i){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->averagingPositions();
    }
    virtual std::string doPrint()const{return std::string("Averaging");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    MOLLYIntegrator* myTheIntegrator;
  };

}
#endif /* MODIFIERAVERAGING_H */
