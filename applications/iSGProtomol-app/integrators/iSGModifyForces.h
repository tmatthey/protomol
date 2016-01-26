/*  -*- c++ -*-  */
#ifndef iSGMODIFYFORCES_H
#define iSGMODIFYFORCES_H

#include "Modifier.h"
#include "iSGIntegrator.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ iSGModifyForces
  class iSGModifyForces : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGModifyForces(iSGIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    iSGModifyForces(iSGIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      //report << hint << "modifyForces" << endr;
      myTheIntegrator->modifyForces();
    }
    virtual std::string doPrint()const{return std::string("iSGModifyForces");}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    iSGIntegrator* myTheIntegrator;
  };
} // end modifier
#endif /* iSGMODIFYFORCES_H */
