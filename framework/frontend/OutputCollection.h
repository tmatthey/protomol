/*  -*- c++ -*-  */
#ifndef OUTPUTCOLLECTION_H
#define OUTPUTCOLLECTION_H

#include <list>
#include "OutputCache.h"
#include "CheckpointInputStream.h"

namespace ProtoMol {
  class Output;
  class Configuration;
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class OutputFactoryDetails;
  class Integrator;
  //________________________________________________________ OutputCollection
  /**
     Container class for Output objects invoked at application level. Holds
     a cache object for reuse of values, which need to be computed from volatile
     values. The cache is cleared before invoking the first Output object. 
  */
  class OutputCollection  {
    friend class OutputFactoryDetails;
  private:
    typedef std::list<Output*> Container;
    typedef std::list<Output*>::iterator iterator;
  public:
    typedef std::list<Output*>::const_iterator const_iterator;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputCollection();
    ~OutputCollection();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputCollection
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void initialize(const Configuration* config, const Integrator* integrator, const GenericTopology* topo,
    		    const Vector3DBlock* pos, const Vector3DBlock*  vel, const ScalarStructure* energies);
    ///< Initialize all Output objects
    void run(int step);
    ///< Invoke all Output objects with run()
    void updateNextStep(int step);
    void finalize(int step);
    ///< Finalize all Outout objects
    int  getNext() const;
    void adoptOutput(Output* output); 
    ///< Add new Output object to the collection
    template<class T>
    void addToCache(const T& obj){myCache->add(obj);}
    ///< Add a structure or value to Cache

    /// Iterators, const
    const_iterator begin() const {return myOutputList.begin();}
    const_iterator end()   const {return myOutputList.end();} 
    void setRestore() {myCache->setRestore();}
    void clearRestore() {myCache->clearRestore();}
    void restoreState(CheckpointInputStream& is);

  private:
    /// Iterators
    iterator       begin()       {return myOutputList.begin();}
    iterator       end()         {return myOutputList.end();}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Container    myOutputList;
    OutputCache* myCache;

    const Configuration*   myConfig; 
    const GenericTopology* myTopology; 
    const Integrator*      myIntegrator;
    const ScalarStructure* myEnergies;
    const Vector3DBlock*   myPositions;
    const Vector3DBlock*   myVelocities; 
  };
}
#endif
