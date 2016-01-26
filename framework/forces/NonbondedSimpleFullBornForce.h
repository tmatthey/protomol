/* -*- c++ -*- */
#ifndef NONBONDEDSIMPLEFULLBORNFORCE_H
#define NONBONDEDSIMPLEFULLBORNFORCE_H

#include "SystemForce.h"
#include "NonbondedSimpleFullSystemForceBase.h"
#include "Parallel.h"
#include "mathutilities.h"
#include "Report.h"
using namespace ProtoMol::Report;
//#include "evaluateBorn.h"
namespace ProtoMol {
  //_________________________________________________________________ NonbondedSimpleFullBornForce

  template<class TOneAtomPair>
  class NonbondedSimpleFullBornForce: public SystemForce, private NonbondedSimpleFullSystemForceBase {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef SemiGenericTopology<typename TOneAtomPair::BoundaryConditions> RealTopologyType;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonbondedSimpleFullBornForce(): SystemForce(),myBlockSize(0),myCached(false){}; 
    NonbondedSimpleFullBornForce(TOneAtomPair oneAtomPair, unsigned int blockSize = defaultBlockSize) : 
      SystemForce(),myOneAtomPair(oneAtomPair),myBlockSize(blockSize),myCached(false){}
    virtual ~NonbondedSimpleFullBornForce(){}; 


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonbondedSimpleFullSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doEvaluate(const GenericTopology*, 
		    const Vector3DBlock*, 
		    Vector3DBlock*, 
		    ScalarStructure*, 
		    int i0, int i1, int j0, int j1);

    void doEvaluate(GenericTopology*, 
		    const Vector3DBlock*, 
		    Vector3DBlock*, 
		    ScalarStructure*, 
		    int i0, int i1, int j0, int j1);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology* topo, 
			  const Vector3DBlock* pos, 
			  Vector3DBlock* f, 
			  ScalarStructure* e){       
      if (!topo->doSCPISM)
	report << error << "To use SCPISM forces, please set doscpism in the configuration file" << endr;
      myCached = true;
      doEvaluate(topo,pos,f,e,0,topo->atoms.size(),0,topo->atoms.size());
    }

    virtual void parallelEvaluate(const GenericTopology* topo, 
				  const Vector3DBlock* pos, 
				  Vector3DBlock* f, 
				  ScalarStructure* e){
      if(!myCached)
	splitRangeArea(static_cast<unsigned int>(Parallel::getAvailableNum()),0,topo->atoms.size(),myFromRange,myToRange);
      myCached = true;

      for(int i=0;i<Parallel::getAvailableNum();i++)
	if(Parallel::next())
	  doEvaluate(topo,pos,f,e,myFromRange[i].first,myFromRange[i].second,myToRange[i].first,myToRange[i].second);
    }

    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual unsigned int numberOfBlocks(const GenericTopology*, const Vector3DBlock*){
      return Parallel::getAvailableNum();
    } 
    virtual std::string getKeyword() const{return keyword;}
    virtual void uncache(){myCached=false;};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:  
    virtual void getParameters(std::vector<Parameter>& parameters) const{
      myOneAtomPair.getParameters(parameters);
      parameters.push_back(Parameter("-blocksize",Value(myBlockSize,ConstraintValueType::Positive()),defaultBlockSize));
    }

    virtual unsigned int getParameterSize() const{return 1+TOneAtomPair::getParameterSize();}

    virtual std::string getIdNoAlias() const {
      return (TOneAtomPair::getId()+ " -algorithm " + keyword);
    }

  private:
    virtual Force* doMake(std::string& errMsg, std::vector<Value> values) const{
      unsigned int blockSize;
      int n = values.size()-1;
      values[n].get(blockSize);
      if(!(values[n].valid()) || blockSize == 0){
      	errMsg += keyword + " algorithm: 0 < blocksize (="+values[n].getString()+").";
      	return NULL;
      }
      return (new NonbondedSimpleFullBornForce(TOneAtomPair::make(errMsg,std::vector<Value>(values.begin(),values.end()-1)),
						blockSize));
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TOneAtomPair myOneAtomPair;
    unsigned int myBlockSize;
    std::vector<PairUInt> myFromRange;
    std::vector<PairUInt> myToRange;
    bool myCached;

  };
  //______________________________________________________________________ INLINES


  template<class TOneAtomPair>
  void NonbondedSimpleFullBornForce<TOneAtomPair>::doEvaluate(const GenericTopology* topo, 
								      const Vector3DBlock* positions, 
								      Vector3DBlock* forces, 
								      ScalarStructure* energies, 
								      int i0, int i1, int j0, int j1) {

    const RealTopologyType* realTopo = dynamic_cast<const RealTopologyType*>(topo);

    myOneAtomPair.initialize(realTopo, positions, forces, energies);
  
    for(int blocki = i0; blocki < i1; blocki += myBlockSize){
      int blocki_max = blocki;
      if(blocki_max < j0) blocki_max = j0;
      for(int blockj = blocki_max; blockj < j1; blockj += myBlockSize) {
	int istart = blocki;
	int iend = blocki + myBlockSize;
	if (iend > i1) iend = i1;
	for(int i = istart; i < iend; i++) {
	  int jstart = blockj;
	  if (jstart <= i) jstart = i+1;
	  int jend = blockj + myBlockSize;
	  if (jend > j1) jend = j1;
	  for(int j = jstart; j < jend; j++) {
	    myOneAtomPair.doOneAtomPair(i,j);
	  }
	}
      }
    }
  }

  template<class TOneAtomPair>
  void NonbondedSimpleFullBornForce<TOneAtomPair>::doEvaluate(GenericTopology* topo, 
								      const Vector3DBlock* positions, 
								      Vector3DBlock* forces, 
								      ScalarStructure* energies, 
								      int i0, int i1, int j0, int j1) {

    doEvaluate(const_cast<GenericTopology*>(static_cast<const GenericTopology*>(topo)),
	       positions, forces, energies);
    evaluateBorn(myOneAtomPair, topo, forces, energies);
  }

}
#endif /* NONBONDEDSIMPLEFULLSYSTEMFORCE_H */
