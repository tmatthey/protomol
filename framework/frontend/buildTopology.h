/*  -*- c++ -*-  */
#include "ExclusionType.h"

namespace ProtoMol {

  class GenericTopology;
  class PSF;
  class PAR;

  //_________________________________________________________________ buildTopology

  void buildTopology(GenericTopology* topo,const PSF& psf, const PAR& par, bool dihedralMultPSF);
  
  void buildMoleculeTable(GenericTopology *topo);
  void buildExclusionTable(GenericTopology* topo, const ExclusionType& exclusionType);

}
