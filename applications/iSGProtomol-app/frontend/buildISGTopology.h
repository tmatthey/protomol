#include "ExclusionType.h"

namespace ProtoMol {

  class GenericTopology;
  class PSF;
  class iSGPAR;

  //_________________________________________________________________ buildISGTopology

  void buildISGTopology(GenericTopology* topo,const PSF& psf, const iSGPAR& par);

  void buildMoleculeTable(GenericTopology *topo);
  void buildExclusionTable(GenericTopology* topo, const ExclusionType& exclusionType);
  void buildMoleculeBondingLists(GenericTopology* topo);

}
