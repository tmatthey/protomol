/*  -*- c++ -*-  */
#ifndef STAGE_H
#define STAGE_H

#include <string>
#include <vector>

#include "Real.h"
#include "Report.h"
#include "Array.h"
#include "Molecule.h"
#include "PSF.h"
#include "Vector3DBlock.h"

namespace ProtoMol {
  //_________________________________________________________________STAGE
  class STAGE{
  public:
    //
    // Container class for the detailed transformation path for each type of
    // molecular species.  This object contains the list of atomic charges
    // for each stage of a particular insertion/deletion and the stage numbers for which
    // each atom is to be transformed.

    //______________________________________________________________________AtomType
    struct AtomType{
      // This structure holds data for an atom type.  It consists of the mass, charge,
      // and alphaLJ parameter for the atom type, as well as the specific charge the
      // atomtype should have at each stage of the transformation.

      // constructors
      AtomType(const AtomType &AT) {
        type_name = AT.type_name;
        mass = AT.mass;
        charge = AT.charge;
        stage = AT.stage;
        alphaLJ = AT.alphaLJ;

        //unsigned int myCharges = AT.insert_charge.size();
        //insert_charge.resize(myCharges);
        //delete_charge.resize(myCharges);
        for (unsigned int i=0; i<AT.insert_charge.size(); i++) {
          insert_charge.push_back(AT.insert_charge[i]);
          delete_charge.push_back(AT.delete_charge[i]);
        }
      }
      AtomType(std::string N,
               Real M,
               Real Q,
               unsigned int S,
               unsigned int NS): type_name(N), mass(M), charge(Q), stage(S), alphaLJ(0.5) {
        insert_charge.resize(NS);
        delete_charge.resize(NS);
      }
      AtomType(std::string N,
               Real M,
               Real Q,
               unsigned int S,
               Real A,
               unsigned int NS): type_name(N), mass(M), charge(Q), stage(S), alphaLJ(A) {
        insert_charge.resize(NS);
        delete_charge.resize(NS);
      }

      /// the type name of this atomtype
      std::string type_name;
      /// permanent mass for this atomtype
      Real mass;
      /// permanent charge for this atomtype
      Real charge;
      /// the exact stage # in which this atomtype´s LJ interactions are switched on/off
      unsigned int stage;
      /// the value of the soft-core alphaLJ parameter to use
      Real alphaLJ;
      /// the charge this atomtype should have for each stage of an insertion/deletion
      std::vector<Real> delete_charge;
      std::vector<Real> insert_charge;

    };

    //______________________________________________________________________osmoticType
    struct osmoticType{
      // This structure holds the data for the atom types belonging to a particular osmotic
      // molecule.  It consists of AtomTypes, where AtomType is defined above, a molecule type
      // name, the number of stages to use for insertion/deletion, and the target fugacity.

      // constructor
      osmoticType(std::string N,
                  Real F,
                  unsigned int S): fugacity(F), NumberOfStages(S) {myMolecule.name = N;}

      /// the target fugacity for this molecule type
      Real fugacity;

      /// the number of stages to use for insertion/deletion
      unsigned int NumberOfStages;

      /// the list of AtomTypes on this molecule
      std::vector<AtomType> atomTypes;

      /// the information for one molecule of this component
      Molecule myMolecule;

      /// PSF object to hold the structure of this component
      PSF myStructure;

      /// stock xyz coordinates for a molecule of this component
      Vector3DBlock myCoordinates;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class STAGE
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // STAGE container
    std::vector<osmoticType> components;  //  MoleculeType container

  };

  //____________________________________________________________________________INLINES

}
#endif /* STAGE_H */
