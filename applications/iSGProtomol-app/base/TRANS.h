/*  -*- c++ -*-  */
#ifndef TRANS_H
#define TRANS_H

#include <string>
#include <vector>

#include "Real.h"
#include "Report.h"
#include "Array.h"

namespace ProtoMol {
  //_________________________________________________________________TRANS
  class TRANS{
  public:
    //
    // Container class for the detailed transformation path for each type of  
    // molecular identity change.  This object contains the list of atomic charges
    // for each stage of a particular identity transformation, the alphaLJ parameters
    // for each atomtype and transformation, the atomic masses and stage numbers of each
    // atom type's identity, and also the ideal gas stage chemical potential differences
    // for each stage of a particular transformation.

    //______________________________________________________________________AtomType
    struct AtomType{
      // This structure holds data for an atom type.  It consists of a vector containing
      // all possible masses, and a vector containing all possible charges, and a vector
      // containing the alphaLJ parameters for a particular transformation type.

      AtomType(int I, int S, std::string name){
	mass.resize(I);
        charge.resize(I);
        old_charge.resize(ArraySizes(I)(I)(S));
        new_charge.resize(ArraySizes(I)(I)(S));
        alphaLJ.resize(ArraySizes(I)(I));
        for (int o=0; o<I; o++) {
          for (int n=0; n<I; n++) {
            for (int s=0; s<S; s++) {
              old_charge[o][n][s] = 0.0;
              new_charge[o][n][s] = 0.0;
            }
            alphaLJ[o][n] = 0.5;
          }
        }    
        stage.resize(I);
	type_name = name;
      }
      
      AtomType(Real M, Real Q, int S) {
	mass.push_back(M);
        charge.push_back(Q);
        stage.push_back(S);
      }

      // mass, charge, alphaLJ, and atomtype name
      std::vector<Real> mass;
      std::vector<Real> charge;
      Array<Real,3> old_charge;
      Array<Real,3> new_charge;
      Array<Real,2> alphaLJ;
      std::vector<int> stage;
      std::string type_name;

    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class TRANS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // TRANS container
    std::vector<AtomType> atomTypes;  //  AtomType container
    Array<Real,3> DeltaMuIG;          //  ideal gas chemical potential differences.  Units: (kcal/mol)
    unsigned int NumberOfStages;

  };

  //____________________________________________________________________________INLINES

}
#endif /* MQF_H */
