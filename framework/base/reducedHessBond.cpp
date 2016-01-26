#include "reducedHessBond.h"

namespace ProtoMol {
  Matrix3by3 reducedHessBond(const Vector3D& atom_i,
			     const Vector3D& atom_j,
			     const Real _k,
			     const Real _r0)
  {
    Vector3D rij(atom_j - atom_i);
    // the relative position vectors  

    Real a =  rij.normSquared();
    Real na =  sqrt(a);
  
    Real tm1 = 2.0 * _k * (na - _r0) / na;

    Real tm2 = 2.0 * _k * _r0 / na;

    return (Matrix3by3(tm1,0,0,0,tm1,0,0,0,tm1)+Matrix3by3(rij,rij*tm2/a));
  }
}
