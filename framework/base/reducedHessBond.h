/*  -*- c++ -*-  */
#ifndef REDUCEDHESSBOND_H
#define REDUCEDHESSBOND_H

#include "Vector3D.h"
#include "Matrix3by3.h"

namespace ProtoMol
{
	Matrix3by3 reducedHessBond(const Vector3D& atom_i,
	                           const Vector3D& atom_j,
	                           const Real _k,
	                           const Real _r0);
}
#endif
