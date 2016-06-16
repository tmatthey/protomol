/*  -*- c++ -*-  */
#ifndef MODIFIERMETASHAKE_H
#define MODIFIERMETASHAKE_H

#include "ModifierMetaRattleShake.h"
#include "Vector3DBlock.h"

namespace ProtoMol
{
	class Integrator;

	//_________________________________________________________________ ModifierMetaShake
	class ModifierMetaShake : public ModifierMetaRattleShake
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierMetaShake(Real eps, int maxIter, int order);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class ModifierMetaShakeShake
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual Real calcError() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	};
}
#endif /* MODIFIERMETASHAKE_H */
