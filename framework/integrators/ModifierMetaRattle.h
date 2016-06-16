/*  -*- c++ -*-  */
#ifndef MODIFIERMETARATTLE_H
#define MODIFIERMETARATTLE_H

#include "ModifierMetaRattleShake.h"
#include "Vector3DBlock.h"

namespace ProtoMol
{
	class Integrator;

	//_________________________________________________________________ ModifierMetaRattle
	class ModifierMetaRattle : public ModifierMetaRattleShake
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		ModifierMetaRattle(Real eps, int maxIter, int order);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class ModifierMetaRattleShake
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual Real calcError() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	};
}
#endif /* MODIFIERMETARATTLE_H */
