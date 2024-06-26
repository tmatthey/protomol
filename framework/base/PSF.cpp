#include "PSF.h"

namespace ProtoMol
{
	//_____________________________________________________________________ PSF

	void PSF::clear()
	{
		atoms.clear();
		bonds.clear();
		angles.clear();
		dihedrals.clear();
		impropers.clear();
		donors.clear();
		acceptors.clear();
		nonbondeds.clear();
		ngrp.clear();
	}
}
