//Update coordinates for GUI

#ifdef FAHCORE
//Don't want to include cosm.h here, it's painful to include 
//the cosm Makesystem with ours on windows and linux
//This will only work with 32 bit machines, if compiling on
//64 bit machine, be careful.
typedef unsigned char u8;
typedef float f32;
typedef unsigned int u32;
typedef signed long long int s64;
typedef unsigned long long int u64;
#include "guisrv.h"
#endif

static int gui_current( FAH_CURRENT *current, int step, int nsteps, float energy, float temp ) {
	current->iterations_done = step;
	current->energy = energy;
	current->temperature = temp;
	return 1;
}

static int gui_coords( FAH_XYZ *atoms, Vector3DBlock *pos ) {
	for (unsigned int i = 0; i < positions.size(); i++ ) {
		atoms[i].x = (*pos)[i].x;
		atoms[i].y = (*pos)[i].y;
		atoms[i].z = (*pos)[i].z;
	}
	//memcpy( atoms, x, sizeof(float) * 3 * natoms );
	return 1;
}

static int gui_bonds( FAH_BOND *bonds, GenericTopology *top ) {
	int i;

	for ( i = 0; i < top->bonds.size; i++ ) {
		//FAH requires a < b
		if (top->bonds[i].atom1 < top->bonds[i].atom2)
		{
			bonds[i].a = top->bonds[i].atom1;
			bonds[i].b = top->bonds[i].atom2;
		}else
		{
			bonds[i].a = top->bonds[i].atom2;
			bonds[i].b = top->bonds[i].atom1;
		}
	}
	return 1;
}

static int gui_atoms( FAH_ATOM *atoms, GenericTopology *top ) {
	int i;
	int element;
	//int atomtype;
	u8 myType[4];
	//float c6, c12;
	float mass;
	float myRadius;

	for ( i = 0; i < (top->atoms).size; i++ ) {
		//Determine element by mass
		//This will not work with united atom models or with heavy hydrogens etc
		element = 0; //unknown
		mass = top->atoms[i].scaledMass;
		if ( mass < 1.2 && mass >= 1.0 ) //hydrogen
			element = 1;
		    myType = "HXX";
		    myRadius = 1.2;
		else if ( mass > 11.8 && mass < 12.2 ) //carbon
			element = 6;
		    myType = "CXX";
		    myRadius = 1.7;
		else if ( mass > 14.0 && mass < 15 )   //nitrogen
			element = 7;
		    myType = "NXX";
		    myRadius = 1.55;
		else if ( mass > 15.5 && mass < 16.5 ) //oxygen
			element = 8;
		    myType = "OXX";
		    myRadius = 1.52;
		else if ( mass > 31.5 && mass < 32.5 ) //sulphur
			element = 16;
		    myType = "SXX";
		    myRadius = 1.85;
		else if ( mass > 29.5 && mass < 30.5 ) //phosphorus
			element = 15;
		    myType = "PXX";
		    myRadius = 1.9;

		atoms[i].type   = element;
		//strncpy( atoms[i].type, *(top->atoms.atomname[i]), 3 );
		atoms[i].charge = top->atoms[i].scaledCharge;

		//Convert c6, c12 to sigma for radius
		//atomtype = top->atoms.atom[i].type;
		//c6  = nbfp[ 2 * ntypes * atomtype + 2 * atomtype ];
		//c12 = nbfp[ 2 * ntypes * atomtype + 2 * atomtype + 1];
		//if ( c6 > 1e-8 ) {
		atoms[i].radius = myRadius;
		//}
		//else
		//	atoms[i].radius = 0.08;

		//Scale down for visualization
		atoms[i].radius *= 0.5;

		//Scale to angstroms
		//atoms[i].radius *= 10.0;

		fprintf( stderr,  "returning charge=%f radius=%f\n", atoms[i].charge, atoms[i].radius );
	}
	return 1;
}


static int gui_callback( int request, gsFAHData *fahdata, va_list context )
{
		GenericTopology* guitopo;
		Vector3DBlock* guipos;
		//rvec *x;
		//int ntypes;
		//float *nbfp;
		int step;
		int nsteps;
		float energy;
		float temp;

	guitopo = va_arg( context, GenericTopology* ); //PRB protomol topology
	guipos = va_arg( context, Vector3DBlock* );    //PRB protomol positions
	step   = va_arg( context, int );
	nsteps = va_arg( context, int );
	energy = va_arg( context, double );
	temp   = va_arg( context, double );
		
	switch( request ) {
		case 1:
			//fprintf( stderr, "name=%s, natoms=%d\n", fahdata->info.name, fahdata->info.atom_count, fahdata->info.bond_count );
			//Everything else has already been set through gsInitData
			fahdata->info.iterations = nsteps;
			gui_atoms( fahdata->atoms, guitopo );
			gui_bonds( fahdata->bonds, guitopo );
			break;
		case 2:
			gui_current( &fahdata->current, step, nsteps, energy, temp );
			gui_coords( fahdata->xyz, guipos );
			break;
	}
	return 1;
}