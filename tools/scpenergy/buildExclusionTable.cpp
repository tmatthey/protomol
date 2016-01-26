#include "buildExclusionTable.h"
#include "GenericTopology.h"
#include "Report.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

    void buildExclusionTable(GenericTopology* topo, const ExclusionType& exclusionType) {

        if(!exclusionType.valid())
            report << error <<"[buildExclusionTable()] Exclusion type not defined/valid."<<endr;


        topo->exclude = exclusionType;

        //  Resize array.
	topo->exclusions.resize(topo->atoms.size());

        //  If exclusionType is equal to NONE, return.
        if( exclusionType == ExclusionType::NONE )
            return;


        const int numBonds     = topo->bonds.size(),
                  numAngles    = topo->angles.size(),
                  numDihedrals = topo->dihedrals.size();



        //  Add excluded bonds.
        for( int i = 0; i < numBonds; i++ )
	  {
            topo->exclusions.add( topo->bonds[i].atom1, topo->bonds[i].atom2,
                    EXCLUSION_FULL );
	  }

        if( exclusionType != ExclusionType::ONE2 ) {

            //  Add excluded angles.
            for( int i = 0; i < numAngles; i++ )
                topo->exclusions.add( topo->angles[i].atom1,
                        topo->angles[i].atom3, EXCLUSION_FULL );

            if( exclusionType != ExclusionType::ONE3 ) {
                //  Add excluded dihedrals.
                for( int i = 0; i < numDihedrals; i++ ) {

                    if( exclusionType == ExclusionType::ONE4 )
                        topo->exclusions.add( topo->dihedrals[i].atom1,
                                topo->dihedrals[i].atom4, EXCLUSION_FULL );
                    else
                        topo->exclusions.add( topo->dihedrals[i].atom1,
                                topo->dihedrals[i].atom4, EXCLUSION_MODIFIED );

                }

            }


        }

        topo->exclusions.optimize();

    }


}
