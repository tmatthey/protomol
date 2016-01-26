#include "HMCIntegrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"

using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
    //____________________________________________________________ HMCIntegrator

    const string HMCIntegrator::keyword("HybridMC");

    const int HMCIntegrator::myNumParameters(2);

    //  ---------------------------------------------------------------------  //

    HMCIntegrator::HMCIntegrator()    
            : MCIntegrator() {}

    //  ---------------------------------------------------------------------  //

    HMCIntegrator::HMCIntegrator(int cycles,
            Real initialTemperature,
            ForceGroup *overloadedForces,
            StandardIntegrator *nextIntegrator)

            : MCIntegrator( cycles, initialTemperature, overloadedForces,
                    nextIntegrator) {
            }

    //  ---------------------------------------------------------------------  //

    void HMCIntegrator::initialize(GenericTopology *topo,
            Vector3DBlock *positions,
            Vector3DBlock *velocities,
            ScalarStructure *energies){
        MCIntegrator::initialize(topo,positions,velocities,energies);
    }

    //  ---------------------------------------------------------------------  //

    MTSIntegrator*  HMCIntegrator::doMake(string& , const vector<Value>& values, 
            ForceGroup* fg, StandardIntegrator *nextIntegrator)const{

        return new HMCIntegrator(values[0],values[1],fg,nextIntegrator);

    }

    //  ---------------------------------------------------------------------  //

    void HMCIntegrator::perturbSystem() {
        randomVelocity( getInitialTemperature(), myTopo, myVelocities );
        buildMolecularMomentum(myVelocities,myTopo);
    }

    //  ---------------------------------------------------------------------  //

    void HMCIntegrator::run( int numTimesteps ) {

        Real probability = 0.,
             MDinitEner  = 0.,
             MDfinEner   = 0.;

        report.precision(8);
        report.setf( std::ios::fixed );

        for( int i = 0; i < numTimesteps; i++ ) {

            saveValues();

            preStepModify();

            perturbSystem();

            //  The metropolis test is based in part on the change in KE from the new
            //  velocities.
            MDinitEner = myEnergies->potentialEnergy() +
                         kineticEnergy( myTopo, myVelocities );

            walk( myCycleLength ); 

            MDfinEner = myEnergies->potentialEnergy() +
                        kineticEnergy( myTopo, myVelocities );

            if( metropolisTest( MDfinEner, MDinitEner, probability ) ) {

                report << debug(1) << "Positions accepted. Difference: ";
                report.width(11);
                report << MDfinEner - MDinitEner << " Probability: ";
                report.width(11);
                report << probability << endr;

            }
            else {

                restoreValues();

                report << debug(1) << "Positions rejected. Difference: ";
                report.width(11);
                report << MDfinEner - MDinitEner << " Probability: ";
                report.width(11);
                report << probability << endr;

            } 

            postStepModify();

        }  

    }

}

