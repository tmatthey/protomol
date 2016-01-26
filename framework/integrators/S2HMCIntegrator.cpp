#include "S2HMCIntegrator.h"
#include "Report.h"
#include "mathutilities.h"
#include "pmconstants.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"
#include "topologyutilities.h"
#include "GenericTopology.h"
#include "StandardIntegrator.h"

using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
    //____________________________________________________________ S2HMCIntegrator

    const string S2HMCIntegrator::keyword( "S2HMCIntegrator" );

    const int S2HMCIntegrator::myNumParameters(2);

    //  --------------------------------------------------------------------  //

    S2HMCIntegrator::S2HMCIntegrator()    
            : MCIntegrator() {}

    //  --------------------------------------------------------------------  //

    S2HMCIntegrator::S2HMCIntegrator(int cycles,
            Real initialTemperature,
            ForceGroup *overloadedForces,
            StandardIntegrator *nextIntegrator)

            : MCIntegrator( cycles, initialTemperature, overloadedForces,
                    nextIntegrator) {
            }

    //  --------------------------------------------------------------------  //

    void S2HMCIntegrator::initialize(GenericTopology *topo,
            Vector3DBlock *positions,
            Vector3DBlock *velocities,
            ScalarStructure *energies){

        MCIntegrator::initialize(topo,positions,velocities,energies);

        //  ----------------------------------------------------------------  //
        //  Remove the shadow modifier created if the keyword "shadowEnergy"  //
        //  was given in the config file.  We use the 4th order 'seperable    //
        //  approximation, then the reported shadow can be used for re-       //
        //  weighting.                                                        //
        //  ----------------------------------------------------------------  //

        while( removeModifier( std::string( "Modifier Shadow" ) ) )

        //  ----------------------------------------------------------------  //
        //  Invert the masses to avoid multiple divisions.
        //  ----------------------------------------------------------------  //

        inverseMass.resize( myPositions->size() );

        for( unsigned int i = 0; i < inverseMass.size(); i++ )
            inverseMass[i] = 1. / myTopo->atoms[i].scaledMass;


    }

    //  --------------------------------------------------------------------  //

    MTSIntegrator*  S2HMCIntegrator::doMake(string& , const vector<Value>&
            values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const{

        return new S2HMCIntegrator(values[0],values[1],fg,nextIntegrator);

    }

    //  --------------------------------------------------------------------  //

    void S2HMCIntegrator::perturbSystem() {
        randomVelocity( getInitialTemperature(), myTopo, myVelocities );
        buildMolecularMomentum(myVelocities,myTopo);
    }

    //  --------------------------------------------------------------------  //

    void S2HMCIntegrator::run( int numTimesteps ) {


        //  ----------------------------------------------------------------  //
        //  Useful constants.                                                 //
        //  ----------------------------------------------------------------  //

        Real dt = next()->getTimestep() * Constant::INV_TIMEFACTOR;

        const Real dt_24  = dt / 24.,
                   dt2_24 = dt * dt_24;

        //  Fix time -------------------------------------------------------  //
        Real actTime;
        actTime = myTopo->time + myCycleLength * numTimesteps * next()->getTimestep();


        //  ----------------------------------------------------------------  //

        for( int i = 0; i < numTimesteps; i++ ) {

            saveValues();

            preStepModify();


            //  ------------------------------------------------------------  //
            //  Generate new set of random momenta, P:
            //  ------------------------------------------------------------  //

            perturbSystem();  //  TODO: Do I need to save these?


            //  ------------------------------------------------------------  //
            //  Compute shadow: S(q,p) = H(q,p) + dt^2/24 Uq M^-1 Uq 
            //  ------------------------------------------------------------  //

            Real initKE = kineticEnergy( myTopo, myVelocities );

            Real initShadow = myEnergies->potentialEnergy() + initKE;

            Vector3DBlock Uq( *( next()->getForces() ) );

            Real Uq_M_Uq = 0.;

            for( unsigned int i = 0; i < myPositions->size(); i++ )
                Uq_M_Uq += Uq[i].normSquared() * inverseMass[i];

            initShadow += dt2_24 * Uq_M_Uq;


            //  ------------------------------------------------------------  //
            //  Iteratively solve for \hat p, then compute \hat q.            //
            //  ------------------------------------------------------------  //

            solveForP();

            ((StandardIntegrator*)next())->calculateForces();

            //  ------------------------------------------------------------  //
            //  Run MD
            //  ------------------------------------------------------------  //

            walk( myCycleLength ); 


            //  ------------------------------------------------------------  //
            //  Iteratively solve for q, then compute p.                      //
            //  ------------------------------------------------------------  //

            solveForQ();

            ((StandardIntegrator*)next())->calculateForces();


            //  ------------------------------------------------------------  //
            //  Compute final shadow.
            //  ------------------------------------------------------------  //

            Real finShadow = myEnergies->potentialEnergy() +
                             kineticEnergy( myTopo, myVelocities );

            Uq.intoAssign( *( next()->getForces() ) );

            Uq_M_Uq = 0.;

            for( unsigned int i = 0; i < myPositions->size(); i++ )
                Uq_M_Uq += Uq[i].normSquared() * inverseMass[i];

            finShadow += dt2_24 * Uq_M_Uq;


            //  ------------------------------------------------------------  //
            //  Metropolis test.                                              //
            //  ------------------------------------------------------------  //

            report.precision(18);
            report.setf( std::ios::fixed );

            Real probability = 0.;

            if( metropolisTest( finShadow, initShadow, probability ) ) {

                (*myEnergies)[ScalarStructure::SHADOW] = finShadow;
                report << debug(1) << "Positions accepted. Difference: ";
                report.width(11);
                report << finShadow - initShadow << " Probability: ";
                report.width(11);
                report << probability << endr;

            }
            else {
                initShadow -= initKE;

                restoreValues();

                initShadow += kineticEnergy( myTopo, myVelocities );
                (*myEnergies)[ScalarStructure::SHADOW] = initShadow;

                report << debug(1) << "Positions rejected. Difference: ";
                report.width(11);
                report << finShadow - initShadow << " Probability: ";
                report.width(11);
                report << probability << endr;
            }
            
            postStepModify();

        }  
        //  Fix time -------------------------------------------------------  //
        myTopo->time = actTime;
        //

    }


    //  --------------------------------------------------------------------  //

    void S2HMCIntegrator::solveForP() {


        //  ----------------------------------------------------------------  //
        //  Useful constants.                                                 //
        //  ----------------------------------------------------------------  //

        const Real dt     = next()->getTimestep() * Constant::INV_TIMEFACTOR,
                   dt_24  = dt / 24.,
                   dt2_24 = dt * dt_24;


        //  ----------------------------------------------------------------  //
        //  Compute Uq( q - dt V )                                            //
        //  ----------------------------------------------------------------  //

        Vector3DBlock origPos( *myPositions );

        myPositions->intoWeightedAdd( -dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();

        Vector3DBlock newVel( *( next()->getForces() ) );


        //  ----------------------------------------------------------------  //
        //  Compute Uq( q + dt V )                                            //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( origPos );

        myPositions->intoWeightedAdd( dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();


        //  ----------------------------------------------------------------  //
        //  V = v - dt / ( 24 M ) * ( Uq( q + dt V ) - Uq( q - dt V ) )       //
        //  ----------------------------------------------------------------  //

        newVel.intoSubtract( *( next()->getForces() ) );

        for( unsigned int i = 0; i < newVel.size(); i++ )
            newVel[i] = (*myVelocities)[i] - newVel[i] * dt_24 * inverseMass[i];


        //  ----------------------------------------------------------------  //
        //  Iteratively solve for V using previous approximation.             //
        //  ----------------------------------------------------------------  //

        Real norm = 1;

        while ( norm > 10e-8 ) {

            Vector3DBlock prevNewVel( newVel );


            //  ------------------------------------------------------------  //
            //  Compute U_q( q - dt V )                                       //
            //  ------------------------------------------------------------  //

            myPositions->intoAssign( origPos );

            myPositions->intoWeightedAdd( -dt, prevNewVel );

            ((StandardIntegrator*)next())->calculateForces();

            newVel.intoAssign( *( next()->getForces() ) );


            //  ------------------------------------------------------------  //
            //  Compute U_q( q + dt V )                                       //
            //  ------------------------------------------------------------  //

            myPositions->intoAssign( origPos );

            myPositions->intoWeightedAdd( dt, prevNewVel );

            ((StandardIntegrator*)next())->calculateForces();


            //  ------------------------------------------------------------  //
            //  V = v - dt / ( 24 M ) * ( Uq( q + dt V ) - Uq( q - dt V ) )   //
            //  ------------------------------------------------------------  //

            newVel.intoSubtract( *( next()->getForces() ) );

            for( unsigned int i = 0; i < newVel.size(); i++ )
                newVel[i] = (*myVelocities)[i] - newVel[i] * dt_24 * inverseMass[i];


            //  ------------------------------------------------------------  //
            //  Compute norm.                                                 //
            //  ------------------------------------------------------------  //

            norm = 0.;

            for( unsigned int i = 0; i < prevNewVel.size(); i++ )
                norm += Vector3D( newVel[i] - prevNewVel[i] ).normSquared();

        }


        //  ----------------------------------------------------------------  //
        //  Now, solve for Q.                                                 // 
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( origPos );

        myPositions->intoWeightedAdd( dt, newVel );

        ((StandardIntegrator*)next())->calculateForces();

        Vector3DBlock newPos( *( next()->getForces() ) );


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( origPos );

        myPositions->intoWeightedAdd( -dt, newVel );

        ((StandardIntegrator*)next())->calculateForces();

        newPos.intoAdd( *( next()->getForces() ) );

        for( unsigned int i = 0; i < newPos.size(); i++ )
            newPos[i] = origPos[i] - newPos[i] * dt2_24 * inverseMass[i];


        //  ----------------------------------------------------------------  //
        //  Save new (q,p).                                                   //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( newPos );

        myVelocities->intoAssign( newVel );


    }


    //  --------------------------------------------------------------------  //

    void S2HMCIntegrator::solveForQ() {


        //  ----------------------------------------------------------------  //
        //  Useful constants.                                                 //
        //  ----------------------------------------------------------------  //

        const Real dt     = next()->getTimestep() * Constant::INV_TIMEFACTOR,
                   dt_24  = dt / 24.,
                   dt2_24 = dt_24 * dt;


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        Vector3DBlock origPos( *myPositions );

        myPositions->intoWeightedAdd( dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();

        Vector3DBlock newPos( *( next()->getForces() ) );


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( origPos );

        myPositions->intoWeightedAdd( -dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();

        newPos.intoAdd( *( next()->getForces() ) );

        for( unsigned int i = 0; i < newPos.size(); i++ )
            newPos[i] = origPos[i] + newPos[i] * dt2_24 * inverseMass[i];


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        Real norm = 1;

        while ( norm > 10e-8 ) {

            //  ------------------------------------------------------------  //
            //                                                                //
            //  ------------------------------------------------------------  //

            Vector3DBlock prevNewPos( newPos );


            //  ------------------------------------------------------------  //
            //  Compute U_q( q - dt V )                                       //
            //  ------------------------------------------------------------  //

            myPositions->intoAssign( prevNewPos );

            myPositions->intoWeightedAdd( dt, *myVelocities );

            ((StandardIntegrator*)next())->calculateForces();

            newPos.intoAssign( *( next()->getForces() ) );


            //  ------------------------------------------------------------  //
            //  Compute U_q( q + dt V )                                       //
            //  ------------------------------------------------------------  //

            myPositions->intoAssign( prevNewPos );

            myPositions->intoWeightedAdd( -dt, *myVelocities );

            ((StandardIntegrator*)next())->calculateForces();

            newPos.intoAdd( *( next()->getForces() ) );

            for( unsigned int i = 0; i < newPos.size(); i++ )
                newPos[i] = origPos[i] + newPos[i] * dt2_24 * inverseMass[i];


            //  ------------------------------------------------------------  //
            //  Compute new norm.                                             //
            //  ------------------------------------------------------------  //

            norm = 0.;

            for( unsigned int i = 0; i < newPos.size(); i++ )
                norm += Vector3D( newPos[i] - prevNewPos[i] ).normSquared();

        }


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( newPos );

        myPositions->intoWeightedAdd( -dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();

        Vector3DBlock newVel( *( next()->getForces() ) );


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( newPos );

        myPositions->intoWeightedAdd( dt, *myVelocities );

        ((StandardIntegrator*)next())->calculateForces();

        newVel.intoSubtract( *( next()->getForces() ) );

        for( unsigned int i = 0; i < newVel.size(); i++ )
            newVel[i] = (*myVelocities)[i] + newVel[i] * dt_24 * inverseMass[i];


        //  ----------------------------------------------------------------  //
        //  Save new (q,p).                                                   //
        //  ----------------------------------------------------------------  //

        myPositions->intoAssign( newPos );

        myVelocities->intoAssign( newVel );


    }

    
}

