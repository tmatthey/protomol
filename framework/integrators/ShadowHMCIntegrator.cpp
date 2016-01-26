#include "cmath"
#include "ScalarStructure.h"
#include "ShadowHMCIntegrator.h"
#include "STSIntegrator.h"

using namespace ProtoMol::Report;

using std::vector;
using std::string;
using std::deque;

namespace ProtoMol {

    //  ____________________________________________________ ShadowHMCIntegrator

    const string ShadowHMCIntegrator::keyword( "ShadowHMC" );


    //  --------------------------------------------------------------------  //

    ShadowHMCIntegrator::ShadowHMCIntegrator() : MCIntegrator(), myOrder(0),
                                                 myShadowK(0), myC(0),
                                                 savedPE(0), shadowMod(NULL),
                                                 myOptimize(false)
                                                 {}


    //  --------------------------------------------------------------------  //

    ShadowHMCIntegrator::ShadowHMCIntegrator( int cycleLen, Real initTemp,
                                              int shadowOrder, Real c,
                                              bool optimize,
                                              ForceGroup *overloadedForces,
                                              StandardIntegrator *nextIntegr )
            : MCIntegrator( cycleLen, initTemp, overloadedForces, nextIntegr ),
              myOrder(shadowOrder), myShadowK(shadowOrder/2), myC(c),
              savedPE(0), shadowMod(NULL), myOptimize(optimize) //, myRatio(ratio)
              {}


    //  --------------------------------------------------------------------  //

    ShadowHMCIntegrator::~ShadowHMCIntegrator() {}


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::initialize( GenericTopology *topo,
                                          Vector3DBlock   *positions,
                                          Vector3DBlock   *velocities,
                                          ScalarStructure *energies ) {

        MCIntegrator::initialize( topo, positions, velocities, energies );


        //  ----------------------------------------------------------------  //
        //  Remove the shadow modifier created if the keyword "shadowEnergy"  //
        //  was given in the config file.  We create our own modifier later.  //
        //  ----------------------------------------------------------------  //

        while( removeModifier( std::string( "Modifier Shadow" ) ) )
            ;  //  report << plain << "Removed modifier!" << endr;


        //  ----------------------------------------------------------------  //
        //  Create a shadow modifier and attach it to the next integrator.    //
        //  ----------------------------------------------------------------  //

        shadowMod = new ModifierShadow( myOrder, 1, next() );

        if( shadowMod == NULL )
            report << error << "Could not create shadow modifier.\n";
        else
            next()->adoptPostStepModifier( shadowMod );


        //  ----------------------------------------------------------------  //
        //  If the first step isn't taken, then the shadow will have a value  //
        //  of 0, which does weird things to my reweighted values.  Just set  //
        //  it to be the total energy to begin with.                          //
        //  ----------------------------------------------------------------  //

        (*myEnergies)[ScalarStructure::SHADOW] = myEnergies->potentialEnergy() +
                                                 kineticEnergy( myTopo, myVelocities );


        //  ----------------------------------------------------------------  //
        //  If optimization is on, make sure a valid ratio is also given.     //
        //  ----------------------------------------------------------------  //

        if( myOptimize ) {

            report << debug(1) << "Optimizing c" << endr;

            optimizeC();

        }

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::run( int numTimesteps ) {

        Real initTotEnergy    = 0.,
             finTotEnergy     = 0.,
             initShadowEnergy = 0.,
             finShadowEnergy  = 0.,
             MDinitEner       = 0.,
             MDfinEner        = 0.,
             energyDiff       = 0.,
             probability      = 0.;

        bool acceptMomenta = false;

        report.precision(8);
        report.setf( std::ios::fixed );

        for( int i = 0; i < numTimesteps; i++ ) {

            preStepModify();

            //  Save current pos, vel, forces, and energies in case we reject.                        
            saveValues();

			// save last shadow in case step not accepted //
			Real LastOtherEnergy = (*myEnergies)[ScalarStructure::OTHER];


            while( !acceptMomenta ) {

                //  Calculate new random velocities.                            
                perturbSystem();

                //  Clear the history queues and reset beta.                    
                shadowMod->resetHistory();

                //  Run enough 'previous' steps of next integrator to be able to
                //  calculate the shadow at step '0'.
                runPreSteps( myShadowK + 1 );

                //  Save the initial total and shadow energies.
                initTotEnergy = myEnergies->potentialEnergy() +
                                kineticEnergy( myTopo, myVelocities );

                initShadowEnergy = (*myEnergies)[ScalarStructure::SHADOW];


                //  Accept these momenta with probability proportional to:
                //  min{ 1, exp^( -\beta( H_[2k] - c - H ) ) }. 
                energyDiff = initShadowEnergy - myC - initTotEnergy;

                //  Decide whether or not to accept new momenta.
                if( metropolisTest( energyDiff, 0., probability ) ) {

                    acceptMomenta = true;

                    report << debug(1) << "  Momenta accepted. Difference: ";
                    report.width(11);
                    report << energyDiff << " Probability: ";
                    report.width(11);
                    report << probability << endr;

                }
                else {

                    //  If rejected, restore saved values and try again.
                    restoreValues();

                    report << debug(1) << "  Momenta rejected. Difference: ";
                    report.width(11);
                    report << energyDiff << " Probability: ";
                    report.width(11);
                    report << probability << endr;

                }


            }

            //  If numTimesteps > 1, we want to have to reaccept the momenta.
            acceptMomenta = false;


            //  Choose max{ H_[2k] - c, H } to be the initial 'energy'.         
            MDinitEner = std::max( initShadowEnergy - myC, initTotEnergy );


            //  FIXME: Don't need to calculate shadow every step ...
            next()->run( myCycleLength );


            //  Make note of the final total and shadow energies.               
            finTotEnergy = myEnergies->potentialEnergy() +
                           kineticEnergy( myTopo, myVelocities );

            finShadowEnergy = (*myEnergies)[ScalarStructure::SHADOW];


            //  Choose max{ H_[2k] - c, H } to be the final 'energy'.           
            MDfinEner = std::max( finShadowEnergy - myC, finTotEnergy );

            //  Based on the difference in the max energies, decide to accept   
            //  the new positions or not.                                       
            if( metropolisTest( MDfinEner, MDinitEner, probability  ) ) {

				// save sampling Hamiltonian value
				(*myEnergies)[ScalarStructure::OTHER] = MDfinEner;

                report << debug(1) << "Positions accepted. Difference: ";
                report.width(11);
                report << MDfinEner - MDinitEner << " Probability: ";
                report.width(11);
                report << probability << endr;

            }
            else {

                restoreValues();

				// restore last sampling Hamiltonian value
				(*myEnergies)[ScalarStructure::OTHER] = LastOtherEnergy;

                report << debug(1) << "Positions rejected. Difference: ";
                report.width(11);
                report << MDfinEner - MDinitEner << " Probability: ";
                report.width(11);
                report << probability << endr;

            }

            postStepModify();


        }


    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::runPreSteps( int numTimesteps ) {

        //  Save all values for current timestep (probably time 0).
        ScalarStructure saveEnergies( *myEnergies );
        Vector3DBlock   savePositions( *myPositions );
        Vector3DBlock   saveVelocities( *myVelocities );

        vector< Vector3DBlock > savedForces;
        vector< Real > savedPE;

        //  For each integrator level, save these values.
        Integrator *nextInt = next();

        for( int i = 0; i < level(); i++ ) {

            savedForces.push_back( *( nextInt->getForces() ) );
            savedPE.push_back( nextInt->myPotEnergy );

            if( level() - i > 1 )
                nextInt = nextInt->next();

        }


        //  This is a total kludge, but seems to work for the time being.
        next()->backward();
        next()->run( numTimesteps );

        next()->forward();
        next()->run( numTimesteps );


        //  Save shadow before restoring energies to step 0.
        Real savedShadow = (*myEnergies)[ScalarStructure::SHADOW];


        //  Restore original values to be exactly as they were at step 0.
        myEnergies->intoAssign( saveEnergies );
        myVelocities->intoAssign( saveVelocities );
        myPositions->intoAssign( savePositions );
        (*myEnergies)[ScalarStructure::SHADOW] = savedShadow;


        //  For each integrator level, restore these values.
        nextInt = next();

        for( int i = 0; i < level(); i++ ) {

            nextInt->getForces()->intoAssign( savedForces[i] );
            nextInt->myPotEnergy = savedPE[i];
            if( level() - i > 1 )
                nextInt = nextInt->next();

        }

        nextInt->myBeta = 0.;

    }


    //  --------------------------------------------------------------------  //

     MTSIntegrator * ShadowHMCIntegrator::doMake( string &/*errMsg*/,
                                                 const vector<Value>& values,
                                                 ForceGroup* forceGrp,
                                                 StandardIntegrator *nextIntegr ) const {

        return new ShadowHMCIntegrator( values[0], values[1], values[2],
                                        values[3], values[4], 
                                        forceGrp, nextIntegr );

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::getParameters(vector< Parameter> &parameters) const {

        MCIntegrator::getParameters(parameters);

        parameters.push_back( Parameter( "order2k", Value( myOrder ), 8.0, Text("Shadow calculation order") ));
        parameters.push_back( Parameter( "c", Value( myC ), 0.0, Text("SHMC offset paramiter 'c'") ) );
        parameters.push_back( Parameter( "optimize", Value( myOptimize ), false, Text("Automatic 'c' calculation") ) );

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::perturbSystem() {

        randomVelocity( getInitialTemperature(), myTopo, myVelocities );

        buildMolecularMomentum( myVelocities, myTopo );

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::saveValues() {

        savedPE = next()->myPotEnergy;

        myOldForces.intoAssign( *( next()->getForces() ) );

        MCIntegrator::saveValues();

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::restoreValues() {

        next()->myPotEnergy = savedPE;

        next()->getForces()->intoAssign( myOldForces );

        MCIntegrator::restoreValues();

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::calcDelGstats( Real  muH, Real  varH,
                                             Real &muG, Real &varG  ) {

        Real sigmaH = sqrt( varH ),
             beta_kb = 1. / ( Constant::BOLTZMANN * getInitialTemperature() );

        Real f2 = erfc( ( myC + beta_kb * varH - muH ) ) *
                  exp( beta_kb * ( myC + 0.5 * beta_kb * varH - muH ) );

        Real f1 = erfc( ( muH - myC ) / ( M_SQRT2 * sigmaH ) ) / f2;

        Real f3 = M_SQRT2 * beta_kb * varH * sigmaH *
                  exp( -( muH - myC) * ( muH - myC ) / ( 2 * varH ) ) /
                  ( 0.5 * M_2_SQRTPI *
                    erf( ( muH - myC ) / ( M_SQRT2 * sigmaH ) ) - 1 - f2 );


        muG =  muH  - ( beta_kb * varH ) / ( 1 + f1 );

        varG = varH + ( beta_kb * beta_kb * varH * varH * f1 ) /
                      ( 1 + 2 * f1 + f1 * f1 ) + f3;

    }


    //  --------------------------------------------------------------------  //

	void ShadowHMCIntegrator::optimizeC(){ // Real ratio ) {

        //  ----------------------------------------------------------------  //
        //  FIXME:  These should be user options.                             //
        //  ----------------------------------------------------------------  //

        int numOptSteps = 100;


        Real muDelH,
             varDelH,
             beta_kb = 1. / ( Constant::BOLTZMANN * getInitialTemperature() );


        //  ----------------------------------------------------------------  //
        //  First, calcuate the timestep according to the desired ratio of    //
        //  the variances.  cf. Eq.(xx) from SwHI0x.  For stability reasons,  //
        //  limit the maximum timestep to be 1.5 fs.                          //
        //  ----------------------------------------------------------------  //
		//  NOW IGNORED AS VARIANCE NOT INCREASED BY SHMC

        //calcMGStepVar( muDelH, varDelH, numOptSteps );

        //newTimestep = next()->getTimestep() *
        //              sqrt( sqrt( log( ratio ) ) / beta_kb / sqrt( varDelH ) );

        //if( newTimestep > 1.5 )
        //    newTimestep = 1.5;

        //dynamic_cast< STSIntegrator * >( next() )->setTimestep( newTimestep );


        //  ----------------------------------------------------------------  //
        //  FIXME:  Debugging only.                                           //
        //  ----------------------------------------------------------------  //

        /*
        report.precision(14);
        report.setf( std::ios::fixed );

        report << debug(1) << " mu( delta H ): ";
        report.width(18);
        report << muDelH << endr;

        report << debug(1) << "var( delta H ): ";
        report.width(18);
        report << varDelH << endr;

        report << debug(1) << "Sugg. timestep: ";
        report.width(18);
        report << newTimestep << endr;
        report << debug(1) << endr;
        */


        //  ----------------------------------------------------------------  //
        //  Now, run again at the new timestep to get the new mu and sigma.   //
        //  ----------------------------------------------------------------  //

        calcMGStepVar( muDelH, varDelH, numOptSteps );


        //  ----------------------------------------------------------------  //
        //  FIXME:  Debugging only.                                           //
        //  ----------------------------------------------------------------  //

        /*
        report << debug(1) << " mu( delta H ): ";
        report.width(18);
        report << muDelH << endr;

        report << debug(1) << "var( delta H ): ";
        report.width(18);
        report << varDelH << endr;
        report << debug(1) << endr;
        */


        Real muDelG          = 0.,
             varDelG         = 0.,
             mgRatio         = 0.,
             mdRatio         = 0.,
             hPointRatio     = 0.,
             sigmaDelH       = sqrt( varDelH ),
             aveNumSteps     = 10e9,  //  Arbitrarily large number.
             minNumSteps     = 10e9,  //  Arbitrarily large number.
             suggestedC      = 0,
             mgStageCost     = 2 * myShadowK + 2,    //  FIXME
             mdStageCost     = myCycleLength;


        //  ----------------------------------------------------------------  //
        //  Run several test cases to try out different c's.
        //  ----------------------------------------------------------------  //

        for( myC = muDelH - 4 * sigmaDelH;
             myC < muDelH - beta_kb * sigmaDelH;
             myC += 0.5 * sigmaDelH ) {


            //  ------------------------------------------------------------  //
            //  Given muDelH and varDelH, I can compute the same values of    //
            //  only the accepted momenta steps, muDelG and varDelG.          //
            //  ------------------------------------------------------------  //

            calcMGStepVar( muDelH, varDelH, numOptSteps );
            calcDelGstats( muDelH, varDelH, muDelG, varDelG );


            //  ------------------------------------------------------------  //
            //  Compute the expected ave number of MD steps per SHMC step.    //
            //  hPointRatio is an approx., P_HL, instead of the exact eq.     //
            //  ------------------------------------------------------------  //

            mgRatio = exp( beta_kb *
                           ( myC + 0.5 * beta_kb * varDelH - muDelH ) );

            hPointRatio = 0.5 * erfc( ( muDelG - myC ) /
                                      ( M_SQRT2 * sqrt( varDelG) ) );

            mdRatio = ( 1 - hPointRatio ) + 0.5 * hPointRatio * hPointRatio;

            aveNumSteps = ( mgStageCost / mgRatio ) +
                          ( mdStageCost / mdRatio );


            //  ------------------------------------------------------------  //
            //  Store the minimum average number of steps and it's c value.   //
            //  ------------------------------------------------------------  //

            if( aveNumSteps < minNumSteps ) {

                minNumSteps = aveNumSteps;
                suggestedC = myC;

            }


            //  ------------------------------------------------------------  //
            //  FIXME:  Debugging only.                                       //
            //  ------------------------------------------------------------  //

            /*
            report << debug(1) << endr;
            report << debug(1) << "           myC: ";
            report.width(18);
            report << myC << endr;

            report << debug(1) << "mgRatio       : ";
            report.width(18);
            report << mgRatio << endr;

            report << debug(1) << "mdRatio       : ";
            report.width(18);
            report << mdRatio << endr;

            report << debug(1) << " mu( delta H ): ";
            report.width(18);
            report << muDelH << endr;

            report << debug(1) << "var( delta H ): ";
            report.width(18);
            report << varDelH << endr;

            report << debug(1) << " mu( delta G ): ";
            report.width(18);
            report << muDelG << endr;

            report << debug(1) << "var( delta G ): ";
            report.width(18);
            report << varDelG << endr;

            report << debug(1) << "  Ave. # steps: ";
            report.width(18);
            report << aveNumSteps << endr;
            report << debug(1) << endr;
            */

        }


        myC = suggestedC;

    }


    //  --------------------------------------------------------------------  //

    void ShadowHMCIntegrator::calcMGStepVar( Real &mu, Real &var, int numSteps ) {

        Real totEnergy    = 0.,
             shadEnergy   = 0.,
             savedDeltaH  = 0.,
             savedDeltaH2 = 0.;

        for( int i = 0; i < numSteps; i++ ) {

            //  ------------------------------------------------------------  //
            //  Get new momenta; Run previous steps to calc H2k at time 0.    //
            //  ------------------------------------------------------------  //

            shadowMod->resetHistory();
            perturbSystem();
            runPreSteps( myShadowK + 1 );


            //  ------------------------------------------------------------  //
            //  Save delta H = H2k - H for calculating mu and sigma.          //
            //  ------------------------------------------------------------  //

            totEnergy = myEnergies->potentialEnergy() +
                        kineticEnergy( myTopo, myVelocities );

            shadEnergy = (*myEnergies)[ScalarStructure::SHADOW];

            savedDeltaH  += shadEnergy - totEnergy;

            savedDeltaH2 += ( shadEnergy - totEnergy ) *
                            ( shadEnergy - totEnergy );

        }

        mu  = savedDeltaH  / numSteps;
        var = savedDeltaH2 / numSteps - mu * mu;

    }


    //  --------------------------------------------------------------------  //

}

