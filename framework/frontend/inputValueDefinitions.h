/*  -*- c++ -*-  */
#ifndef INPUTVALUEDEFINITIONS_H
#define INPUTVALUEDEFINITIONS_H

#include "InputValue.h"

namespace ProtoMol
{
	/*
	 * It's good practice to only add keywords which all
	 * applications, especially protomol, benefit. If you
	 * need to add new keywords which only make sense
	 * in your application add them in front of main().
	 * 
	 */

	//________________________________________________________ InputValue<>
	declareInputValue( InputAMBER,STRING,NOTEMPTY);
	declareInputValue( InputBoundaryConditions,STRING,NOTEMPTY);
	declareInputValue( InputCellManager,STRING,NOTEMPTY);
	declareInputValue( InputCheckpointFreq, INT, NOTNEGATIVE);
	declareInputValue( InputCheckpointFile, STRING, NOTEMPTY);
	declareInputValue( InputConfig,STRING,NOTEMPTY);
	declareInputValue( InputCoords,STRING,NOTEMPTY);
	declareInputValue( InputDebug,INT,NOCONSTRAINTS);
	declareInputValue( InputDihedralMultPSF,BOOL,NOCONSTRAINTS);
	declareInputValue( InputDoSCPISM, BOOL, NOCONSTRAINTS);
	declareInputValue( InputEigTextFile,STRING,NOTEMPTY);
	declareInputValue( InputEigenValues,STRING,NOTEMPTY);
	declareInputValue( InputEigenVectors,STRING,NOTEMPTY);
	declareInputValue( InputFirststep,INT,NOCONSTRAINTS);
	declareInputValue( InputIntegrator,INTEGRATOR,NOTEMPTY);
	declareInputValue( InputMaxPackages,INT,NOCONSTRAINTS);
	declareInputValue( InputMinimalImage,BOOL,NOCONSTRAINTS);
	declareInputValue( InputMolVirialCalc,BOOL,NOCONSTRAINTS);
	declareInputValue( InputNumsteps,INT,NOTNEGATIVE);
	declareInputValue( InputOutput,BOOL,NOCONSTRAINTS);
	declareInputValue( InputOutputfreq,INT,NOTNEGATIVE);
	declareInputValue( InputPAR,STRING,NOTEMPTY);
	declareInputValue( InputPDBScaling,BOOL,NOCONSTRAINTS);
	declareInputValue( InputPSF,STRING,NOTEMPTY);
	declareInputValue( InputParallelMode,STRING,NOTEMPTY);
	declareInputValue( InputParallelPipe,INT,NOCONSTRAINTS);
	declareInputValue( InputPositions,STRING,NOTEMPTY);
	declareInputValue( InputRattle,BOOL,NOCONSTRAINTS);
	declareInputValue( InputRattleEpsilon,REAL,NOTNEGATIVE);
	declareInputValue( InputRattleMaxIter,INT,NOTNEGATIVE);
	declareInputValue( InputReducedImage,BOOL,NOCONSTRAINTS);
	declareInputValue( InputRemoveAngularMomentum,INT,NOCONSTRAINTS);
	declareInputValue( InputRemoveLinearMomentum,INT,NOCONSTRAINTS);
	declareInputValue( InputRestoreState, STRING, NOTEMPTY );
	declareInputValue( InputSeed,INT,NOTNEGATIVE);
	declareInputValue( InputShadow, BOOL, NOCONSTRAINTS );
	declareInputValue( InputShadowFreq, INT, NOTNEGATIVE );
	declareInputValue( InputShadowOrder, INT, NOTNEGATIVE );
	declareInputValue( InputShake,BOOL,NOCONSTRAINTS);
	declareInputValue( InputShakeEpsilon,REAL,NOTNEGATIVE);
	declareInputValue( InputShakeMaxIter,INT,NOTNEGATIVE);
	declareInputValue( InputTemperature,REAL,NOTNEGATIVE);
	declareInputValue( InputUseBarrier,BOOL,NOCONSTRAINTS);
	declareInputValue( InputVelocities,STRING,NOTEMPTY);
	declareInputValue( InputVirialCalc,BOOL,NOCONSTRAINTS);
}

#endif
