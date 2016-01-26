#!/usr/bin/gawk -f

BEGIN { count = 0; }

/ / { ValArray[1] += $4; ValSqArray[1] += $4 * $4; count++; }

END {

    #print "DihEnergy \tValue      \tValueStDev     ";
    #print "----------\t-----------\t---------------";

    aveVal = ( ValArray[1] / count );
    aveValStDev = sqrt( ValSqArray[1] / count - aveVal * aveVal );

    print "DihE", "\t", aveVal, "\t", aveValStDev;

}
