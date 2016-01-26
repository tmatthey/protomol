#!/usr/bin/awk -f

#  --------------------------------------------------------------------------  #
#  This script takes as input an allEnergies output file from ProtoMol and     #
#  prints out the average, std dev, max, min, and range.  The script will      #
#  ignore any lines that begin with a '#'.  The number of lines averaged is    #
#  listed in the top left corner.                                              #
#                                                                              #
#  usage: ./aveEnergies.awk allEnergiesFile                                    #
#  --------------------------------------------------------------------------  #

BEGIN {

    #  ----------------------------------------------------------------------  #
    #  Variables to decide if we should store first values as min/max.         #
    #  ----------------------------------------------------------------------  #

    initMax = 1
    initMin = 1

    #  ----------------------------------------------------------------------  #
    #  Change column titles here.                                              #
    #  ----------------------------------------------------------------------  #

    split( "Time Potential Kinetic Total Temperature Bond Angle Dihedral \
            Improper VdW Coulomb Other Volume Shadow ? ? ? ?", titles )
    
}

#  --------------------------------------------------------------------------  #
#  For each non-commented line, store the values in each column.  Save the     #
#  first set of values as both min and max.                                    #
#  --------------------------------------------------------------------------  #

/^ *#/ { commentLines++ }

$1 !~ /^ *#/ {

    #  ----------------------------------------------------------------------  #
    #  Don't bother storing the first column which contains the step number.   #
    #  ----------------------------------------------------------------------  #

    for( i = 2; i <= NF; i++ ) {

        colSum[i] += $i

        colSumSq[i] += ( $i * $i )

        if( initMin == 1 ) {
            min[i] = $i
        }
        else if( $i < min[i] ) {
            min[i] = $i
        }

        if( initMax == 1 ) {
            max[i] = $i
        }
        else if( $i > max[i] ) {
            max[i] = $i
        }

    }

    initMin = 0
    initMax = 0

}

END {

    #  ----------------------------------------------------------------------  #
    #  NF = number of fields (columns)                                         #
    #  NR = number of records (rows)                                           #
    #  ----------------------------------------------------------------------  #

    numCount = NR - commentLines


    #  ----------------------------------------------------------------------  #
    #  Print out the averages, etc.                                            #
    #  ----------------------------------------------------------------------  #

    printf "\n(%05d) %8s %16s %16s %16s %16s %16s\n", \
           numCount, "Col", "Mean", "Std Dev", "Max", "Min", "Range"

    printf "%16s %16s %16s %16s %16s %16s\n", \
           "---", "----", "-------", "---", "---", "-----"

    for( i = 2; i <= NF; i++ ) {

        mean = ( colSum[i] / numCount )

        #  ------------------------------------------------------------------  #
        #  Due to roundoff, it is possible to have a negative variance.  In    #
        #  this case, the std dev should be 0.  Otherwise, take the sqrt().    #
        #  ------------------------------------------------------------------  #

        variance = ( colSumSq[i] / numCount ) - ( mean * mean )

        stdDev = 0

        if( variance > 0 )
            stdDev = sqrt( variance )

        printf "%16s %16.6f %16.6f %16.6f %16.6f %16.6f\n", \
               titles[i], mean, stdDev, max[i], min[i], max[i]-min[i]

    }

}

