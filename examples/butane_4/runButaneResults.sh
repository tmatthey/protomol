#!/bin/bash -f
                                                                                     
#  -------------------------------------------------------------------------  #
#  Define some local variables.
#  -------------------------------------------------------------------------  #
                                                                                     
mol="UA_butane"
                                    
#  -------------------------------------------------------------------------  #
#  Determine the full path to the tests directory.
#  -------------------------------------------------------------------------  #

basedir=`pwd`
         
#  -------------------------------------------------------------------------  #
#  Define all graph label options here.
#  -------------------------------------------------------------------------  #

temperature="500"
alg="hmc" 
steps="100000"     
                                  
#  -------------------------------------------------------------  #
#  Process the dihedrals output file
#  -------------------------------------------------------------  #

  # Use awk script to get the dihedral energy average
  # and create a dihedral value data file for the xmgrace histogram

  gawk -f procDihedralFile.awk UA_butane.$alg.dihedrals

  rm -f dhmc.out.dihedrals

  gawk '{print $3}' UA_butane.$alg.dihedrals > $alg.out.dihedrals

  # Use xmgrace to plot the dihedral values in sequence
  # and plot a histogram showing the dihedral value distribution
	   
  mytitle="title \"Method $alg - Leap Frog 1.0fs - TEMP $temperature - $steps STEPS\""
	   
  xmgrace -hdevice PNG -printfile dihedrals.png $alg.out.dihedrals \
          -pexec autoscale -pexec "$mytitle"  \
	  -pexec 'WORLD YMAX 6.3' -pexec 'WORLD YMIN 0.0' &
          #-param ../MDgraphtemplate -This is useful once you choose a final format

  xmgrace -hdevice EPS -printfile dihedrals.eps $alg.out.dihedrals \
          -pexec "HISTOGRAM(S0, MESH(0,6.3,100), OFF, OFF)" \
	  -pexec "S1.y = S1.y / SUM(S1.y)" -pexec "KILL G0.S0" \
	  -pexec autoscale \
	  -pexec "$mytitle" -pexec 'xaxis label "Dihedral Value"' \
	  -pexec 'yaxis label "% of Dihedrals at Given Value"' &
          #-param ../MDgraphtemplate -This is useful once you choose a final format