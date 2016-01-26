#!/bin/bash -f
                                                                                     
# vim: set cindent:
                                                                                     
#  -------------------------------------------------------------------------  #
#  Define some local variables.
#  -------------------------------------------------------------------------  #
                                                                                     
mol="UA_butane"
date="Mar30"
outfreq="10"
compute="true"
print="false"
subscript="PM"
                                    
#  -------------------------------------------------------------------------  #
#  Determine the full path to the tests directory.
#  -------------------------------------------------------------------------  # 
basedir=`pwd`          
#  -------------------------------------------------------------------------  #
#  Define all test options here.  Each option will have it's own directory.
#  -------------------------------------------------------------------------  #

temperatures="300 500"
#temperatures="100 300 500 700"
#temperatures="700"

alg="hmc dhmc" 
#alg="HMC DHMC CDHMC"
#alg="CDHMC"
#alg="MD HMC DHMC CDHMC"

steps="1000000"    
                                  
#  -------------------------------------------------------------------------  #
#  Loop through each option.
#  -------------------------------------------------------------------------  #
                                                                                    
for temperature in $temperatures
  do

  for method in $alg
    do
        
  for nsteps in $steps
    do

      #  -------------------------------------------------------------  #
      #  Set the directory to work in.
      #  -------------------------------------------------------------  #
                                                                                     
      testdir="$basedir/$date/$temperature/$method"
                                                                   
      #  -------------------------------------------------------------  #
      #  Create the test directories and copy inputs if needed.
      #  -------------------------------------------------------------  #
                                                                                     
      if [ ! -d $testdir ];
          then
                                                                                     
          mkdir -p $testdir
                                                                                     
          echo "creating directory:  $testdir"
                                                                                     
	  cp $mol.pdb $testdir
	  cp $mol.psf $testdir
	  cp $mol.par $testdir
          cp dihedralSetDHMC $testdir
	  cp angleSetDHMC $testdir
                                                                           
      fi                                                 
                                                                                     
      #  -------------------------------------------------------------  #
      #  Copy the config files over.
      #  -------------------------------------------------------------  #

      if [ $method == "MD" ]; 
        then                                                                      
        sed -e "s/thetemperature/$temperature/g" $mol.$method.conf > $testdir/$mol.$method.conf
      fi
      if [ $method == "hmc" ];
        then
        sed -e "s/thetemperature/$temperature/g" $mol.$method.conf > $testdir/$mol.$method.conf
      fi
      if [ $method == "dhmc" ];
        then
        sed -e "s/thetemperature/$temperature/g" $mol.$method.conf > $testdir/$mol.$method.conf
      fi
                                                                               
      #  -------------------------------------------------------------  #
      #  Run the tests.
      #  -------------------------------------------------------------  #
                                                                                     
      cd $testdir

	if [ $compute != "false" ];
	    then
            echo "Running from directory:  $testdir"
	    scriptfile="$subscript$temperature$method.pbs"
	    touch $scriptfile
	    echo "#!/bin/bash" >> $scriptfile
	    echo "#PBS -N $scriptfile" >> $scriptfile
	    echo "#PBS -l nodes=1:ppn=1" >> $scriptfile
	    echo "cd \$PBS_O_WORKDIR" >> $scriptfile
            echo "mkdir /var/scratch/pbrenne1/$temperature$method" >> $scriptfile
            echo "mv * /var/scratch/pbrenne1/$temperature$method" >> $scriptfile
            echo "cd /var/scratch/pbrenne1/$temperature$method" >> $scriptfile
            echo "CBprotomol $mol.$method.conf --numsteps $nsteps --outputfreq $outfreq >& screenoutput$method.txt" >> $scriptfile
            echo "mv * \$PBS_O_WORKDIR" >> $scriptfile
	    echo "cd \$PBS_O_WORKDIR" >> $scriptfile
            echo "rm -rf /var/scratch/pbrenne1/$temperature$method" >> $scriptfile
	    chmod a+x $scriptfile
	    qsub $scriptfile
	fi                                                                                     

      #  -------------------------------------------------------------  #
      #  Print the test results.
      #  -------------------------------------------------------------  #
                                                   

	if [ $print != "false" ];
	    then

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
	fi
                                                              
      sleep 1
                                                                                     
      cd $basedir
  done
  done
done  #  End foreach().