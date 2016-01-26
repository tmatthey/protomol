#! /bin/sh

failbit=0

# MINIMIZATION
minimization() {
    echo "******** RUNNING MINIMIZATION ********"
    if [ -f "Minimization/output.regression" ]; then # REMOVE FILE
       /bin/rm Minimization/output.regression
    fi
    if [ -f "Minimization/wwd.min.energy" ]; then # REMOVE FILE
       /bin/rm Minimization/wwd.min.energy
    fi
    if [ -f "Minimization/wwd.min.pdb" ]; then # REMOVE FILE
       /bin/rm Minimization/wwd.min.pdb
    fi
    ../../applications/protomol-app/protomol Minimization/wwd.min.conf >& Minimization/output.regression
    echo "Checking setup..."
    x=`diff -I 2007 -I Timing Minimization/output.regression Minimization/output.correct`
    if [ ! -f "Minimization/output.regression" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff Minimization/wwd.energy Minimization/wwd.min.energy`
    if [ ! -f "Minimization/wwd.min.energy" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN wwd.min.energy, COMPARE TO wwd.energy"
       failbit=1
    fi
    echo "Checking final positions..."
    x=`diff -I REMARKS Minimization/wwd.pdb Minimization/wwd.min.pdb`
    if [ ! -f "Minimization/wwd.min.pdb" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION POSITIONS ARE IN wwd.min.pdb, COMPARE TO wwd.pdb"
       failbit=1
    fi
}


# HESSIAN
hessian() {
    echo "******** RUNNING HESSIAN TESTS ********"
    if [ -f "Hessian/output.regression" ]; then # REMOVE FILE
       /bin/rm Hessian/output.regression
    fi
    if [ -f "Hessian/wwd.hess.energy" ]; then # REMOVE FILE
       /bin/rm Hessian/wwd.hess.energy
    fi
    if [ -f "Hessian/wwdaevect.bin.dat" ]; then # REMOVE FILE
       /bin/rm Hessian/wwdaevect.bin.dat
    fi
    ../../applications/protomol-app/protomol Hessian/wwd.hess.conf >& Hessian/output.regression
    echo "Checking setup..."
    x=`diff -I 2007 -I Timing Hessian/output.regression Hessian/output.correct`
    if [ ! -f "Hessian/output.regression" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff Hessian/wwd.energy Hessian/wwd.hess.energy`
    if [ ! -f "Hessian/wwd.hess.energy" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN wwd.hess.energy, COMPARE TO wwd.energy"
       failbit=1
    fi
    echo "Checking final eigenvectors..."
    x=`diff Hessian/wwdaevect.bin.dat NormalMode/wwdaevect.bin.dat`
    if [ ! -f "Hessian/wwdaevect.bin.dat" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  EIGENVECTORS ARE DIFFERENT."
       failbit=1
    fi
}

# NORMAL MODE RUN
normalmode() {
    echo "******** RUNNING NORMALMODE ********"
    if [ -f "NormalMode/output.regression" ]; then # REMOVE FILE
       /bin/rm NormalMode/output.regression
    fi
    if [ -f "NormalMode/wwd.anm.energy" ]; then # REMOVE FILE
       /bin/rm NormalMode/wwd.anm.energy
    fi
    if [ -f "NormalMode/wwd.anm.pdb" ]; then # REMOVE FILE
       /bin/rm NormalMode/wwd.anm.pdb
    fi
    ../../applications/protomol-app/protomol NormalMode/wwd.anm.conf >& NormalMode/output.regression
    echo "Checking setup..."
    x=`diff -I 2007 -I Timing NormalMode/output.regression NormalMode/output.correct`
    if [ ! -f "NormalMode/output.regression" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff NormalMode/wwd.energy NormalMode/wwd.anm.energy`
    if [ ! -f "NormalMode/wwd.anm.energy" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN wwd.anm.energy, COMPARE TO wwd.energy"
       failbit=1
    fi
    echo "Checking final positions..."
    x=`diff -I REMARKS NormalMode/wwd.pdb NormalMode/wwd.anm.pdb`
    if [ ! -f "NormalMode/wwd.anm.pdb" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION POSITIONS ARE IN wwd.anm.pdb, COMPARE TO wwd.pdb"
       failbit=1
    fi
}


echo "********* RUNNING NORMAL MODE REGRESSION TESTS *********"
if [ $1 == "all" ]; then
   echo "********* RUNNING ALL TESTS *********"
   minimization
   hessian
   normalmode
elif [ $1 == "cgminimize" ]; then
   minimization
elif [ $1 == "hessian" ]; then
   hessian
elif [ $1 == "normalmode" ]; then
   normalmode
fi

if [ $failbit == 1 ]; then
   echo "ONE OR MORE REGRESSION TESTS FAILED.  SEE ABOVE."
else
   echo "ALL REGRESSION TESTS PASSED."
fi
