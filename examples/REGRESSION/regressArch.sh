#! /bin/sh

failbit=0
# MINIMIZATION
minimization() {
    echo "******** RUNNING MINIMIZATION ********"
    if [ -f "Minimization/arch1/output.regression" ]; then # REMOVE FILE
       /bin/rm Minimization/arch1/output.regression
    fi
    if [ -f "Minimization/arch2/output.regression" ]; then # REMOVE FILE
       /bin/rm Minimization/arch2/output.regression
    fi
    if [ -f "Minimization/arch1/wwd.min.energy" ]; then # REMOVE FILE
       /bin/rm Minimization/arch1/wwd.min.energy
    fi
    if [ -f "Minimization/arch2/wwd.min.energy" ]; then # REMOVE FILE
       /bin/rm Minimization/arch2/wwd.min.energy
    fi
    if [ -f "Minimization/arch1/wwd.min.pdb" ]; then # REMOVE FILE
       /bin/rm Minimization/arch1/wwd.min.pdb
    fi
    if [ -f "Minimization/arch2/wwd.min.pdb" ]; then # REMOVE FILE
       /bin/rm Minimization/arch2/wwd.min.pdb
    fi
    $1 Minimization/wwd.min.arch1.conf >& Minimization/arch1/output.regression
    $2 Minimization/wwd.min.arch2.conf >& Minimization/arch2/output.regression
    echo "Checking setup..."
    x=`diff -I Timing Minimization/arch1/output.regression Minimization/arch2/output.regression`
    if [ ! -f "Minimization/arch1/output.regression" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Minimization/arch2/output.regression" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN.  CHECK EXECUTABLE PATH."
       failbit=1
    elif [ -n ${x:0:0} ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff Minimization/arch1/wwd.min.energy Minimization/arch2/wwd.min.energy`
    if [ ! -f "Minimization/arch1/wwd.min.energy" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Minimization/arch2/wwd.min.energy" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN.  CHECK EXECUTABLE PATH."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN Minimization/arch1/wwd.min.energy, COMPARE TO Minimization/arch2/wwd.min.energy"
       failbit=1
    fi
    echo "Checking final positions..."
    x=`diff Minimization/arch1/wwd.min.pdb Minimization/arch2/wwd.min.pdb`
    if [ ! -f "Minimization/arch1/wwd.min.pdb" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Minimization/arch2/wwd.min.pdb" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION POSITIONS ARE IN Minimization/arch1/wwd.min.pdb, COMPARE TO Minimization/arch2/wwd.min.pdb"
       failbit=1
    fi
}


# HESSIAN
hessian() {
    echo "******** RUNNING HESSIAN TESTS ********"
    if [ -f "Hessian/arch1/output.regression" ]; then # REMOVE FILE
       /bin/rm Hessian/arch1/output.regression
    fi
    if [ -f "Hessian/arch2/output.regression" ]; then # REMOVE FILE
       /bin/rm Hessian/arch2/output.regression
    fi
    if [ -f "Hessian/arch1/wwd.hess.energy" ]; then # REMOVE FILE
       /bin/rm Hessian/arch1/wwd.hess.energy
    fi
    if [ -f "Hessian/arch2/wwd.hess.energy" ]; then # REMOVE FILE
       /bin/rm Hessian/arch2/wwd.hess.energy
    fi
    if [ -f "Hessian/arch1/wwdaevect.bin.dat" ]; then # REMOVE FILE
       /bin/rm Hessian/arch1/wwdaevect.bin.dat
    fi
    if [ -f "Hessian/arch2/wwdaevect.bin.dat" ]; then # REMOVE FILE
       /bin/rm Hessian/arch2/wwdaevect.bin.dat
    fi
    $1 Hessian/wwd.hess.arch1.conf >& Hessian/arch1/output.regression
    $2 Hessian/wwd.hess.arch2.conf >& Hessian/arch2/output.regression
    echo "Checking setup..."
    x=`diff -I Timing Hessian/arch1/output.regression Hessian/arch2/output.regression`
    if [ ! -f "Hessian/arch1/output.regression" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Hessian/arch2/output.regression" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n ${x:0:0} ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff Hessian/arch1/wwd.hess.energy Hessian/arch2/wwd.hess.energy`
    if [ ! -f "Hessian/arch1/wwd.hess.energy" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Hessian/arch2/wwd.hess.energy" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN Hessian/arch1/wwd.hess.energy, COMPARE TO Hessian/arch2/wwd.hess.energy"
       failbit=1
    fi
    echo "Checking final eigenvectors..."
    x=`diff Hessian/arch1/wwdaevect.bin.dat Hessian/arch2/wwdaevect.bin.dat`
    if [ ! -f "Hessian/arch1/wwdaevect.bin.dat" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "Hessian/arch2/wwdaevect.bin.dat" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION EIGENVECTORS ARE DIFFERENT."
       failbit=1
    fi
}

# NORMAL MODE RUN
normalmode() {
    echo "******** RUNNING NORMALMODE ********"
    if [ -f "NormalMode/arch1/output.regression" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch1/output.regression
    fi
    if [ -f "NormalMode/arch2/output.regression" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch2/output.regression
    fi
    if [ -f "NormalMode/arch1/wwd.anm.energy" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch1/wwd.anm.energy
    fi
    if [ -f "NormalMode/arch2/wwd.anm.energy" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch2/wwd.anm.energy
    fi
    if [ -f "NormalMode/arch1/wwd.anm.pdb" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch1/wwd.anm.pdb
    fi
    if [ -f "NormalMode/arch2/wwd.anm.pdb" ]; then # REMOVE FILE
       /bin/rm NormalMode/arch2/wwd.anm.pdb
    fi
    $1 NormalMode/wwd.anm.arch1.conf >& NormalMode/arch1/output.regression
    $2 NormalMode/wwd.anm.arch2.conf >& NormalMode/arch2/output.regression
    echo "Checking setup..."
    x=`diff -I Timing NormalMode/arch1/output.regression NormalMode/arch2/output.regression`
    if [ ! -f "NormalMode/arch1/output.regression" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "NormalMode/arch2/output.regression" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n ${x:0:0} ]; then 
       echo "PASSED"
    else
       echo "FAILED.  SETUP IS WRONG."
       failbit=1
    fi
    echo "Checking energies..."
    x=`diff NormalMode/arch1/wwd.anm.energy NormalMode/arch2/wwd.anm.energy`
    if [ ! -f "NormalMode/arch1/wwd.anm.energy" ]; then 
       echo "FAILED.  YOUR SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "NormalMode/arch2/wwd.anm.energy" ]; then 
       echo "FAILED.  BENCHMARK SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x ]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION ENERGIES ARE IN NormalMode/arch1/wwd.anm.energy, COMPARE TO NormalMode/arch2/wwd.energy"
       failbit=1
    fi
    echo "Checking final positions..."
    x=`diff -I REMARKS NormalMode/arch1/wwd.anm.pdb NormalMode/arch2/wwd.anm.pdb`
    if [ ! -f "NormalMode/arch1/wwd.anm.pdb" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ ! -f "NormalMode/arch2/wwd.anm.pdb" ]; then 
       echo "FAILED.  SIMULATION DID NOT RUN."
       failbit=1
    elif [ -n $x]; then
       echo "PASSED"
    else
       echo "FAILED.  REGRESSION POSITIONS ARE IN NormalMode/arch1/wwd.anm.pdb, COMPARE TO NormalMode/arch2/wwd.anm.pdb"
       failbit=1
    fi
}


echo "********* RUNNING NORMAL MODE REGRESSION TESTS *********"
if [ $1 == "all" ]; then
   echo "********* RUNNING ALL TESTS *********"
   minimization $2 $3
   hessian $2 $3
   normalmode $2 $3
elif [ $1 == "cgminimize" ]; then
   minimization $2 $3
elif [ $1 == "hessian" ]; then
   hessian $2 $3
elif [ $1 == "normalmode" ]; then
   normalmode $2 $3
fi

if [ $failbit == 1 ]; then
   echo "ONE OR MORE REGRESSION TESTS FAILED.  SEE ABOVE."
else
   echo "ALL REGRESSION TESTS PASSED."
fi
