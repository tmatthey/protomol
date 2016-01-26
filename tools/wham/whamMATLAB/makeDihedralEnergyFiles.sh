#!/bin/bash

cd ..
tempdirs=`find . -name '*k' -type d -maxdepth 1`

for tempdir in $tempdirs
  do
  echo -e $tempdir
  cd $tempdir
  rm -rf $tempdir.potential.out $tempdir.dihedral.out $tempdir.dihedralE.out 
  gawk '{print $2}' energy.out > $tempdir.potential.out
  gawk '{print $3}' dihedral.out > $tempdir.dihedral.out
  gawk '{print $4}' dihedral.out > $tempdir.dihedralE.out  
  cd ..
done