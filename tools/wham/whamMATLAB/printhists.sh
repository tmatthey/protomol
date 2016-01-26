#!/bin/bash

cd ..
tempdirs=`find *K -maxdepth 0`

for tempdir in $tempdirs
  do
  echo -e $tempdir
  cd $tempdir
  gawk '{print $3}' dihedral.out > $tempdir.dihedrals.out
  xmgrace $tempdir.dihedrals.out &  
  cd ..
done