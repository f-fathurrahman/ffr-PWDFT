#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo
  echo "ERROR"
  echo "Need two parameters: main file and compiler name (and or options)"
  echo "Example: ./build.sh myfile.f90 \"pgf90 -O2\""
  echo
  exit 1
fi

INC="-I../src/"
LIB="../src/libmain.a -L/home/efefer/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lfftw3"
LIBXC="-L$HOME/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc"

bas=`basename $1 .f90`

# remove the previous executable
rm -vf $bas.x

$2 $INC $1 $LIB $LIBXC -o $bas.x

#gfortran -Wall -O3 -ffree-form $INC $1 $LIB -o $bas.x
echo "Test executable: $bas.x"

