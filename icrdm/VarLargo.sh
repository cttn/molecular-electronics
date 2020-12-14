#!/bin/bash

#Comienza el loop en el numero de anillos.
for i in $(seq 40 40 400);
  do
  mkdir $i
  cp *.f90 $i
  cd $i

gfortran decbipol-ng19.f90 -o decbipol
./decbipol

     echo "module diag_icrdm
  integer, parameter :: numanillos=$i
  integer, parameter :: num_conf=1,inc=0
!El programa empieza a nombrar las configuraciones desde inc+1
end module diag_icrdm
program dummy
end program dummy" > diag_icrdm_mod.f90

gfortran diag_icrdm_mod.f90

gfortran diagbipol-ng19.f90 -o diagb -llapack
./diagb

gfortran icrdm-ng19.f90 -o icrdm
./icrdm

  pwd
  cd ..
  done
