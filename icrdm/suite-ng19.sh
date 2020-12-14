#! /bin/bash
# Bash script

gfortran diag_icrdm_mod.f90

gfortran decbipol-ng19.f90 -o decbipol
./decbipol


gfortran diagbipol-ng19.f90 -o diagb -llapack
./diagb

#pwd

gfortran icrdm-ng19.f90 -o icrdm
./icrdm


echo '  '
echo 'Terminado. :-)'
