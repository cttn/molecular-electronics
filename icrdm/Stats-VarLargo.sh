#!/bin/bash

hist_out='histotef-ng19_tef.dat' #salida de cada histotef
results='VarLargo-ng19.dat'      #acumulacion de resultados
inicial=40
saltos=40
final=400


#Comienza el loop en el numero de anillos.
for i in $(seq $inicial $saltos $final);
  do
  cd $i
  cp ../histotef-ng19.f90 histotef-ng19.f90

#  gfortran diag_icrdm_mod.f90
  gfortran histotef-ng19.f90 -o hist
  ./hist

  cat $hist_out >> ../$results
  cd ..
done 
