#!/bin/bash

hist_out='histotef-ng19_tef.dat' #salida de cada histotef
results='VarG-ng19.dat'      #acumulacion de resultados
inicial=0
saltos=1
final=1000


# Borra archivos previos para no mezclar datos
if [ -f $results ]
  then
   rm $results
fi

#Comienza el loop en el indice (alpha) de los puntos de Gamma.
for i in $(seq $inicial $saltos $final);
  do
     echo ' '
     echo "$i de $final ----------------------------"
     echo "module histo_mod_G
  integer, parameter :: i_alpha=$i
end module histo_mod_G
program dummy_hist
end program dummy_hist" > histo_mod_G.f90

  gfortran histo_mod_G.f90
  gfortran histotef-ng19G.f90 -o hist
  ./hist

  cat $hist_out >> $results

done 
