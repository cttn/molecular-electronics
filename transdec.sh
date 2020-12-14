#!/bin/bash

# Algunas variables 
programa='transdec-0.8.3.f90'
ejecut='tran'
log='tran.log'
num_lines=1

# Indica que programa estoy corriendo
echo '   ' >> $log
date >> $log

echo '   '
echo "  ###################################"
echo "      Corriendo $programa"
echo "  ###################################"
echo '   '

# Compila y ejecuta Fortran. Registra la salida a terminal
gfortran $programa -o $ejecut
script -a $log -c ./$ejecut -f -q

# Notifica lo último que aparece en la terminal
lastline="$(tail -n $num_lines $log)"
notify-send $programa -i /home/carlos/Imágenes/Icons/gears-trans.ico  "El programa finalizó su ejecución:  $lastline"


