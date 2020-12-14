#!/bin/bash

## Compila y ejecuta numNMRPAni.f90 para diferentes valores de eta.

MOD='mods'
numFile='numNMRPAni.f90'
tmpFile='r.f90'
tmpRes='etaRange.dat'

numLine='28'
#etaRange="4.0d0 4.5d0 12.0d0 15.0d0 20.0d0 23.0d0 24.0d0 25.0d0"
etaRange="1.0d0 3.0d0 6.0d0 8.0d0 12.0d0 16.0d0 20.0d0"

### Título
echo " Corriendo $numFile para un conjunto de valores de eta:"

### Verificando que en la linea indicada se defina eta
match=$(awk "NR==$numLine" $numFile | grep -o "eta =")
if [ ! "$match" == "eta =" ]
then
	echo "    "
	echo "Hubo un error: El numero de línea indicado: $numLine, no es el apropiado"
	exit 1
fi

### Core
for eta in $etaRange 
do
	### Modifica la linea del archivo que define el valor de eta
	line="$numLine"'s/= .*/='$eta'/'
	sed "$line" numNMRPAni.f90 > $tmpFile

	### Compila y ejecuta el codigo de fortran modificado
    gfortran $MOD/modMT.f90 $MOD/modstatPAni.f90 $MOD/modNMRPAni.f90 $MOD/modIOUtil.f90 $tmpFile -o xR && ./xR

	### Cambia el nombre a los resultados para ordenarlos
	if [ -f "$tmpRes" ]
	then
	  mv $tmpRes 'numPAni_'$eta'.dat'
    else
	  echo "Hubo algun error, $tmpRes desaparecido..."
	  exit 1
	fi

	echo "      ---> $eta Listo"
done


