#!/bin/bash

MOD='mods'

gfortran $MOD/modMT.f90 $MOD/modstatPAni.f90 $MOD/modNMRPAni.f90 $MOD/modIOUtil.f90 numNMRPAni.f90 \
	-o nump && ./nump



