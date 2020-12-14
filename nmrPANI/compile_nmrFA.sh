#!/bin/bash

gfortran MOD_MT.f90 MOD_nmrFA.f08 -llapack -o nmrFA $1
