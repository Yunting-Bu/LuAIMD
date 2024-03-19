#!/bin/bash

cd scr

gfortran -c machina_basic.f90
gfortran -c read_input.f90
gfortran -c basis.f90
gfortran -c cint.f90
gfortran -c output.f90
gfortran -c scf.f90
gfortran -c gradient.f90
gfortran -c md.f90
gfortran -c main.f90

gfortran -fopenmp -g -o LuAIMD \
machina_basic.o read_input.o basis.o cint.o scf.o gradient.o md.o output.o main.o \
lapack.a blas.a \
/home/byt/software/libcint/lib/libcint.so  # Replace PATH with the actual path.

rm *.o
rm *.mod

mv LuAIMD ..

cd ..

echo "Completed!"
