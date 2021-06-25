clear 
close all
clc

mex -O -largeArrayDims -output computeGREEN_K_f90_mexed integration_fortran_K.f90 
