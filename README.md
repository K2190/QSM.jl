# QSM MATLAB to Julia using 

This folder contains the Julia code converted from MATLAB code taken from (https://github.com/kamesy/QSM.m)

## Folder Structure  
julia<br>
-- examples<br>
-- src<br>
---- bgremove<br>
---- inversion<br>
------ masks_weights<br>
------ mex<br>
---- io<br>
------ dicomSorter<br>
---- unwrap<br>
------ mex<br>
---- utils<br>
------ error_metrics<br>
------ fd<br>
------ kernels<br>
------ poisson_solver<br>
-- third_party<br>
---- lsmr<br>
---- lsqrSOL<br>
---- NifTI<br>
sharedLibLinux : These files contains the C files used in Julia ccall <br>
-- *.c<br>
-- Makefile<br>
install_pkgs.jl

## Julia Packages 
Use install_pkgs.jl to install packages used in the code


## Developing Shared Library

MATLAB builds individual files and it uses the MEX approach to link the MATLAB .m files with the .C files. In Julia, all the C files are collected together in a single folder and a single shared library is generated. 

#### Netlib BLAS/LAPACK
C files uses APIs from BLAS and LAPACK. Download and build these required libs using the instructions given in <br>
http://www.netlib.org/blas/<br>
http://www.netlib.org/lapack/<br>



#### Build and Combining all libraries

Modify the Makefile based on the location of above libs.

#### make
This will generate all teh required objects and .so files

Following command used to combine all the external libs and *.o files libjulia_qsm.so. This so is referred in Julia ccall.

#### gcc -shared *.o libblas.a liblapack.a -lm -lgomp -o libjulia_qsm.so

Add the path of this generated libjulia_qsm.so file in ../julia/examples/run_qsm.jl (Line 10)

## Running the examples

#### Steps to run the examples

1. Start a terminal and setup the required environment for MATLAB
2. Start visual studio code. (Assuming that required Julia extension are already installed)
3. In the julia file in .../julia/examples/run_qsm.jl, setup the path to the models given in .../julia/models
4. Two set of models are given. One for 50x50x50 and other for 160x160x160
5. Set the output file name to be saved
6. Run .../julia/examples/run_qsm.jl
