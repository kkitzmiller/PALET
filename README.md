# PALET
PALET is a software package used to solve a system of linear equations exactly using the PALET algorithm.

Files in this package:

[README.md](/README.md)	     This file

[PALET.cpp](/PALET.cpp)	     C++ file containing the PALET algorithm

[makefile](/makefile)	     The makefile which compiles PALET.cpp

Dependencies: [GMP](https://github.com/yuhangwang/GMP), [MPFR](https://github.com/alisw/MPFR), [givaro](https://github.com/linbox-team/givaro), [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS), [fflas-ffpack](https://github.com/linbox-team/fflas-ffpack)

To compile PALET, please update the makefile with the correct paths to its dependencies. Then type “make” in this folder. This will create the PALET executable.

To run, use ./PALET [Matrixfile.txt] where [Matrixfile.txt] is the path to the coefficient matrix of the system you would like to solve. The RHS is the unit vector e1.
