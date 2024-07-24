# PALET
PALET is a software package used to solve a system of linear equations exactly using the PALET algorithm.

Files and folders in this distribution:
	README.md	     This file
	PALET.cpp	     C++ file containing the PALET algorithm
	makefile	     The makefile which compiles PALET.cpp
	Matrixfiles	   Folder of random/magic matrices used to test PALET

Dependencies:
	GMP
	MPFR
	GIVARO
	OPENBLAS
	FFLAS-FFPACK

To compile PALET, please update the makefile with the correct paths to its dependencies. Then type “make” in this folder. This will create the PALET executable.

To run, use ./PALET [Matrixfile.txt]

Example usage: ./PALET Matrixfiles/Magic501.txt

This will solve the system Ax=b where A is the matrix stored in Magic501.txt and b is the e1 unit vector.
