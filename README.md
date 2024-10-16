This repository contains the PALET algorithm which is used to solve systems of
linear equations exactly.

Files and folders in this distribution:

    README.md               this file

    PALET.cpp               C++ file containing the PALET algorithm

    makefile                The makefile which compiles PALET.cpp

Dependencies:

    GMP, MPFR, givaro, OpenBLAS, fflas-ffpack

Usage:

    To compile PALET, please update the makefile with the correct paths to its 
    dependencies. Then type “make” in this folder. This will create the PALET 
    executable.

    To run, use ./PALET [Matrixfile.txt] where [Matrixfile.txt] is the path to 
    the coefficient matrix of the system you would like to solve. The RHS is the 
    unit vector e1.

Authors:

    Kelsey Kitzmiller, Chris Lourenco, Erick Moreno-Centeno

Contact:

    Please contact Kelsey Kitzmiller (k.kitzmiller@tamu.edu)
    
Copyright: 

    This software is copyright by Kelsey Kitzmiller, Chris Lourenco, 
    Erick Moreno-Centeno. All Rights Reserved.

Disclaimer:

    This code is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
