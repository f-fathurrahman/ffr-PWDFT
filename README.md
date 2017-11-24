
**This is a work in progress**

A poor man's program to carry out electronic structure calculations
based on plane wave basis set and pseudopotential.

This program is written using simple Fortran 90. I have tried to make
the code as simple and transparent as possible.

To compile the program you need to edit the top part of the Makefile
to include compiler or platform specific details. This program is currently
tested using 5 compilers:

- GNU Fortran compiler
- G95 compiler
- Intel Fortran compiler
- PGI Fortran compiler
- Sun Fortran compiler

As of now, the compilation will not produce any executables yet. It
will only produce a static library file named `libmain.a`.
Several examples of how to use the library are given in the directory
`tests`.

