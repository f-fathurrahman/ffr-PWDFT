
**This is a work in progress**

A poor man's program to carry out electronic structure calculations
based on plane wave basis set and pseudopotential.

This program is written using simple Fortran 90. I have tried to make
the code as simple and transparent as possible.

To compile the program you need to edit the top part of the Makefile
to include compiler or platform specific details.

Currently it has been tested using the following compilers:
- [`gfortran`](https://gcc.gnu.org/fortran/)
- [`g95`](http://www.g95.org)
- [`ifort`](https://software.intel.com/en-us/fortran-compilers)
- [`pgf90`](https://www.pgroup.com/products/community.htm)
- [`sunf95`](http://www.oracle.com/technetwork/server-storage/developerstudio/downloads/index.html)
- [`flang`](https://github.com/flang-compiler/flang)

As of now, the compilation will not produce any executables yet. It
will only produce a static library file named `libmain.a`.
Several examples of how to use the library are given in the directory
`tests`.

See also: [PWDFT.jl](https://github.com/f-fathurrahman/PWDFT.jl)
