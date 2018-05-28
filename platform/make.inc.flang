MAKE = make
AR = ar

F77 = flang -Mpreprocess
F77_OPTS = -Wall -O3

F90 = flang -Mpreprocess
F90_OPTS = -Wall -O3

LIB_LINALG = -L/home/efefer/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LIB_FFTW3  = -lfftw3
LIB_XC     = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS       = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)
