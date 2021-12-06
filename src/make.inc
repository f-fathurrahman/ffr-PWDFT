MAKE = make
AR = ar

F77 = nvfortran
F77_OPTS = -Wall -O3 -cpp

F90 = nvfortran
F90_OPTS = -Wall -O3 -cpp

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIB_XC     = -L/home/efefer/mysoftwares/libxc-5.1.5/lib -lxcf90 -lxc
LIBS       = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)
