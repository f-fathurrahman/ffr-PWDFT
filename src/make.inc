MAKE = make
AR   = ar

F90      = g95
F90_OPTS = -Wall -ftrace=full -cpp -O3 -fno-second-underscore

F77 = g95
F77_OPTS = -Wall -ftrace=full -cpp -O3 -fno-second-underscore

LIB_LINALG = -lblas -llapack
LIB_FFTW3  = -lfftw3
LIB_XC     = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS       = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)

