MAKE = make
AR = ar

F77 = pgf77
F77_OPTS = -Minform=inform -fast

F90 = pgf90
F90_OPTS = -Minform=inform -fast -Mpreprocess

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3

LIB_XC     = -L/home/efefer/mysoftwares/libxc-3.0.0/lib -lxcf90 -lxc
LIBS       = $(LIB_LINALG) $(LIB_FFTW3) $(LIB_XC)
