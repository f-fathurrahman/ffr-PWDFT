MAKE = make
AR   = ar

F90      = g95
F90_OPTS = -Wall -ftrace=full -cpp -O3

LIB_LINALG = -lblas -llapack
LIB_FFTW3  = -lfftw3
LIBS       = $(LIB_LINALG) $(LIB_FFTW3)
