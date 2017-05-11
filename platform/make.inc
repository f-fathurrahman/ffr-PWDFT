MAKE = make
AR = ar
F90 = ifort

# Available macro:
#   -D_VERBOSE_ARNOLDI
#   -D_TIME_APPLY_HAM
F90_OPTS = -fpp -D_VERBOSE_ARNOLDI # -check bounds -check uninit

#LIBS =  -Wl,--start-group $(MKLROOT)/lib/32/libmkl_intel.a $(MKLROOT)/lib/32/libmkl_sequential.a $(MKLROOT)/lib/32/libmkl_core.a -Wl,--end-group -lpthread
LIB_MKL = -L/opt/intel/Compiler/11.1/038/mkl/lib/32/ -lmkl_intel -lmkl_core -lmkl_sequential
LIB_FFTW3 = -lfftw3
LIBS = $(LIB_MKL) $(LIB_FFTW3)

CUDA_LIBS = -L/home/marisa1709efefer/mysoftwares/cuda/lib -lcudartemu -lcufftemu -lcublasemu

