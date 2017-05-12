#include platform/make.inc.ifort
#include platform/make.inc.gfortran
include platform/make.inc.g95
#include platform/make.inc.pgi

SRC = \
m_constants.f90 \
m_options.f90 \
m_cell.f90 \
m_atoms.f90 \
m_realspace.f90 \
m_PWGrid.f90 \
init_atoms_xyz.f90 \
inv_m3x3.f90 \
det_m3x3.f90 \
init_PWGrid.f90 \
info_PWGrid.f90 \
mm_to_nn.f90 \
init_gvec.f90 \
init_gvecw.f90 \
dealloc_PWGrid.f90 \
fft_fftw3.f90 \
init_rgrid.f90 \
dealloc_realspace.f90 \
xsf.f90 \
Poisson_solve_fft.f90


OBJ = $(SRC:.f90=.o) $(SRC:.f=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F90) -O3 -c $<


# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

clean:
	rm -rf *.o *.mod libmain.a *.x


