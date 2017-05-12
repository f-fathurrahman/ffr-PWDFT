#include platform/make.inc.ifort
#include platform/make.inc.gfortran
#include platform/make.inc.g95
include platform/make.inc.pgi

SRC = \
m_constants.f90 \
m_options.f90 \
m_cell.f90 \
m_atoms.f90 \
m_realspace.f90 \
m_PWGrid.f90 \
m_states.f90 \
m_Ps_HGH.f90 \
m_PsPot.f90 \
fft_param.f90 \
fft_support.f90 \
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
Poisson_solve_fft.f90 \
init_strfact.f90 \
init_states.f90 \
init_PsPot.f90 \
op_K.f90 \
m_hamiltonian.f90 \
alloc_hamiltonian.f90 \
dealloc_hamiltonian.f90 \
op_V_loc.f90 \
op_H.f90 \
calc_grad.f90 \
m_energies.f90 \
info_energies.f90 \
update_potentials.f90 \
LDA_VWN.f90 \
prec_Gv2.f90 \
ortho_gram_schmidt.f90 \
ortho_check.f90 \
calc_Rhoe_R.f90 \
random_wfc.f90 \
calc_energies.f90 \
init_V_ps_loc_harmonic.f90 \
info_energies.f90 \
KS_solve_Emin_pcg.f90 \
calc_Ewald.f90 \
init_V_ps_loc_G.f90 \
info_atoms.f90 \
info_PsPot.f90 \
dealloc_atoms.f90 \
dealloc_PsPot.f90


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


