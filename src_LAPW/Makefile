
include make.inc.gfortran

SRC_main = main_lapwdft.f90

SRC = m_lattice.f90 \
m_atoms.f90 \
m_muffin_tins.f90 \
m_spin.f90 \
m_electric_vector_pot.f90 \
m_symmetry.f90 \
m_gvectors.f90 \
m_kpoints.f90 \
m_gkvectors.f90 \
m_qpoints.f90 \
m_sht.f90 \
m_density_pot_xc.f90 \
m_mixing.f90 \
m_charge_moment_current.f90 \
m_apwlo.f90 \
m_hamiltonian.f90 \
m_states.f90 \
m_core_states.f90 \
m_energy.f90 \
m_force_stress.f90 \
m_convergence.f90 \
m_dos_optics_response.f90 \
m_plotting.f90 \
m_oep_hf.f90 \
m_bse.f90 \
m_timing.f90 \
m_constants.f90 \
m_misc.f90 \
read_input.f90 \
read_species_files.f90 \
add_lorb_cnd.f90

OBJ = $(SRC:.f90=.o)
EXE = lapwdft.x

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
lib:	$(OBJ)
	ar rcs liblapwdft.a $(OBJ)

main:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_W90) $(LIB_LPK) $(LIB_FFT) \
    $(LIB_libxc) $(LIB_SYS)

clean:
	rm -f *.o *.mod 

mpi_stub.o: mpi_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c mpi_stub.f90

libxcifc_stub.o: libxcifc_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c libxcifc_stub.f90

mkl_stub.o: mkl_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c mkl_stub.f90

oblas_stub.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c oblas_stub.f90

blis_stub.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c blis_stub.f90

zfftifc.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c zfftifc.f90
