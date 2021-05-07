
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
modmain.f90 \
oblas_stub.f90 \
zfftifc.f90 \
libxcifc_stub.f90 \
modxcifc.f90 \
read_input.f90 \
read_species_files.f90 \
add_lorb_cnd.f90 \
default_atoms.f90 default_apwlo.f90 default_muffin_tins.f90 \
default_lattice.f90 default_plotting.f90 default_mixing.f90 \
default_misc.f90 default_symmetry.f90 default_force_stress.f90 \
default_hamiltonian.f90 default_states.f90 default_spin.f90 \
default_qpoints.f90 default_kpoints.f90 default_gkvectors.f90 \
default_convergence.f90 default_oep_hf.f90 \
default_charge_moment_current.f90 \
default_electric_vector_pot.f90 \
info_apwlo.f90 info_lattice.f90 info_muffin_tins.f90 \
init_am_variables.f90 init_chgden_pot_xc.f90  init_Gk_vectors.f90 \
init_radial_meshes.f90 init_APW_LO.f90  init_crystal_structure.f90  \
init_gvector_arrays.f90 init_spin_variables.f90 init_atoms_cores.f90 \
init_eigensystems.f90 init_idx_atom_species.f90 init_timing.f90 \
init_charges_states.f90 init_forces.f90 init_kpoints.f90 \
init_vector_field_E_A.f90 \
r3mv.f90 gauntyry.f90 genppts.f90 \
potnucl.f90 gencfun.f90 \
gensfacgp.f90 genffacgp.f90 genjlgprmt.f90 \
genylmg.f90 genylmv.f90 sbessel.f90 radnucl.f90 \
vecfbz.f90 r3frac.f90 gaunt.f90 checkmt.f90 \
factr.f90 factnm.f90 wigner3j.f90 mtdmin.f90 \
r3minv.f90 allatoms.f90 gengclg.f90 reciplat.f90 \
genidxlo.f90 genrmesh.f90 gengkvec.f90 findngkmax.f90 \
xc_am05.f90  xc_c_tb09.f90  xc_pbe.f90  xc_pwca.f90 \
xc_pzca.f90  xc_vbh.f90  xc_wc06.f90  xc_xalpha.f90 \
c_pbe.f90  c_pbe_gcor.f90  x_pbe.f90  x_wc06.f90 \
rfint.f90 rfsht.f90 rbsht.f90 gradrf.f90 wsplintp.f90 wsplint.f90 \
r3cross.f90 atom.f90 gengvec.f90 sortidx.f90 \
gridsize.f90 splintwp.f90 fderiv.f90 rdirac.f90 rdiracint.f90 splint.f90 \
polynm.f90 gradrfmt.f90 rfmtsm.f90 ztorfmt.f90 \
rtozfmt.f90 spline.f90 nfftifc.f90 rfmtlm.f90 \
ztorflm.f90 fsmooth.f90 gradzfmt.f90 rtozflm.f90 clebgor.f90 \
symmetry.f90 findsymlat.f90 findsymcrys.f90 findsymsite.f90 findsym.f90 \
r3mtm.f90 i3mdet.f90 r3mm.f90 r3mtv.f90 i3minv.f90

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
	$(F90) $(F90_OPTS) main_lapwdft.f90 -o $(EXE) $(OBJ) $(LIB_W90) $(LIB_LPK) $(LIB_FFT) \
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
