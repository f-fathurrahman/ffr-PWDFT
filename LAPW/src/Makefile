
#include make.inc.gfortran
include make.inc.gfortran.dynlib
#include make.inc.flang

SRC_main = main_lapwdft.f90

SRC = m_lattice.f90 \
m_atoms.f90 \
m_atomic_species.f90 \
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
modrandom.f90 \
oblas_stub.f90 \
zfftifc_fftw.f90 \
libxcifc.f90 \
modxcifc.f90 \
read_input.f90 \
read_species_files.f90 \
add_lorb_cnd.f90 \
default_atoms.f90 default_apwlo.f90 default_muffin_tins.f90 \
default_lattice.f90 default_plotting.f90 default_mixing.f90 \
default_misc.f90 default_symmetry.f90 default_force_stress.f90 \
default_hamiltonian.f90 default_states.f90 default_spin.f90 \
default_qpoints.f90 default_kpoints.f90 default_gvectors.f90 default_gkvectors.f90 \
default_convergence.f90 default_oep_hf.f90 \
default_charge_moment_current.f90 \
default_electric_vector_pot.f90 \
default_core_states.f90 \
default_density_pot_xc.f90 \
default_dos_optics_response.f90 \
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
my_findsymlat.f90 my_findsymcrys.f90 my_findsym.f90 my_findsymsite.f90 \
r3mtm.f90 i3mdet.f90 r3mm.f90 r3mtv.f90 i3minv.f90 axangrot.f90 \
init0.f90 nuclei.f90 energynn.f90 \
genshtmat.f90 zpotcoul.f90 sphcover.f90 genrlmv.f90 sctovec.f90 \
rminv.f90 zminv.f90 timesec.f90 \
mixerifc.f90 mixadapt.f90 mixlinear.f90 mixbroyden.f90 \
init1.f90 gndstate.f90 \
rhoinit.f90 maginit.f90 \
potks.f90 genvsig.f90 \
mixpack.f90 gencore.f90 \
linengy.f90 genapwlofr.f90 \
genapwfr.f90 genlofr.f90 gensocfr.f90 \
genevfsv.f90 occupy.f90 \
putoccsv.f90 rhomag.f90 \
energy.f90 writeeval.f90 writefermi.f90 writeforces.f90 writesym.f90 \
writelinen.f90 writeengy.f90 writechg.f90 writestate.f90 writemom.f90 \
potcoul.f90 rfirsm.f90 genbs.f90 genws.f90 potxc.f90 potxc_default.f90 \
oepmain.f90 gentau.f90 projsbf.f90 rfpack.f90 findband.f90 olprad.f90 hmlrad.f90 \
rschrodint.f90 \
eveqn.f90 \
putevalfv.f90 putevalsv.f90 \
getevecfv.f90 getevecsv.f90 \
stheta.f90 \
sdelta.f90 \
findswidth.f90 \
match.f90 \
rhomagk.f90 \
rhomagsh.f90 \
symrf.f90 \
rhocore.f90 \
charge.f90 \
rhonorm.f90 \
moment.f90 \
symrvf.f90 \
rfinp.f90 \
energykncr.f90 \
rfmtctof.f90 \
getkmat.f90 \
exxengy.f90 \
putevecfv.f90 \
putevecsv.f90 \
rfmtpack.f90 \
genzvclmt.f90 \
potefield.f90 \
bdipole.f90 \
rfmtftoc.f90 \
potxcmt.f90 \
potxcir.f90 \
oepvcl.f90 \
oepresk.f90 \
rfinpc.f90 \
rfint0.f90 \
initoep.f90 \
gentauk.f90 \
gentaucr.f90 \
eveqnit.f90 \
eveqnfv.f90 \
eveqnss.f90 \
eveqnsv.f90 \
rotzflm.f90 \
findkpt.f90 \
rotaxang.f90 \
axangsu2.f90 \
stheta_sq.f90 \
stheta_lr.f90 \
stheta_mp.f90 \
stheta_fd.f90 \
sbesseldm.f90 \
wavefmt.f90 \
zbsht.f90 \
symrfmt.f90 \
symrfir.f90 \
chargemt.f90 \
symrvfmt.f90 \
symrvfir.f90 \
rfinterp.f90 \
exxengyk.f90 \
wavefcr.f90 \
zfsht.f90 \
zpotclmt.f90 \
zfmtinp.f90 \
sdelta_sq.f90 \
sdelta_lr.f90 \
sdelta_mp.f90 \
sdelta_fd.f90 \
sphcrd.f90 \
curlrvf.f90 \
curden.f90 \
ggamt_1.f90 ggamt_sp_1.f90 ggamt_sp_2a.f90 ggamt_sp_2b.f90 ggamt_2a.f90 ggamt_2b.f90 \
ggair_1.f90 ggair_sp_2b.f90 ggair_sp_2a.f90 ggair_2a.f90 ggair_2b.f90 ggair_sp_1.f90 \
oepvclk.f90 \
getevalsv.f90 \
genwfsv.f90 \
genzrm.f90 \
rzfmtinp.f90 \
genzfrm.f90 \
rzfinp.f90 \
genzrho.f90 \
rfmtinp.f90 \
hmlaa.f90 hmlalo.f90 hmllolo.f90 hmlistl.f90 \
olpaa.f90 \
olpalo.f90 \
olplolo.f90 \
olpistl.f90 \
getevalfv.f90 \
hmlfv.f90 \
olpfv.f90 \
eveqnfvz.f90 \
eveqnfvr.f90 \
eveqnz.f90 lopzflm.f90 roteuler.f90 ylmrot.f90 erf.f90 hermite.f90 rotrfmt.f90 \
gengclgq.f90 zfinp.f90 curdenk.f90 grad2rfmt.f90 zmctmu.f90 ylmroty.f90 rotrflm.f90 rlmrot.f90 \
force.f90 forcek.f90 symveca.f90 \
info_apwlo.f90 info_crystal.f90 info_muffin_tins.f90 \
info_symmetry.f90 info_gvectors.f90 \
plot3d.f90 plotpt3d.f90 rfplot.f90 \
my_rdiracint.f90 my_rhoinit.f90 \
my_potks.f90 my_potks_no_symm.f90 my_potcoul.f90 my_genzvclmt.f90 my_zpotclmt.f90  my_zpotcoul.f90 \
my_potxc.f90 my_potxcir.f90  my_potxcmt.f90 \
my_rbsht.f90 my_rfip.f90 my_rfmtpack.f90 my_rfplot.f90 my_rfsht.f90 \
my_symrf.f90 my_symrfir.f90  my_symrfmt.f90 \
my_rotrfmt.f90 my_rotrflm.f90 \
z_to_rf_lm.f90 z_to_rf_mt.f90 rf_interp.f90  rf_mt_c_to_f.f90 \
r_to_zf_lm.f90 r_to_zf_mt.f90 \
my_potxcmt_lda_nospin.f90 \
my_potxcir_lda_nospin.f90 \
my_genshtmat.f90 \
print_tv0symc.f90


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
libdyn:	$(OBJ)
	gfortran -shared -o liblapwdft.so $(OBJ) $(LIBS)

lib:	$(OBJ)
	ar rcs liblapwdft.a $(OBJ)

main:	$(OBJ)
	$(F90) $(F90_OPTS) main_lapwdft.f90 -o $(EXE) $(OBJ) $(LIB_LPK) $(LIB_FFT) $(LIB_LIBXC) $(LIB_SYS)

clean:
	rm -f *.o *.mod liblapwdft.a liblapwdft.so

mpi_stub.o: mpi_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c mpi_stub.f90

#libxcifc_stub.o: libxcifc_stub.f90
#	$(F90) $(F90_OPTS_NO_WARN) -c libxcifc_stub.f90

libxcifc.o: libxcifc.f90
	$(F90) -I$(LIBXC_HOME)/include $(F90_OPTS_NO_WARN) -c libxcifc.f90

mkl_stub.o: mkl_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c mkl_stub.f90

oblas_stub.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c oblas_stub.f90

blis_stub.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c blis_stub.f90

zfftifc.o: oblas_stub.f90
	$(F90) $(F90_OPTS_NO_WARN) -c zfftifc.f90
