SUBROUTINE init0()
  USE m_states, ONLY: stype, swidth, sdescr, efermi
  USE m_density_pot_xc, ONLY: tempk
  USE m_misc, ONLY: tlast  
  USE m_mixing, ONLY: mixtype, mixdescr
  USE m_constants, ONLY: kboltz
  USE m_plotting
  USE m_convergence, ONLY: iscl
  USE m_timing, ONLY: timeinit
  IMPLICIT NONE 
  REAL(8) :: ts0, ts1

  CALL init_timing()
  
  CALL timesec(ts0)

  CALL init_am_variables()
  CALL init_idx_atom_species()
  CALL init_spin_variables()
  
  CALL init_crystal_structure()
  
  CALL init_vector_field_E_A()
  
  CALL symmetry() ! crystal symmetry set up
  
  CALL init_radial_meshes()

  CALL init_charges_states()
  
  CALL init_gvector_arrays()
  
  CALL init_atoms_cores()
  
  CALL init_chgden_pot_xc()
  
  CALL init_forces()

  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  !
  ! determine nuclear radii and volumes
  CALL nuclei()
  !
  ! determine the nuclear-nuclear energy
  CALL energynn()
  !
  ! get smearing function description
  CALL getsdata(stype,sdescr)
  !
  ! get mixing type description
  CALL getmixdata(mixtype,mixdescr)
  !
  ! generate the spherical harmonic transform (SHT) matrices
  CALL genshtmat()
  !
  ! allocate 1D plotting arrays
  IF(allocated(dvp1d)) DEALLOCATE(dvp1d)
  ALLOCATE(dvp1d(nvp1d))
  IF(allocated(vplp1d)) DEALLOCATE(vplp1d)
  ALLOCATE(vplp1d(3,npp1d))
  IF(allocated(dpp1d)) DEALLOCATE(dpp1d)
  ALLOCATE(dpp1d(npp1d))
  !
  ! initial self-consistent loop number
  iscl = 1
  tlast = .false.
  !
  ! set the Fermi energy to zero
  efermi = 0.d0
  !
  ! set the temperature from the smearing width
  tempk = swidth/kboltz

  CALL timesec(ts1)
  CALL timesec(ts1)
  timeinit = timeinit + ts1 - ts0

END SUBROUTINE
