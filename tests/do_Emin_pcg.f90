PROGRAM do_Emin_pcg

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_states, ONLY : Nstates, Focc, &
                       evals => KS_evals, &
                       evecs => KS_evecs
  USE m_PWGrid, ONLY : Ngwx

  IMPLICIT NONE 
  INTEGER :: Narg
  CHARACTER(64) :: filexyz, arg_1
  REAL(8) :: ecutwfc_Ry
  REAL(8) :: LL(3,3)
  INTEGER :: iargc  ! pgf90 
  INTEGER :: tstart, tstop, clocks_per_second

  CALL system_clock( tstart, clocks_per_second )

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given:'
    WRITE(*,*) '       ecutwfc_Ry and path to structure file'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_1 )
  READ(arg_1, *) ecutwfc_Ry

  CALL getarg( 2, filexyz )

  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../pseudopotentials/pade_gth/'

  CALL init_PsPot()

  LL(1,:) = (/ 16.d0, 0.d0, 0.d0 /)
  LL(2,:) = (/ 0.d0, 16.d0, 0.d0 /)
  LL(3,:) = (/ 0.d0, 0.d0, 16.d0 /)

  CALL init_PWGrid( 0.5d0*ecutwfc_Ry, LL )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_PWGrid()

!  CALL init_betaNL()

  ! Initialize occupation numbers
  CALL init_states()

  ! Structure factor, shifted to FFT grid
  CALL init_strfact()

  ! Ewald energy
  CALL calc_Ewald()

  ! Memory for potentials
  CALL alloc_hamiltonian()

  ! Local pseudopotential
  CALL init_V_ps_loc_G()

  ! Manually allocate KS eigenvectors and eigenvalues
  ALLOCATE( evecs(Ngwx,Nstates), evals(Nstates) )

  CALL random_wfc( Ngwx, Nstates, evecs )
  CALL z_ortho_check( Ngwx, Nstates, 1.d0, evecs )

  CALL KS_solve_Emin_pcg( 3.d-5, 1000, .FALSE. )
  !CALL KS_solve_Emin_pcg( 3.d-5, 1000, .TRUE. )

  CALL info_energies()

  !
  DEALLOCATE( evecs, evals )
  DEALLOCATE( Focc )

  CALL dealloc_hamiltonian()
  CALL dealloc_PWGrid()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

  CALL system_clock( tstop )

  WRITE(*,*)
  WRITE(*,*) 'Total elapsed time: ', dble(tstop-tstart)/clocks_per_second, ' seconds.'
  WRITE(*,*)

END PROGRAM

