PROGRAM test_3d_harmonic
  USE m_PWGrid, ONLY : Ngwx
  USE m_cell, ONLY : CellVolume
  USE m_realspace, ONLY : Npoints
  USE m_states, ONLY : Focc, Nstates, v => KS_evecs
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_constants, ONLY : ZZERO
  IMPLICIT NONE 
  REAL(8) :: LL(3,3)
  REAL(8) :: ecutwfc_Ry
  REAL(8) :: omega
  REAL(8) :: center(3)
  COMPLEX(8), ALLOCATABLE :: Kpsi(:,:), Vpsi(:,:), Hpsi(:,:), grad(:,:)
  INTEGER :: ist

  ecutwfc_Ry = 20.d0
  LL(1,:) = (/ 6.d0, 0.d0, 0.d0 /)
  LL(2,:) = (/ 0.d0, 6.d0, 0.d0 /)
  LL(3,:) = (/ 0.d0, 0.d0, 6.d0 /)

  CALL init_PWGrid( 0.5d0*ecutwfc_Ry, LL )
  CALL info_PWGrid()

  Nstates = 4
  ALLOCATE( v(Ngwx, Nstates) )
  !CALL random_wfc( Ngwx, Nstates, v )
  v(:,:) = ZZERO
  DO ist = 1,Nstates
    v(ist,ist) = cmplx(1.d0,0.d0,kind=8)
  ENDDO 
  CALL z_ortho_check( Ngwx, Nstates, 1.d0, v )
  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  CALL alloc_hamiltonian()

  omega = 2.d0
  center(1) = 0.5d0*sum( LL(:,1) )
  center(2) = 0.5d0*sum( LL(:,2) )
  center(3) = 0.5d0*sum( LL(:,3) )
  !
  CALL init_rgrid()
  CALL init_V_ps_loc_harmonic( omega, center )

  WRITE(*,*) 'V_ps_loc: ', sum(V_ps_loc)*CellVolume/Npoints
  !STOP 

  ALLOCATE( Kpsi(Ngwx,Nstates) )
  ALLOCATE( Vpsi(Ngwx,Nstates) )
  ALLOCATE( Hpsi(Ngwx,Nstates) )
  ALLOCATE( grad(Ngwx,Nstates) )

  CALL op_K( Nstates, v, Kpsi )
  WRITE(*,*) 'sum(Kpsi) = ', sum(Kpsi)
  
  CALL op_V_loc( Nstates, V_ps_loc, v, Vpsi )
  WRITE(*,*) 'sum(Vpsi) = ', sum(Vpsi)
  
  CALL op_H( Nstates, v, Hpsi )
  WRITE(*,*) 'sum(Hpsi) = ', sum(Hpsi)

  CALL calc_grad( Nstates, v, grad )
  WRITE(*,*) 'sum(grad) = ', sum(grad)

  !CALL KS_solve_Emin_pcg( 3.d-5, 100, .FALSE. )

  CALL calc_Rhoe_R( Focc, v )
  CALL calc_energies( v )
  CALL info_energies()

  DEALLOCATE( Focc )
  DEALLOCATE( v )
  CALL dealloc_realspace()
  CALL dealloc_hamiltonian()
  CALL dealloc_PWGrid()

END PROGRAM 

