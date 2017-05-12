PROGRAM test_calc_rhoe_R
  USE m_PWGrid, ONLY : Ngwx
  USE m_realspace, ONLY : Npoints
  USE m_states, ONLY : Focc, Nstates
  USE m_hamiltonian, ONLY : Rhoe_R
  IMPLICIT NONE 
  REAL(8) :: LL(3,3)
  REAL(8) :: ecutwfc_Ry
  COMPLEX(8), ALLOCATABLE :: v(:,:)

  ecutwfc_Ry = 5.d0
  LL(1,:) = (/ 16.d0, 0.d0, 0.d0 /)
  LL(2,:) = (/ 0.d0, 16.d0, 0.d0 /)
  LL(3,:) = (/ 0.d0, 0.d0, 16.d0 /)

  CALL init_PWGrid( 0.5d0*ecutwfc_Ry, LL )
  CALL info_PWGrid()

  Nstates = 4
  ALLOCATE( v(Ngwx, Nstates) )
  CALL random_wfc( Ngwx, Nstates, v )

  CALL z_ortho_check( Ngwx, Nstates, 1.d0, v )

  ALLOCATE( Focc(Nstates) )
  Focc(:) = 1.d0

  ALLOCATE( Rhoe_R(Npoints) )
  CALL calc_rhoe_R( Focc, v )

  DEALLOCATE( Rhoe_R )
  DEALLOCATE( Focc )
  DEALLOCATE( v )
  CALL dealloc_PWGrid()

END PROGRAM 

