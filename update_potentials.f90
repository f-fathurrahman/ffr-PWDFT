SUBROUTINE update_potentials()

  USE m_realspace, ONLY : Npoints
  USE m_hamiltonian, ONLY : Rhoe_R, V_Hartree, V_xc

  IMPLICIT NONE 
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  
  ALLOCATE( epsxc(Npoints) )
  ALLOCATE( depsxc(Npoints) )

  CALL Poisson_solve_fft( Rhoe_R, V_Hartree )

  CALL excVWN( Npoints, Rhoe_R, epsxc )
  CALL excpVWN( Npoints, Rhoe_R, depsxc )

  V_xc(:) = epsxc(:) + Rhoe_R(:)*depsxc(:)

  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END SUBROUTINE 
