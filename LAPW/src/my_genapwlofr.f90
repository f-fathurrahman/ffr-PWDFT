SUBROUTINE my_genapwlofr()
  IMPLICIT NONE 
  CALL my_genapwfr() ! generate the APW radial functions
  CALL my_genlofr() ! generate the local-orbital radial functions
  CALL my_olprad() ! compute the overlap radial integrals
  CALL my_hmlrad() ! compute the Hamiltonian radial integrals
  RETURN 
END SUBROUTINE 