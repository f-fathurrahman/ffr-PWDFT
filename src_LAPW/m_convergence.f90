MODULE m_convergence

!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
INTEGER maxscl
! current self-consistent loop number
INTEGER iscl
! Kohn-Sham potential convergence tolerance
REAL(8) epspot
! energy convergence tolerance
REAL(8) epsengy
! force convergence tolerance
REAL(8) epsforce
! stress tensor convergence tolerance
REAL(8) epsstress

END MODULE 
