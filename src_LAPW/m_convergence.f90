module m_convergence

!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
integer maxscl
! current self-consistent loop number
integer iscl
! Kohn-Sham potential convergence tolerance
real(8) epspot
! energy convergence tolerance
real(8) epsengy
! force convergence tolerance
real(8) epsforce
! stress tensor convergence tolerance
real(8) epsstress

end module