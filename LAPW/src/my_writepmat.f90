subroutine my_writepmat
use modmain
implicit none

! initialise universal variables
call init0
call init1

! read in the density and potentials from file
call readstate

! read Fermi energy from file
call readfermi

! find the new linearisation energies
call linengy

! generate the APW radial functions
call genapwfr

! generate the local-orbital radial functions
call genlofr

! write the momentum matrix elements in the second-variational basis to file
call my_genpmat()

write(*,*)
write(*,'("Info(writepmat):")')
write(*,'(" momentum matrix elements written to file PMAT.OUT")')

return
end subroutine
