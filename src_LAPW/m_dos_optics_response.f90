module m_dos_optics_response

!----------------------------------------------------------!
!     density of states, optics and response variables     !
!----------------------------------------------------------!
! number of energy intervals in the DOS/optics function plot
integer nwplot
! fine k-point grid size for integration of functions in the Brillouin zone
integer ngrkf
! smoothing level for DOS/optics function plot
integer nswplot
! energy interval for DOS/optics function plot
real(8) wplot(2)
! maximum angular momentum for the partial DOS plot
integer lmaxdos
! dosocc is .true. if the DOS is to be weighted by the occupancy
logical dosocc
! dosmsum is .true. if the partial DOS is to be summed over m
logical dosmsum
! dosssum is .true. if the partial DOS is to be summed over spin
logical dosssum
! number of optical matrix components required
integer noptcomp
! required optical matrix components
integer optcomp(3,27)
! intraband is .true. if the intraband term is to be added to the optical matrix
logical intraband
! lmirep is .true. if the (l,m) band characters should correspond to the
! irreducible representations of the site symmetries
logical lmirep
! spin-quantisation axis in Cartesian coordinates used when plotting the
! spin-resolved DOS (z-axis by default)
real(8) sqados(3)
! q-vector in lattice and Cartesian coordinates for calculating the matrix
! elements < i,k+q | exp(iq.r) | j,k >
real(8) vecql(3),vecqc(3)
! maximum initial-state energy allowed in ELNES transitions
real(8) emaxelnes
! structure factor energy window
real(8) wsfac(2)

end module