SUBROUTINE info_apwlo()
  USE m_atoms
  USE m_atomic_species
  USE m_apwlo
  USE m_muffin_tins
  IMPLICIT NONE 

  integer :: is, l, iord, ia, ias

  write(*,*)
  write(*,*) '----------'
  write(*,*) 'Info APWLO'
  write(*,*) '----------'

  WRITE(*,*) 
  WRITE(*,'(1x,A,F18.10)') 'deapwlo (energy step used for numerical calculation of energy derivatives) = ', deapwlo

  WRITE(*,*)
  WRITE(*,'(1x,A,I2)') 'maxapword (maximum allowable APW order (parameter) = ', 4

  WRITE(*,*)
  WRITE(*,*) 'apword (APW order) for each species and angular momentum'
  DO is = 1,nspecies
    WRITE(*,*) 'Species: ', trim(spfname(is))
    DO l = 0,lmaxapw
      WRITE(*,'(1x,I3,I4)') l, apword(l,is)
    ENDDO
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'maximum of apword over all angular momenta and species'
  WRITE(*,*) 'apwordmax = ', apwordmax

  write(*,*)
  write(*,*) 'total number of APW coefficients (l, m and order) for each species'
  do is = 1,nspecies
    write(*,*) 'lmoapw = ', lmoapw(is)
  enddo

  write(*,*)
  write(*,*) 'polynomial order used for APW radial derivatives'
  write(*,*) 'npapw = ', npapw
  
  write(*,*)
  write(*,*) 'apwe0 (APW initial linearization energies)'
  do is = 1,Nspecies
    write(*,*) 'Species ', trim(spfname(is))
    do l = 0,lmaxapw
      do iord = 1,apword(l,is)
        write(*,'(1x,2I3,F18.10)') iord, l, apwe0(iord,l,is)
      enddo
    enddo
  enddo

  write(*,*)
  write(*,*) 'apwe (APW linearization energies)'
  do is = 1,nspecies
    write(*,*) 'Species ', trim(spfname(is))
    do ia = 1,natoms(is)
      ias = idxas(ia,is)
      write(*,*) 'ias = ', ias
      do l = 0,lmaxapw
        do iord = 1,apword(l,is)
          write(*,'(1x,2I3,F18.10)') iord, l, apwe(iord,l,ias)
        enddo
      enddo
    enddo
  enddo

  write(*,*)
  write(*,*) 'APW derivative order'
  do is = 1,nspecies
    write(*,*) 'Species ', trim(spfname(is))
    do l = 0,lmaxapw
      do iord = 1,apword(l,is)
        write(*,*) 'apwdm = ', apwdm(iord,l,is)
      enddo
    enddo
  enddo

  write(*,*)
  write(*,*) 'apwve is .true. if the linearisation energies are allowed to vary'
  do is = 1,nspecies
    write(*,*) 'Species ', trim(spfname(is))
    do l = 0,lmaxapw
      do iord = 1,apword(l,is)
        write(*,*) 'apwdm = ', apwve(iord,l,is)
      enddo
    enddo
  enddo

  write(*,*)
  write(*,*) 'APW radial functions dimensions = '
  !allocate(apwfr(nrmtmax,2,apwordmax,0:lmaxapw,natmtot))
  write(*,*) 'shape apwfr = ', shape(apwfr)
  write(*,*) 'nrmtmax = ', nrmtmax
  write(*,*) 'nrmt = ', nrmt
  !real(8), allocatable :: apwfr(:,:,:,:,:)

  write(*,*) 'nlorb = ', nlorb
  write(*,*) 'nlotot = ', nlotot

  !! derivate of radial functions at the muffin-tin surface
!real(8), allocatable :: apwdfr(:,:,:)
!! maximum number of local-orbitals
!integer, parameter :: maxlorb=200
!! maximum allowable local-orbital order
!integer, parameter :: maxlorbord=5
!! number of local-orbitals
!integer nlorb(maxspecies)
!! maximum nlorb over all species
!integer nlomax
!! total number of local-orbitals
!integer nlotot
!! local-orbital order
!integer lorbord(maxlorb,maxspecies)
!! maximum lorbord over all species
!integer lorbordmax
!! polynomial order used for local-orbital radial derivatives
!integer nplorb
!! local-orbital angular momentum
!integer lorbl(maxlorb,maxspecies)
!! maximum lorbl over all species
!integer lolmax
!! (lolmax+1)^2
!integer lolmmax
!! local-orbital initial energies
!real(8) lorbe0(maxlorbord,maxlorb,maxspecies)
!! local-orbital energies
!real(8), allocatable :: lorbe(:,:,:)
!! local-orbital derivative order
!integer lorbdm(maxlorbord,maxlorb,maxspecies)
!! lorbve is .true. if the linearisation energies are allowed to vary
!logical lorbve(maxlorbord,maxlorb,maxspecies)
!! local-orbital radial functions
!real(8), allocatable :: lofr(:,:,:,:)
!! band energy search tolerance
!real(8) epsband
!! maximum allowed change in energy during band energy search; enforced only if
!! default energy is less than zero
!real(8) demaxbnd
!! minimum default linearisation energy over all APWs and local-orbitals
!real(8) e0min
!! if autolinengy is .true. then the fixed linearisation energies are set to the
!! Fermi energy minus dlefe
!logical autolinengy
!! difference between linearisation and Fermi energies when autolinengy is .true.
!real(8) dlefe
!! lorbcnd is .true. if conduction state local-orbitals should be added
!logical lorbcnd
!! conduction state local-orbital order
!integer lorbordc
!! excess order of the APW and local-orbital functions
!integer nxoapwlo
!! excess local orbitals
!integer nxlo

end subroutine
