!--------------------
subroutine gndstate()
!--------------------
use modmain
use moddftu
use modulr

implicit none
! local variables
logical, parameter :: mp_mpi=.true.
logical exist
integer ik,nwork,n
real(8) dv,etp,de,timetot
! allocatable arrays
real(8), allocatable :: v(:),work(:)

! initialise global variables
call init0()
call init1()

! initialise q-vector-dependent variables if required
if (xctype(1) < 0) call init2()

! apply strain to the G, k, G+k and q-vectors if required
call straingkq()

if (task==0) trdstate=.false.
if (task==1) trdstate=.true.


! write the real and reciprocal lattice vectors to file
call writelat()

! write symmetry matrices to file
call writesym()

! output the k-point set to file
call writekpts()

! write lattice vectors and atomic positions to file
open(50,file='GEOMETRY'//trim(filext),form='FORMATTED')
call writegeom(50)
close(50)

! write interatomic distances to file
open(50,file='IADIST'//trim(filext),form='FORMATTED')
call writeiad(50)
close(50)

! open INFO.OUT file
open(60,file='INFO'//trim(filext),form='FORMATTED')

! open TOTENERGY.OUT
open(61,file='TOTENERGY'//trim(filext),form='FORMATTED')

! open FERMIDOS.OUT
open(62,file='FERMIDOS'//trim(filext),form='FORMATTED')

! open MOMENT.OUT if required
if (spinpol) then
  open(63,file='MOMENT'//trim(filext),form='FORMATTED')
endif

! open GAP.OUT
open(64,file='GAP'//trim(filext),form='FORMATTED')

! open RMSDVS.OUT
open(65,file='RMSDVS'//trim(filext),form='FORMATTED')

! open DTOTENERGY.OUT
open(66,file='DTOTENERGY'//trim(filext),form='FORMATTED')


! open TMDFTU.OUT
if( tmwrite ) then
  open(67,file='TMDFTU'//trim(filext),form='FORMATTED')
endif

! open MOMENTM.OUT
if(spinpol) then
  open(68,file='MOMENTM'//trim(filext),form='FORMATTED')
endif

! open RESIDUAL.OUT
if( xctype(1) .lt. 0) then
  open(69,file='RESIDUAL'//trim(filext),form='FORMATTED')
endif

! write out general information to INFO.OUT
call writeinfo(60)
write(60,*)


iscl=0
if (trdstate) then
  stop 'trdstate is not supported'
  ! read the Kohn-Sham potential and fields from file
  !call readstate()
  !write(60,'("Potential read in from STATE.OUT")')
  !if (autolinengy) call readfermi()
else
  ! initialise the density and magnetisation from atomic data
  call rhoinit()
  call maginit()
  ! generate the Kohn-Sham potential and magnetic field
  call potks(.true.)
  write(60,'("Kohn-Sham potential initialised from atomic data")')
endif

flush(60)


call genvsig()

! size of mixing vector
n = npmtmax*natmtot+ngtot
if (spinpol) n=n+(npcmtmax*natmtot+ngtot)*ndmag
if (tvmatmt) n=n+2*((lmmaxdm*nspinor)**2)*natmtot

! allocate mixing array
allocate(v(n))

! determine the size of the mixer work array
nwork = -1
call mixerifc(mixtype,n,v,dv,nwork,v)
allocate(work(nwork))

! initialise the mixer
iscl = 0
call mixpack(.true.,n,v)
call mixerifc(mixtype,n,v,dv,nwork,work)

! set the stop signal to .false.
tstop = .false.

! set last self-consistent loop flag
tlast = .false.
etp = 0.d0

! begin the self-consistent loop
write(60,*)
write(60,'("+------------------------------+")')
write(60,'("| Self-consistent loop started |")')
write(60,'("+------------------------------+")')



do iscl=1,maxscl
  
  write(60,*)
  write(60,'("+--------------------+")')
  write(60,'("| Loop number : ",I4," |")') iscl
  write(60,'("+--------------------+")')

  if (iscl >= maxscl) then
    write(60,*)
    write(60,'("Reached self-consistent loops maximum")')
    write(*,*)
    write(*,'("Warning(gndstate): failed to reach self-consistency after ",I4," loops")') iscl
    tlast=.true.
  endif

  ! generate the core wavefunctions and densities
  call gencore()
  
  ! find the new linearisation energies
  call linengy()
  
  ! write out the linearisation energies
  call writelinen()
  
  ! generate the APW and local-orbi\tal radial functions and integrals
  call genapwlofr()
  
  ! generate the spin-orbit coupling radial functions
  call gensocfr()
  
  ! generate the first- and second-variational eigenvectors and eigenvalues
  call genevfsv()
  
  ! find the occupation numbers and Fermi energy
  call occupy()
  !
  if( autoswidth ) then
    write(60,*)
    write(60,'("New smearing width : ",G18.10)') swidth
  endif
  ! write the occupation numbers to file
  do ik = 1,nkpt
      call putoccsv(filext,ik,occsv(:,ik))
  enddo

  ! write eigenvalues to file
  call writeeval()
  
  ! write the Fermi energy to file
  call writefermi()

  ! generate the density and magnetisation
  call rhomag()
  
  ! DFT+U or fixed tensor moment calculation
  if( (dftu /= 0) .or. (ftmtype /= 0) ) then
    !
    ! generate the muffin-tin density matrix used for computing the potential matrix
    call gendmatmt()
    !
    ! write the FTM tensor moments to file
    if(ftmtype /= 0) call writeftm()
    !
    ! generate the DFT+U or FTM muffin-tin potential matrices
    call genvmatmt()
  endif
  
  ! ffr: this is disabled
  if(dftu /= 0) then
    ! write the DFT+U matrices to file
    !call writedftu()
    !! calculate and write tensor moments to file
    !if (tmwrite) then
    !  if (spinorb) then
    !    call writetm3du(67)
    !  else
    !    call writetm2du(67)
    !  endif
    !endif
  endif

  ! compute the Kohn-Sham potentials and magnetic fields
  call potks(.true.)

  if( (xcgrad == 3) .and. (c_tb09 /= 0.d0)) then
      write(60,*)
      write(60,'("Tran-Blaha ''09 constant c : ",G18.10)') c_tb09
  endif

  ! pack interstitial and muffin-tin potential and field into one array
  call mixpack(.true.,n,v)

  ! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)

  ! unpack potential and field
  call mixpack(.false.,n,v)
  
  ! calculate and add the fixed spin moment effective field (after mixing)
  call fsmbfield()
  call addbfsm()

  ! Fourier transform Kohn-Sham potential to G-space
  call genvsig()
  
  ! reduce the external magnetic fields if required
  if (reducebf.lt.1.d0) then
    bfieldc(:)=bfieldc(:)*reducebf
    bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
  endif

  ! compute the energy components
  call energy()


  ! output energy components
  call writeengy(60)
  write(60,*)
  write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
  write(60,'(" (states/Hartree/unit cell)")')
  write(60,*)
  write(60,'("Estimated indirect band gap : ",G18.10)') bandgap(1)
  write(60,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
  write(60,'("Estimated direct band gap   : ",G18.10)') bandgap(2)
  write(60,'(" at k-point ",I6)') ikgap(3)
  ! write total energy to TOTENERGY.OUT
  write(61,'(G22.12)') engytot
  flush(61)
  ! write DOS at Fermi energy to FERMIDOS.OUT
  write(62,'(G18.10)') fermidos
  flush(62)
  ! output charges and moments
  call writechg(60)
  if (spinpol) then
    call writemom(60)
    ! write total moment to MOMENT.OUT
    write(63,'(3G18.10)') momtot(1:ndmag)
    flush(63)
    ! write total moment magnitude to MOMENTM.OUT
    write(68,'(G18.10)') momtotm
    flush(68)
  endif
  ! write estimated Kohn-Sham indirect band gap
  write(64,'(G22.12)') bandgap(1)
  flush(64)
  ! output effective fields for fixed spin moment calculations
  if (fsmtype /= 0) call writefsm(60)
  ! check for WRITE file
  inquire(file='WRITE',exist=exist)
  if (exist) then
    write(60,*)
    write(60,'("WRITE file exists - writing STATE.OUT")')
    call writestate()
    open(50,file='WRITE')
    close(50,status='DELETE')
  endif
  ! write STATE.OUT file if required
  if( nwrite >= 1 ) then
    if (mod(iscl,nwrite) == 0) then
      call writestate()
      write(60,*)
      write(60,'("Wrote STATE.OUT")')
    endif
  endif
  ! write OEP residual
  if (xctype(1) < 0) then
    write(60,*)
    write(60,'("Magnitude of OEP residual : ",G18.10)') resoep
    write(69,'(G18.10)') resoep
    flush(69)
  endif
  

  ! exit self-consistent loop if required
  if (tlast) goto 10
  
  ! check for convergence
  if(iscl >= 2) then
    write(60,*)
    write(60,'("RMS change in Kohn-Sham potential (target) : ",G18.10," (",G18.10,")")') dv,epspot
    write(65,'(G18.10)') dv
    flush(65)
    de = abs(engytot-etp)
    !
    write(60,'("Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') de,epsengy
    write(66,'(G18.10)') de
    flush(66)
    !
    if ((dv.lt.epspot).and.(de.lt.epsengy)) then
      write(60,*)
      write(60,'("Convergence targets achieved")')
      tlast=.true.
    endif
  endif
  ! average the current and previous total energies and store
  if (iscl.gt.1) then
    etp=0.75d0*engytot+0.25d0*etp
  else
    etp=engytot
  endif
  ! check for STOP file (only master process)
  inquire(file='STOP',exist=exist)
  if (exist) then
    write(60,*)
    write(60,'("STOP file exists - stopping self-consistent loop")')
    open(50,file='STOP')
    close(50,status='DELETE')
    tstop=.true.
    tlast=.true.
  endif

  ! output the current total CPU time
  timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
  write(60,*)
  write(60,'("Time (CPU seconds) : ",F12.2)') timetot


! end the self-consistent loop
enddo


10 continue

write(60,*)
write(60,'("+------------------------------+")')
write(60,'("| Self-consistent loop stopped |")')
write(60,'("+------------------------------+")')
! write density and potentials to file only if maxscl > 1
if (maxscl .gt. 1) then
  call writestate()
  write(60,*)
  write(60,'("Wrote STATE.OUT")')
endif

! compute forces if required
if (tforce) then
  call force()
  ! output forces to INFO.OUT
  call writeforces(60)
endif

! compute the current density and total current if required
if (tcden) then
  call curden(afieldc)
  write(60,*)
  write(60,'("Total current per unit cell")')
  write(60,'(3G18.10)') curtot
  write(60,'(" magnitude : ",G18.10)') curtotm
endif

! total time used
timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor

! output timing information
write(60,*)
write(60,'("Timings (CPU seconds) :")')
write(60,'(" initialisation",T40,": ",F12.2)') timeinit
write(60,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
write(60,'(" first-variational eigenvalue equation",T40,": ",F12.2)') timefv
if (tevecsv) then
  write(60,'(" second-variational calculation",T40,": ",F12.2)') timesv
endif
write(60,'(" charge density calculation",T40,": ",F12.2)') timerho
write(60,'(" potential calculation",T40,": ",F12.2)') timepot
if (tforce) then
  write(60,'(" force calculation",T40,": ",F12.2)') timefor
endif
write(60,'(" total",T40,": ",F12.2)') timetot
write(60,*)
write(60,'("+----------------------------+")')
write(60,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
write(60,'("+----------------------------+")')
! close the INFO.OUT file
close(60)
! close the TOTENERGY.OUT file
close(61)
! close the FERMIDOS.OUT file
close(62)
! close the MOMENT.OUT and MOMENTM.OUT files
if (spinpol) then
  close(63); close(68)
endif
! close the GAP.OUT file
close(64)
! close the RMSDVS.OUT file
close(65)
! close the DTOTENERGY.OUT file
close(66)
! close TMDFTU.OUT file
if (tmwrite) close(67)

! close the RESIDUAL.OUT file
if (xctype(1) < 0) close(69)

deallocate(v,work)

return

end subroutine

