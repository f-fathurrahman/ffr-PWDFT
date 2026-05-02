!--------------------------
subroutine my_readspecies()
!--------------------------
use modmain
implicit none
! local variables
integer :: is,ist,ios
integer :: nlx,ilx,lx,ilo
integer :: io,jo,ko,l,i,j

e0min = 0.d0

do is = 1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),status='OLD',form='FORMATTED', iostat=ios)
  if (ios /= 0) then
    write(*,*)
    write(*,'("Error(readspecies): error opening species file ",A)') trim(sppath)//trim(spfname(is))
    write(*,*)
    stop
  endif
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) rminsp(is), rmt(is), rmaxsp(is), nrmt(is)
  !
  ! Some sanity checks
  !
  if (rminsp(is) <= 0.d0) then
    write(*,*)
    write(*,'("Error(readspecies): rminsp <= 0 : ",G18.10)') rminsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  !
  if (rmt(is) <= rminsp(is)) then
    write(*,*)
    write(*,'("Error(readspecies): rmt <= rminsp : ",2G18.10)') rmt(is), rminsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  !
  if (rmaxsp(is) < rmt(is)) then
    write(*,*)
    write(*,'("Error(readspecies): rmaxsp < rmt : ",2G18.10)') rmaxsp(is), rmt(is)
    write(*,*)
    stop
  endif
  !
  if (nrmt(is) < 20) then
    write(*,*)
    write(*,'("Error(readspecies): nrmt too small : ",I8)') nrmt(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  !
  ! multiply nrmt by the scale factor
  ! ffr: why?
  nrmt(is) = nint(dble(nrmt(is))*nrmtscf)
  !
  ! reduce the minimum radial mesh point by the same factor
  rminsp(is) = rminsp(is)/nrmtscf
  !
  ! Read states
  read(50,*) nstsp(is)
  if( (nstsp(is) <= 0) .or. (nstsp(is) > maxstsp) ) then
    write(*,*)
    write(*,'("Error(readspecies): nstsp out of range : ",I8)') nstsp(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  !
  do ist = 1,nstsp(is)
    read(50,*) nsp(ist,is), lsp(ist,is), ksp(ist,is), occsp(ist,is), spcore(ist,is)
    !
    ! Some sanity checks for states
    !
    if (nsp(ist,is) < 1) then
      write(*,*)
      write(*,'("Error(readspecies): nsp < 1 : ",I8)') nsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    endif
    !
    if (lsp(ist,is) < 0) then
      write(*,*)
      write(*,'("Error(readspecies): lsp < 0 : ",I8)') lsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    endif
    !
    if (ksp(ist,is) < 1) then
      write(*,*)
      write(*,'("Error(readspecies): ksp < 1 : ",I8)') ksp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    endif
    !
    if (occsp(ist,is) < 0.d0) then
      write(*,*)
      write(*,'("Error(readspecies): occsp < 0 : ",G18.10)') occsp(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    endif
  enddo
  !
  ! APW order (ffr: what's this?) derivative order in APW?
  ! apword is given for each l ?
  !
  read(50,*) apword(0,is)
  if (apword(0,is) <= 0) then
    write(*,*)
    write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  if (apword(0,is) > maxapword) then
    write(*,*)
    write(*,'("Error(readspecies): apword too large : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxapword in modmain and recompile code")')
    write(*,*)
    stop
  endif
  !
  ! set the APW orders for l>0
  !
  apword(1:lmaxapw,is) = apword(0,is) !ffr: the same for all l?
  !
  ! this is loop over apword(0,is). Only for l=0?
  do io = 1,apword(0,is)
    read(50,*) apwe0(io,0,is), apwdm(io,0,is), apwve(io,0,is)
    if (apwdm(io,0,is) < 0) then
      write(*,*)
      write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,0,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    endif
    !
    ! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is) = apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is) = apwdm(io,0,is)
    apwve(io,1:lmaxapw,is) = apwve(io,0,is)
    e0min = min(e0min, apwe0(io,0,is))
  enddo
  !
  ! ffr: What's this? For exception? No example for nlx /= 0 ?
  !
  read(50,*) nlx
  if (nlx < 0) then
    write(*,*)
    write(*,'("Error(readspecies): nlx < 0 : ",I8)') nlx
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  !
  ! ffr: Loop over nlx, no examples?
  !      probably we need to generate our own species file?
  !      in writespecies nlx is always 0 ???
  !
  ! ffr: Additional APW for each l can be specified by this?
  !      By default for APW with l > 0, is set to be the same as l=0.
  !
  do ilx = 1,nlx
    read(50,*) lx, io
    if (lx < 0) then
      write(*,*)
      write(*,'("Error(readspecies): lx < 0 : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    endif
    if (lx > lmaxapw) then
      write(*,*)
      write(*,'("Error(readspecies): lx > lmaxapw : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    endif
    apword(lx,is) = io
    if (apword(lx,is) <= 0) then
      write(*,*)
      write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    endif
    if (apword(lx,is) > maxapword) then
      write(*,*)
      write(*,'("Error(readspecies): apword too large : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,'("Adjust maxapword in modmain and recompile code")')
      write(*,*)
      stop
    endif
    do io = 1,apword(lx,is)
      read(50,*) apwe0(io,lx,is), apwdm(io,lx,is), apwve(io,lx,is)
      if (apwdm(io,lx,is) < 0) then
        write(*,*)
        write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" exception number ",I4)') ilx
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      endif
      e0min = min(e0min,apwe0(io,lx,is))
    enddo
  enddo
  ! end loop over nlx. ffr: this is not required?
  !
  !
  ! add excess order to APW functions if required
  ! ffr: no additional data is read from species file
  if(nxoapwlo > 0) then
    write(*,*) 'Pass here 232 in my_readspecies: nxoapwlo = ', nxoapwlo
    !
    ! loop over lmaxapw ?
    write(*,*) 'lmaxapw = ', lmaxapw
    !
    do l = 0,lmaxapw
      jo = apword(l,is) ! old
      ko = jo + nxoapwlo
      !write(*,'(1x,A,2I4)') 'jo, ko = ', jo, ko
      !
      ! check requested APW order, it cannot be higher than maxapword
      if(ko > maxapword) then
        ko = maxapword
      endif
      i = 0
      ! add additional derivative order: modify apwe0, apwdm, apwve
      ! loop over apw order
      do io = jo+1,ko
        i = i + 1
        apwe0(io,l,is) = apwe0(jo,l,is) ! the same apwe0
        apwdm(io,l,is) = apwdm(jo,l,is) + i ! ffr: This is added
        apwve(io,l,is) = apwve(jo,l,is)
      enddo
      ! update apword
      apword(l,is) = ko
    enddo
    !
    l = 0
    write(*,*) 'apword(l,is) is updated to ', apword(l,is)
    do io = 1,apword(l,is)
      write(*,*) io, apwe0(io,l,is), apwdm(io,l,is), apwve(io,l,is)
    enddo
  endif
  !
  ! Read local orbitals
  !
  read(50,*) nlorb(is)
  if (nlorb(is) < 0) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb < 0 : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  endif
  if (nlorb(is) > maxlorb) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb too large : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxlorb in modmain and recompile code")')
    write(*,*)
    stop
  endif
  do ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is), lorbord(ilo,is)
    if (lorbl(ilo,is) < 0) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl < 0 : ",I8)') lorbl(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    endif
    if (lorbl(ilo,is) > lmaxo) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl > lmaxo : ",2I8)') lorbl(ilo,is), &
       lmaxo
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    endif
    if (lorbord(ilo,is) < 2) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord < 2 : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    endif
    if (lorbord(ilo,is) > maxlorbord) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord too large : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'("Adjust maxlorbord in modmain and recompile code")')
      write(*,*)
      stop
    endif
    do io = 1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      if (lorbdm(io,ilo,is) < 0) then
        write(*,*)
        write(*,'("Error(readspecies): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" local-orbital ",I4)') ilo
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      endif
      e0min = min(e0min, lorbe0(io,ilo,is))
    enddo
  enddo
  !
  ! add excess local-orbitals if required
  ! ffr: why excess?
  if( nxlo > 0 ) then
    lx=-1
    do ilo = 1,nlorb(is)
      do io = 1,lorbord(ilo,is)
        if(lorbe0(io,ilo,is) < 0.d0) then
          goto 10 ! ffr: what's this?
        endif
      enddo
      if(lorbl(ilo,is) > lx) then
        lx = lorbl(ilo,is)
      endif
      10 continue
    enddo
    ilo = nlorb(is)
    do i=1,nxlo
      if (ilo == maxlorb) exit
      l = lx + i
      if (l > lmaxo) exit
      ilo = ilo + 1
      lorbl(ilo,is) = l
      lorbord(ilo,is) = apword(l,is)+1
      do io = 1,lorbord(ilo,is)
        lorbe0(io,ilo,is) = apwe0(1,l,is)
        lorbdm(io,ilo,is) = io - 1
        lorbve(io,ilo,is) = apwve(1,l,is)
      enddo
    enddo
    nlorb(is) = ilo
  endif
  !
  ! add excess order to local-orbitals if required
  if (nxoapwlo > 0) then
    do ilo = 1,nlorb(is)
      ! find the maximum energy derivative
      jo = 1
      j = lorbdm(jo,ilo,is)
      do io = 1,lorbord(ilo,is)
        i = lorbdm(io,ilo,is)
        if (i > j) then
          jo = io
          j = i
        endif
      enddo
      ko = lorbord(ilo,is) + nxoapwlo
      if (ko > maxlorbord) ko = maxlorbord
      i = 0
      do io = lorbord(ilo,is)+1,ko
        i = i + 1
        lorbe0(io,ilo,is)=lorbe0(jo,ilo,is)
        lorbdm(io,ilo,is)=lorbdm(jo,ilo,is)+i
        lorbve(io,ilo,is)=lorbve(jo,ilo,is)
      enddo
      lorbord(ilo,is)=ko
    enddo
  endif
  close(50)
enddo

! set all muffin-tin radii to single value if required
if (rmtall > 0.d0) then
  rmt(1:nspecies) = rmtall
endif

! add conduction state local-orbitals if required
call add_lorb_cnd()

! subtract 2 Hartree from the minimum energy
e0min=e0min-2.d0
return

end subroutine

