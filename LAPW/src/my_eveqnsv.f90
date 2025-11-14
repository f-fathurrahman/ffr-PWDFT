
subroutine my_eveqnsv(ik,ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
  use modmain
  use moddftu
  !
  implicit none
  ! arguments
  integer, intent(in) :: ik ! added by ffr
  integer, intent(in) :: ngp,igpig(ngkmax)
  real(8), intent(in) :: vgpc(3,ngkmax)
  complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8), intent(in) :: evalfv(nstfv)
  complex(8), intent(in) :: evecfv(nmatmax,nstfv)
  real(8), intent(out) :: evalsvp(nstsv)
  complex(8), intent(out) :: evecsv(nstsv,nstsv)
  !
  ! local variables
  logical socz
  integer nsc,nsd,ld,ist,jst
  integer ispn,jspn,is,ias
  integer nrc,nrci,nrco,irco,irc
  integer l,lm,nm,npc,npci
  integer igp,i,j,k
  real(8) ca,t1
  complex(8) z1
  ! allocatable arrays
  complex(8), allocatable :: wfmt1(:,:),wfmt2(:),wfmt3(:),wfmt4(:,:)
  complex(8), allocatable :: zlflm(:,:),gwfmt(:,:,:),gzfmt(:,:)
  complex(8), allocatable :: wfir1(:),wfir2(:),wfgp(:,:),gwfgp(:,:)
  ! external functions
  complex(8) zdotc,zfmtinp
  external zdotc,zfmtinp
  
  write(*,*)
  write(*,*) '<div> ENTER my_eveqnsv'
  write(*,*)

  !
  ! no calculation of second-variational eigenvectors
  if(.not. tevecsv) then
    !
    write(*,*) '*** No need for second-variational eigenvectors'
    ! 
    ! evalsvp = evalfv
    do i = 1,nstsv
      evalsvp(i) = evalfv(i)
    enddo
    !
    ! evecsv is identity matrix
    evecsv(:,:) = 0.d0
    do i=1,nstsv
      evecsv(i,i) = 1.d0
    enddo
    return
  endif
  !
  ! coupling constant of the external A-field (1/c)
  ca = 1.d0/solsc
  ! number of spin combinations after application of Hamiltonian
  if( spinpol ) then
    if (ncmag .or. spinorb) then
      nsc = 3
    else
      nsc = 2
    endif
    nsd = 2
  else
    ! XXX ffr: In what case we reach here?
    nsc = 1
    nsd = 1
  endif
  !
  ! special case of spin-orbit coupling and collinear magnetism
  if( spinorb .and. cmagz ) then
    socz = .true.
  else
    socz = .false.
  endif
  write(*,*) 'nsc = ', nsc, ' nsd = ', nsd, ' socz = ', socz

  !
  ld = lmmaxdm*nspinor
  !
  ! zero the second-variational Hamiltonian (stored in the eigenvector array)
  evecsv(:,:) = 0.d0
  !
  !-------------------------!
  !     muffin-tin part     !
  !-------------------------!
  allocate(wfmt1(npcmtmax,nstfv))
  !
  if (xcgrad == 4) allocate(gwfmt(npcmtmax,3,nstfv))

  allocate(wfmt2(npcmtmax), wfmt3(npcmtmax), wfmt4(npcmtmax,nsc))
  if (spinorb) allocate(zlflm(lmmaxo,3))
  if (tafield .or. (xcgrad == 4)) allocate(gzfmt(npcmtmax,3))

  wfmt2(:) = 0.d0
  wfmt3(:) = 0.d0
  wfmt4(:,:) = 0.d0

  !
  ! begin loop over atoms
  do ias = 1,natmtot
    is = idxis(ias)
    nrc = nrcmt(is)
    nrci = nrcmti(is)
    nrco = nrc - nrci
    irco = nrci + 1
    npc = npcmt(is)
    npci = npcmti(is)
    !
    ! compute the first-variational wavefunctions
    do ist = 1,nstfv
      write(*,*) 'sum evecfv(1:nmat(1,ik),ist) = ', sum(evecfv(1:nmat(1,ik),ist))
      call my_wavefmt( lradstp, ias, ngp, apwalm(:,:,:,ias), evecfv(:,ist), wfmt1(:,ist) )
      write(*,*) 'sum wfmt1(1:npc,ist) = ', sum(wfmt1(1:npc,ist))
    enddo
    write(*,*) 'sum wfmt1 = ', sum(wfmt1(1:npc,:))
    write(*,*) 'sum bsmt = ', sum(bsmt(1:npc,ias,ndmag))
    !
    ! begin loop over states
    do jst = 1,nstfv
      if (spinpol) then
        !
        ! convert wavefunction to spherical coordinates
        write(*,*) 'sum wfmt1(:,jst) = ', sum(wfmt1(:,jst))
        call zbsht(nrc, nrci, wfmt1(:,jst), wfmt2)
        write(*,*) 'sum wfmt2 = ', sum(wfmt2)
        !write(1111,*) real(wfmt1(:,jst))
        !write(2222,*) aimag(wfmt1(:,jst))
        !write(3333,*) real(wfmt2)
        !write(4444,*) aimag(wfmt2)
        !stop 'DEBUG 134'
        !
        ! apply Kohn-Sham effective magnetic field
        wfmt3(1:npc) = bsmt(1:npc,ias,ndmag)*wfmt2(1:npc)  ! ffr: THIS IS IMPORTANT
        write(*,*) 'sum wfmt3 = ', sum(wfmt3)
        !
        ! convert to spherical harmonics and store in wfmt4
        call zfsht(nrc, nrci, wfmt3, wfmt4)
        wfmt4(1:npc,2) = -wfmt4(1:npc,1)
        write(*,*) 'sum wfmt4(:,1) = ', sum(wfmt4(:,1))
        write(*,*) 'sum wfmt4(:,2) = ', sum(wfmt4(:,2))
        !
        ! non-collinear magnetic field
        if (ncmag) then
          !write(*,*) 'ncmag is applied'
          wfmt3(1:npc) = cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
          call zfsht(nrc, nrci, wfmt3, wfmt4(:,3))
        endif
        if( socz ) then
          wfmt4(1:npc,3) = 0.d0
        endif
        !
        ! apply spin-orbit coupling if required
        if (spinorb) then
          !write(*,*) 'spinorb is applied'
          ! inner part of muffin-tin
          i = 1
          do irc = 1,nrci
            call lopzflm(lmaxi,wfmt1(i,jst),lmmaxo,zlflm)
            t1 = socfr(irc,ias)
            do j = 1,lmmaxi
              wfmt4(i,1) = wfmt4(i,1) + t1*zlflm(j,3)
              wfmt4(i,2) = wfmt4(i,2) - t1*zlflm(j,3)
              wfmt4(i,3) = wfmt4(i,3) + t1*(zlflm(j,1) + cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
              i = i + 1
            enddo
          enddo
          !
          ! outer part of muffin-tin
          do irc = irco,nrc
            call lopzflm(lmaxo,wfmt1(i,jst),lmmaxo,zlflm)
            t1=socfr(irc,ias)
            do j=1,lmmaxo
              wfmt4(i,1)=wfmt4(i,1)+t1*zlflm(j,3)
              wfmt4(i,2)=wfmt4(i,2)-t1*zlflm(j,3)
              wfmt4(i,3)=wfmt4(i,3)+t1*(zlflm(j,1) &
               +cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
              i=i+1
            enddo
          enddo
        endif
      else
        ! No SOC, no spinpol
        do k=1,nsc
          wfmt4(1:npc,k) = 0.d0
        enddo
      endif
      !
      ! apply muffin-tin potential matrix if required
      if( tvmatmt ) then ! defined in case of DFT+U
        !write(*,*) 'muffin-tin potential matrix is applied'
        do l=0,lmaxdm
          if (tvmmt(l,ias)) then
            nm=2*l+1
            lm=idxlm(l,-l)
            do k=1,nsc
              if (k == 1) then
                ispn=1
                jspn=1
              else if (k == 2) then
                ispn=2
                jspn=2
              else
                ispn=1
                jspn=2
              endif
              if (l.le.lmaxi) then
                call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,ispn,lm,jspn,ias), &
                 ld,wfmt1(lm,jst),lmmaxi,zone,wfmt4(lm,k),lmmaxi)
              endif
              i=npci+lm
              call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
               wfmt1(i,jst),lmmaxo,zone,wfmt4(i,k),lmmaxo)
            enddo
          endif
        enddo
      endif
      !
      ! apply vector potential if required
      if (tafield) then
        call gradzfmt( nrc, nrci, rlcmt(:,1,is), rlcmt(:,-1,is), wfmt1(:,jst), &
                       npcmtmax, gzfmt)
        do i=1,npc
          z1=afieldc(1)*gzfmt(i,1)+afieldc(2)*gzfmt(i,2)+afieldc(3)*gzfmt(i,3)
          z1=ca*cmplx(-aimag(z1),dble(z1),8)
          wfmt4(i,1:nsd)=wfmt4(i,1:nsd)+z1
        enddo
      endif
      !
      ! add to second-variational Hamiltonian matrix
      !
      ! upper diagonal block
      do ist = 1,jst
        z1 = zfmtinp(nrc, nrci, wrcmt(:,is), wfmt1(:,ist), wfmt4)
        write(*,*) 'Upper diagonal block: ist = ', ist, ' z1 from zfmtinp = ', z1
        evecsv(ist,jst) = evecsv(ist,jst) + z1
      enddo
      !
      ! lower diagonal block
      if( nsc >= 2 ) then ! spinpol or ncmag or SOC
        j = jst + nstfv
        do ist = 1,jst
          i = ist + nstfv
          z1 = zfmtinp(nrc, nrci, wrcmt(:,is), wfmt1(:,ist), wfmt4(:,2))
          evecsv(i,j) = evecsv(i,j) + z1
        enddo
      endif
      !
      ! off-diagonal block
      if (nsc == 3) then ! in case of ncmag .or. spinorb
        do ist = 1,nstfv
          z1 = zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist),wfmt4(:,3))
          evecsv(ist,j) = evecsv(ist,j) + z1
        enddo
      endif
      ! end loop over states
    enddo
    !
    ! apply tau-DFT non-multiplicative potential if required
    if (xcgrad == 4) then
      do ist = 1,nstfv
        call gradzfmt( nrc, nrci, rlcmt(:,1,is),rlcmt(:,-1,is),wfmt1(:,ist), &
                       npcmtmax, gwfmt(:,:,ist))
      enddo
      do jst=1,nstfv
        do k=1,3
          call zbsht(nrc,nrci,gwfmt(:,k,jst),wfmt2)
          wfmt2(1:npc)=wsmt(1:npc,ias)*wfmt2(1:npc)
          call zfsht(nrc,nrci,wfmt2,gzfmt(:,k))
        enddo
        do ist=1,nstfv
          z1=0.d0
          do k=1,3
            z1=z1+zfmtinp(nrc,nrci,wrcmt(:,is),gwfmt(:,k,ist),gzfmt(:,k))
          enddo
          z1=0.5d0*z1
          evecsv(ist,jst)=evecsv(ist,jst)+z1
          if (spinpol) then
            i=ist+nstfv
            j=jst+nstfv
            evecsv(i,j)=evecsv(i,j)+z1
          endif
        enddo
      enddo
    endif
    !
  enddo     ! end loop over atoms
  deallocate(wfmt2,wfmt3,wfmt4)
  if (spinorb) deallocate(zlflm)
  if (tafield.or.(xcgrad == 4)) deallocate(gzfmt)
  !
  deallocate(wfmt1)
  if (xcgrad == 4) deallocate(gwfmt)
  

  !---------------------------!
  !     interstitial part     !
  !---------------------------!
  if (spinpol .or. tafield .or. (xcgrad == 4)) then
    if (socz) nsc = 2
    if (xcgrad == 4) allocate(gwfgp(ngkmax,nstfv))
    !
    allocate(wfir1(ngtot),wfir2(ngtot),wfgp(ngkmax,nsc))
    !
    ! begin loop over states
    do jst = 1,nstfv
      wfir1(:)=0.d0
      do igp=1,ngp
        wfir1(igfft(igpig(igp)))=evecfv(igp,jst)
      enddo
      !
      ! Fourier transform wavefunction to real-space
      call zfftifc(3,ngridg,1,wfir1)
      !
      ! multiply with magnetic field and transform to G-space
      if( spinpol ) then
        wfir2(:) = bsir(:,ndmag)*wfir1(:) ! IMPORTANT: multiply in real-space
        call zfftifc(3, ngridg, -1, wfir2)
        do igp = 1,ngp
          wfgp(igp,1) = wfir2(igfft(igpig(igp)))
        enddo
        wfgp(1:ngp,2) = -wfgp(1:ngp,1)
        !
        if( ncmag ) then
          wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
          call zfftifc(3,ngridg,-1,wfir2)
          do igp=1,ngp
            wfgp(igp,3)=wfir2(igfft(igpig(igp)))
          enddo
        endif
      else
        wfgp(1:ngp,1:nsd)=0.d0
      endif
      !
      ! apply vector potential if required
      if (tafield) then
        wfir1(:) = 0.d0
        do igp = 1,ngp
          t1 = -ca*dot_product(afieldc(:),vgpc(:,igp))
          wfir1(igfft(igpig(igp))) = t1*evecfv(igp,jst)
        enddo
        call zfftifc(3,ngridg,1,wfir1)
        wfir1(:) = wfir1(:)*cfunir(:)
        call zfftifc(3,ngridg,-1,wfir1)
        do igp=1,ngp
          z1 = wfir1(igfft(igpig(igp)))
          wfgp(igp,1:nsd) = wfgp(igp,1:nsd) + z1
        enddo
      endif
      !
      ! add to second-variational Hamiltonian matrix
      !
      ! upper diagonal block
      do ist=1,jst
        evecsv(ist,jst) = evecsv(ist,jst) + zdotc(ngp, evecfv(:,ist), 1, wfgp, 1)
      enddo
      !
      ! lower diagonal block
      if( nsc >= 2) then
        j = jst + nstfv
        do ist = 1,jst
          i = ist + nstfv
          evecsv(i,j) = evecsv(i,j) + zdotc(ngp, evecfv(:,ist), 1, wfgp(:,2), 1)
        enddo
      endif
      !
      ! off-diagonal block
      if (nsc == 3) then
        do ist=1,nstfv
          evecsv(ist,j) = evecsv(ist,j) + zdotc(ngp, evecfv(:,ist), 1, wfgp(:,3), 1)
        enddo
      endif
      !
    enddo    ! end loop over states nstfv
    !
    ! apply tau-DFT non-multiplicative potential if required
    if (xcgrad == 4) then
      do k = 1,3
        ! determine the gradient of the wavefunctions
        do ist=1,nstfv
          do igp=1,ngp
            z1=evecfv(igp,ist)
            gwfgp(igp,ist)=vgpc(k,igp)*cmplx(-aimag(z1),dble(z1),8)
          enddo
        enddo
        do jst=1,nstfv
          wfir1(:)=0.d0
          do igp=1,ngp
            wfir1(igfft(igpig(igp)))=gwfgp(igp,jst)
          enddo
          call zfftifc(3,ngridg,1,wfir1)
          wfir1(:)=wsir(:)*wfir1(:)
          call zfftifc(3,ngridg,-1,wfir1)
          do igp=1,ngp
            wfgp(igp,1)=wfir1(igfft(igpig(igp)))
          enddo
          do ist=1,nstfv
            z1=0.5d0*zdotc(ngp,gwfgp(:,ist),1,wfgp,1)
            evecsv(ist,jst)=evecsv(ist,jst)+z1
            if (spinpol) then
              i=ist+nstfv
              j=jst+nstfv
              evecsv(i,j)=evecsv(i,j)+z1
            endif
          enddo
        enddo ! nstfv
      enddo ! k
    endif
    deallocate(wfir1, wfir2, wfgp)
    !
    if (xcgrad == 4) deallocate(gwfgp)
  endif
  !
  ! add the diagonal first-variational part
  i = 0
  do ispn = 1,nspinor
    do ist = 1,nstfv
      i = i + 1
      evecsv(i,i) = evecsv(i,i) + evalfv(ist)
    enddo
  enddo
  write(*,*) '!!! sum evecsv before diagonalization = ', sum(evecsv)
  !
  if( spcpl .or. (.not. spinpol) ) then
    ! spins are coupled; or spin-unpolarised: full diagonalisation
    call eveqnz(nstsv, nstsv, evecsv, evalsvp)
  else
    ! spins not coupled: block diagonalise H
    call eveqnz(nstfv, nstsv, evecsv, evalsvp)
    i = nstfv + 1
    call eveqnz(nstfv, nstsv, evecsv(i,i), evalsvp(i))
    do i = 1,nstfv
      do j = 1,nstfv
        evecsv(i, j+nstfv) = 0.d0
        evecsv(i+nstfv, j) = 0.d0
      enddo
    enddo
  endif

  write(*,*)
  write(*,*) '</div> EXIT my_eveqnsv'
  write(*,*)

  return

end subroutine

