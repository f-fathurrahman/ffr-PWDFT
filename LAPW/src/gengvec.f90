!   Generates a set of ${\bf G}$-vectors used for the Fourier transform of the
!   charge density and potential and sorts them according to length. Integers
!   corresponding to the vectors in lattice coordinates are stored, as well as
!   the map from these INTEGER coordinates to the ${\bf G}$-vector index. A map
!   from the ${\bf G}$-vector set to the standard FFT array structure is also
!   generated. Finally, the number of ${\bf G}$-vectors with magnitude less than
!   {\tt gmaxvr} is determined.
SUBROUTINE gengvec()
  USE m_gkvectors, ONLY: gkmax
  USE m_lattice, ONLY: epslat, avec, bvec
  USE m_gvectors, ONLY: ngridg, igfft, ivg, ivgig, gc, intgv, ngvec, ngtot, gmaxvr, vgc
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ig,jg,i1,i2,i3,j1,j2,j3
  REAL(8) :: v1(3),v2(3),v3(3)
  ! ALLOCATABLE arrays
  INTEGER, ALLOCATABLE :: idx(:),ivg0(:,:)
  REAL(8), ALLOCATABLE :: vgc0(:,:),gc0(:)
  
  IF(gmaxvr < 0.d0) gmaxvr = abs(gmaxvr)*gkmax

  ! ensure |G| cut-off is at least twice |G+k| cut-off
  if(gmaxvr <= 2.d0*gkmax) then
    write(*,*) 'INFO gengvec: gmaxvr will be set to 2*gkmax'
  endif
  gmaxvr=max(gmaxvr,2.d0*gkmax+epslat)
  write(*,*) 'gmaxvr = ', gmaxvr
  write(*,*) 'gkmax  = ', gkmax
  
  ! find the G-vector grid sizes
  CALL gridsize(avec,gmaxvr,ngridg,ngtot,intgv)
  
  ! allocate global G-vector arrays
  IF(allocated(ivg)) DEALLOCATE(ivg)
  ALLOCATE(ivg(3,ngtot))
  IF(allocated(ivgig)) DEALLOCATE(ivgig)
  ALLOCATE(ivgig(intgv(1,1):intgv(2,1),intgv(1,2):intgv(2,2), &
   intgv(1,3):intgv(2,3)))
  IF(allocated(igfft)) DEALLOCATE(igfft)
  ALLOCATE(igfft(ngtot))
  IF(allocated(vgc)) DEALLOCATE(vgc)
  ALLOCATE(vgc(3,ngtot))
  IF(allocated(gc)) DEALLOCATE(gc)
  ALLOCATE(gc(ngtot))

  ! allocate local arrays
  ALLOCATE(idx(ngtot),ivg0(3,ngtot))
  ALLOCATE(vgc0(3,ngtot),gc0(ngtot))
  ig=0
  DO i1=intgv(1,1),intgv(2,1)
    v1(:)=dble(i1)*bvec(:,1)
    DO i2=intgv(1,2),intgv(2,2)
      v2(:)=v1(:)+dble(i2)*bvec(:,2)
      DO i3=intgv(1,3),intgv(2,3)
        v3(:)=v2(:)+dble(i3)*bvec(:,3)
        ig = ig + 1
        ! map from G-vector to (i1,i2,i3) index
        ivg0(1,ig) = i1
        ivg0(2,ig) = i2
        ivg0(3,ig) = i3
        ! G-vector in Cartesian coordinates
        vgc0(:,ig) = v3(:)
        ! length of each G-vector
        gc0(ig) = sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)
      ENDDO 
    ENDDO 
  ENDDO 
  
  ! sort by vector length
  CALL sortidx(ngtot,gc0,idx)
  ! reorder arrays
  DO ig=1,ngtot
    jg=idx(ig)
    ivg(:,ig)=ivg0(:,jg)
    gc(ig)=gc0(jg)
    vgc(:,ig)=vgc0(:,jg)
  ENDDO 
  
  ! find the number of vectors with G < gmaxvr
  ngvec=1
  DO ig=2,ngtot
    IF(gc(ig) > gmaxvr) THEN 
      ngvec = ig-1
      EXIT 
    ENDIF 
  ENDDO 
  
  ! generate index arrays
  DO ig=1,ngtot
    i1=ivg(1,ig)
    i2=ivg(2,ig)
    i3=ivg(3,ig)
    ! map from (i1,i2,i3) to G-vector index
    ivgig(i1,i2,i3)=ig
    ! Fourier transform index
    IF(i1 >= 0) THEN 
      j1=i1
    ELSE 
      j1=ngridg(1)+i1
    ENDIF 
    IF(i2 >= 0) THEN 
      j2=i2
    ELSE 
      j2=ngridg(2)+i2
    ENDIF 
    IF(i3 >= 0) THEN 
      j3=i3
    ELSE 
      j3=ngridg(3)+i3
    ENDIF 
    igfft(ig)=j3*ngridg(2)*ngridg(1)+j2*ngridg(1)+j1+1
  ENDDO 
  DEALLOCATE(idx,ivg0,gc0,vgc0)
  RETURN 
END SUBROUTINE 
