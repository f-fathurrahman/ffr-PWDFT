SUBROUTINE symrfir(rfir)
  USE m_gvectors, ONLY: ngtot, ngridg, ngvec, vgc, ivgig, ivg, igfft
  USE m_symmetry, ONLY: symlat, nsymcrys, tv0symc, isymlat, vtcsymc, lsplsymc
! !INPUT/OUTPUT PARAMETERS:
!   rfir : real intersitial function (inout,real(ngtot))
! !DESCRIPTION:
!   Symmetrises a real scalar interstitial function. The function is first
!   Fourier transformed to $G$-space, and THEN  averaged over each symmetry by
!   rotating the Fourier coefficients and multiplying them by a phase factor
!   corresponding to the symmetry translation.
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(inout) :: rfir(ngtot)
  ! local variables
  LOGICAL tv0
  INTEGER isym,lspl,ilspl
  INTEGER sym(3,3),ig,ifg,jfg
  INTEGER i1,i2,i3,j1,j2,j3
  REAL(8) v1,v2,v3,t1
  COMPLEX(8) z1
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: zfft1(:),zfft2(:)
  ALLOCATE(zfft1(ngtot),zfft2(ngtot))
  ! Fourier transform function to G-space
  zfft1(:)=rfir(:)
  CALL zfftifc(3,ngridg,-1,zfft1)
  zfft2(:)=0.d0
  ! loop over crystal symmetries
  DO isym=1,nsymcrys
  ! zero translation vector flag
    tv0=tv0symc(isym)
  ! translation vector in Cartesian coordinates
    IF(.not.tv0) THEN 
      v1=vtcsymc(1,isym)
      v2=vtcsymc(2,isym)
      v3=vtcsymc(3,isym)
    ENDIF 
  ! index to lattice symmetry of spatial rotation
    lspl=lsplsymc(isym)
  ! inverse rotation required for rotation of G-vectors
    ilspl=isymlat(lspl)
    sym(:,:)=symlat(:,:,ilspl)
    DO ig=1,ngvec
      ifg=igfft(ig)
  ! multiply the transpose of the inverse symmetry matrix with the G-vector
      IF(lspl.eq.1) THEN 
        jfg=ifg
      else
        i1=ivg(1,ig); i2=ivg(2,ig); i3=ivg(3,ig)
        j1=sym(1,1)*i1+sym(2,1)*i2+sym(3,1)*i3
        j2=sym(1,2)*i1+sym(2,2)*i2+sym(3,2)*i3
        j3=sym(1,3)*i1+sym(2,3)*i2+sym(3,3)*i3
        jfg=igfft(ivgig(j1,j2,j3))
      ENDIF 
      IF(tv0) THEN 
  ! zero translation vector
        zfft2(jfg)=zfft2(jfg)+zfft1(ifg)
      else
  ! complex phase factor for translation
        t1=-(vgc(1,ig)*v1+vgc(2,ig)*v2+vgc(3,ig)*v3)
        z1=cmplx(cos(t1),sin(t1),8)
        zfft2(jfg)=zfft2(jfg)+z1*zfft1(ifg)
      ENDIF 
    ENDDO 
  ENDDO 
  ! Fourier transform to real-space and normalise
  CALL zfftifc(3,ngridg,1,zfft2)
  t1=1.d0/dble(nsymcrys)
  rfir(:)=t1*dble(zfft2(:))
  DEALLOCATE(zfft1,zfft2)
  RETURN 
END SUBROUTINE 