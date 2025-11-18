subroutine my_symrvfir(tspin,tnc,rvfir)
  use modmain
  implicit none
  ! arguments
  logical, intent(in) :: tspin,tnc
  real(8), intent(inout) :: rvfir(ngtot,*)
  ! local variables
  logical tv0
  integer nd,isym,lspl,ilspl,lspn
  integer sym(3,3),ig,ifg,jfg
  integer i1,i2,i3,j1,j2,j3,i
  real(8) sc(3,3),v1,v2,v3,t1
  complex(8) zv1(3),zv2(3),z1
  ! allocatable arrays
  complex(8), allocatable :: zfft1(:,:),zfft2(:,:)
  
  ! dimension of the vector field
  if (tnc) then
    nd=3
  else
    nd=1
  end if
  
  allocate(zfft1(ngtot,nd),zfft2(ngtot,nd))
  ! Fourier transform vector function to G-space
  do i=1,nd
    zfft1(:,i)=rvfir(:,i)
    call zfftifc(3,ngridg,-1,zfft1(:,i))
  end do
  zfft2(:,:)=0.d0
  do isym=1,nsymcrys
    ! zero translation vector flag
    tv0 = tv0symc(isym)
    ! translation vector in Cartesian coordinates
    if (.not.tv0) then
      v1=vtcsymc(1,isym)
      v2=vtcsymc(2,isym)
      v3=vtcsymc(3,isym)
    end if
    ! index to spatial rotation lattice symmetry
    lspl=lsplsymc(isym)
    ! inverse rotation required for rotation of G-vectors
    ilspl=isymlat(lspl)
    sym(:,:)=symlat(:,:,ilspl)
    if (tspin) then
      ! global spin proper rotation in Cartesian coordinates
      lspn=lspnsymc(isym)
      sc(:,:)=symlatd(lspn)*symlatc(:,:,lspn)
    else
      ! set spin rotation equal to spatial rotation
      lspn=lspl
      sc(:,:)=symlatc(:,:,lspl)
    end if
    do ig=1,ngvec
      ifg=igfft(ig)
      ! multiply the transpose of the inverse symmetry matrix with the G-vector
      if (lspl.eq.1) then
        jfg=ifg
      else
        i1=ivg(1,ig); i2=ivg(2,ig); i3=ivg(3,ig)
        j1=sym(1,1)*i1+sym(2,1)*i2+sym(3,1)*i3
        j2=sym(1,2)*i1+sym(2,2)*i2+sym(3,2)*i3
        j3=sym(1,3)*i1+sym(2,3)*i2+sym(3,3)*i3
        jfg=igfft(ivgig(j1,j2,j3))
      end if
      ! complex phase factor for translation
      if (tv0) then
        z1=1.d0
      else
        t1=-(vgc(1,ig)*v1+vgc(2,ig)*v2+vgc(3,ig)*v3)
        z1=cmplx(cos(t1),sin(t1),8)
      end if
      ! translation, spatial rotation and global spin rotation
      if (lspn.eq.1) then
        ! global spin symmetry is the identity
        zfft2(jfg,:)=zfft2(jfg,:)+z1*zfft1(ifg,:)
      else
        if (tnc) then
          ! non-collinear case
          zv1(:)=zfft1(ifg,:)
          zv2(1)=sc(1,1)*zv1(1)+sc(1,2)*zv1(2)+sc(1,3)*zv1(3)
          zv2(2)=sc(2,1)*zv1(1)+sc(2,2)*zv1(2)+sc(2,3)*zv1(3)
          zv2(3)=sc(3,1)*zv1(1)+sc(3,2)*zv1(2)+sc(3,3)*zv1(3)
          zfft2(jfg,:)=zfft2(jfg,:)+z1*zv2(:)
        else
          ! collinear case
          zfft2(jfg,1)=zfft2(jfg,1)+sc(3,3)*z1*zfft1(ifg,1)
        end if
      end if
    end do
  end do
  ! Fourier transform to real-space and normalise
  t1=1.d0/dble(nsymcrys)
  do i=1,nd
    call zfftifc(3,ngridg,1,zfft2(:,i))
    rvfir(:,i)=t1*dble(zfft2(:,i))
  end do
  deallocate(zfft1,zfft2)
  return
end subroutine

