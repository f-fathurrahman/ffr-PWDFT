program main
  IMPLICIT NONE 
  INTEGER :: lmaxi, lmaxo

  lmaxi = 1
  lmaxo = 6
  CALL debug_genshtmat(lmaxi, lmaxo)
end



SUBROUTINE debug_genshtmat(lmaxi, lmaxo)
  USE m_sht, ONLY: trotsht
  IMPLICIT NONE
  !
  INTEGER :: lmaxi, lmaxo
  INTEGER :: lmmaxi, lmmaxo 
  ! variables from m_sht
  REAL(8), ALLOCATABLE :: rbshti(:,:)
  REAL(8), ALLOCATABLE :: rfshti(:,:)
  REAL(8), ALLOCATABLE :: rbshto(:,:)
  REAL(8), ALLOCATABLE :: rfshto(:,:)
  !
  COMPLEX(8), ALLOCATABLE :: zbshti(:,:)
  COMPLEX(8), ALLOCATABLE :: zfshti(:,:)
  COMPLEX(8), ALLOCATABLE :: zbshto(:,:)
  COMPLEX(8), ALLOCATABLE :: zfshto(:,:)
  !
  REAL(8) :: rotsht(3,3)
  ! local variables
  INTEGER :: itp, i, j
  REAL(8) :: v(3)
  ! automatic arrays
  REAL(8), ALLOCATABLE :: tp(:,:),vtp(:,:),rlm(:)
  COMPLEX(8), ALLOCATABLE :: ylm(:)
  real(8), allocatable :: tmpmat(:,:)

  lmmaxi = (lmaxi + 1)**2
  lmmaxo = (lmaxo + 1)**2

  WRITE(*,*) 'ENTER debug_genshtmat'

  ALLOCATE( tp(2,lmmaxo) )
  ALLOCATE( vtp(3,lmmaxo) )
  ALLOCATE( rlm(lmmaxo) )
  ALLOCATE( ylm(lmmaxo) )

  !
  ! SHT matrices for lmaxo
  !
  ! allocate real SHT matrices
  ALLOCATE(rbshto(lmmaxo,lmmaxo))
  ALLOCATE(rfshto(lmmaxo,lmmaxo))
  ! allocate complex SHT matrices
  ALLOCATE(zbshto(lmmaxo,lmmaxo))
  ALLOCATE(zfshto(lmmaxo,lmmaxo))
  
  ! generate spherical covering set
  CALL sphcover(lmmaxo,tp)
  ! convert (theta, phi) angles to vectors
  DO itp=1,lmmaxo
    CALL sctovec(tp(:,itp),vtp(:,itp))
  ENDDO 

  ! rotate the spherical covering set if required
  IF(trotsht) THEN 
    DO itp=1,lmmaxo
      v(:)=vtp(:,itp)
      CALL r3mv(rotsht,v,vtp(:,itp))
    ENDDO 
  ENDIF 
  ! generate real and complex spherical harmonics and set the backward SHT arrays
  DO itp=1,lmmaxo
    CALL genrlmv(lmaxo,vtp(:,itp),rlm)
    rbshto(itp,1:lmmaxo)=rlm(1:lmmaxo)
    CALL genylmv(lmaxo,vtp(:,itp),ylm)
    zbshto(itp,1:lmmaxo)=ylm(1:lmmaxo)
  ENDDO 
  ! find the forward SHT arrays
  ! real
  rfshto(:,:)= rbshto(:,:)
  CALL rminv(lmmaxo,rfshto)
  ! complex
  zfshto(:,:) = zbshto(:,:)
  CALL zminv(lmmaxo,zfshto)

  write(*,*) 'lmmaxo = ', lmmaxo
  do i=1,lmmaxo
    do j=1,lmmaxo
      write(100,'(1x,F10.5)',advance="no") rbshto(i,j)
    enddo 
    write(100,*)
  enddo

  !--------------------------------!
  !     SHT matrices for lmaxi     !
  !--------------------------------!
  ! allocate real SHT matrices
  ALLOCATE(rbshti(lmmaxi,lmmaxi))
  ALLOCATE(rfshti(lmmaxi,lmmaxi))
  ! allocate complex SHT matrices
  ALLOCATE(zbshti(lmmaxi,lmmaxi))
  ALLOCATE(zfshti(lmmaxi,lmmaxi))

  ! generate spherical covering set for lmaxi
  CALL sphcover(lmmaxi,tp)
  !
  write(*,*) 'theta and phi for lmmaxi'
  do i = 1,lmmaxi
    write(*,'(1x,I4,2F18.10)') i, tp(:,i)
  enddo
  !
  ! convert (theta, phi) angles to vectors
  DO itp=1,lmmaxi
    CALL sctovec(tp(:,itp),vtp(:,itp))
  ENDDO 
  !
  write(*,*) 'After sctovec for lmmaxi'
  do i = 1,lmmaxi
    write(*,'(1x,I4,3F18.10)') i, vtp(:,i)
  enddo
  ! rotate the spherical covering set if required
  IF(trotsht) THEN 
    DO itp=1,lmmaxi
      v(:)=vtp(:,itp)
      CALL r3mv(rotsht,v,vtp(:,itp))
    ENDDO 
  ENDIF 
  ! generate real and complex spherical harmonics and set the backward SHT arrays
  DO itp=1,lmmaxi
    CALL genrlmv(lmaxi,vtp(:,itp),rlm)
    rbshti(itp,1:lmmaxi)=rlm(1:lmmaxi)
    CALL genylmv(lmaxi,vtp(:,itp),ylm)
    zbshti(itp,1:lmmaxi)=ylm(1:lmmaxi)
  ENDDO 
  ! find the forward SHT arrays
  ! real
  rfshti(:,:) = rbshti(:,:)
  CALL rminv(lmmaxi,rfshti)
  ! complex
  zfshti(:,:) = zbshti(:,:)
  CALL zminv(lmmaxi,zfshti)

  write(*,*) 'lmmaxi = ', lmmaxi
  do i=1,lmmaxi
    do j=1,lmmaxi
      write(101,'(1x,F10.5)',advance="no") rbshti(i,j)
    enddo 
    write(101,*)
  enddo

  allocate(tmpmat(lmmaxi,lmmaxi))
  tmpmat = matmul(rfshti,rbshti)
  !tmpmat = matmul(rbshti,rfshti)

  write(*,*) 'tmpmat'
  do i=1,lmmaxi
    do j=1,lmmaxi
      write(*,'(1x,F10.5)',advance="no") tmpmat(i,j)
    enddo 
    write(*,*)
  enddo

  deallocate(tmpmat)
  DEALLOCATE(tp, vtp, rlm, ylm)

  WRITE(*,*) 'EXIT  debug_genshtmat'

  RETURN 
END SUBROUTINE 
