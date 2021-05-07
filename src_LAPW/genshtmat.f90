SUBROUTINE genshtmat()
  USE m_muffin_tins, ONLY: lmmaxo, lmmaxi, lmaxo, lmaxi
  USE m_sht, ONLY: zfshto, zfshti, trotsht, rotsht, rfshto, &
                   rfshti, zbshto, zbshti, rbshto, rbshti
  IMPLICIT NONE 
  ! local variables
  INTEGER itp
  REAL(8) v(3)
  ! automatic arrays
  REAL(8) tp(2,lmmaxo),vtp(3,lmmaxo),rlm(lmmaxo)
  COMPLEX(8) ylm(lmmaxo)
  
  !--------------------------------!
  !     SHT matrices for lmaxo     !
  !--------------------------------!
  ! allocate real SHT matrices
  IF(allocated(rbshto)) DEALLOCATE(rbshto)
  ALLOCATE(rbshto(lmmaxo,lmmaxo))
  IF(allocated(rfshto)) DEALLOCATE(rfshto)
  ALLOCATE(rfshto(lmmaxo,lmmaxo))
  ! allocate complex SHT matrices
  IF(allocated(zbshto)) DEALLOCATE(zbshto)
  ALLOCATE(zbshto(lmmaxo,lmmaxo))
  IF(allocated(zfshto)) DEALLOCATE(zfshto)
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
  rfshto(:,:)=rbshto(:,:)
  CALL rminv(lmmaxo,rfshto)
  ! complex
  zfshto(:,:)=zbshto(:,:)
  CALL zminv(lmmaxo,zfshto)

  !--------------------------------!
  !     SHT matrices for lmaxi     !
  !--------------------------------!
  ! allocate real SHT matrices
  IF(allocated(rbshti)) DEALLOCATE(rbshti)
  ALLOCATE(rbshti(lmmaxi,lmmaxi))
  IF(allocated(rfshti)) DEALLOCATE(rfshti)
  ALLOCATE(rfshti(lmmaxi,lmmaxi))
  ! allocate complex SHT matrices
  IF(allocated(zbshti)) DEALLOCATE(zbshti)
  ALLOCATE(zbshti(lmmaxi,lmmaxi))
  IF(allocated(zfshti)) DEALLOCATE(zfshti)
  ALLOCATE(zfshti(lmmaxi,lmmaxi))
  ! generate spherical covering set for lmaxi
  CALL sphcover(lmmaxi,tp)
  ! convert (theta, phi) angles to vectors
  DO itp=1,lmmaxi
    CALL sctovec(tp(:,itp),vtp(:,itp))
  ENDDO 
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
  rfshti(:,:)=rbshti(:,:)
  CALL rminv(lmmaxi,rfshti)
  ! complex
  zfshti(:,:)=zbshti(:,:)
  CALL zminv(lmmaxi,zfshti)
  RETURN 
END SUBROUTINE 
