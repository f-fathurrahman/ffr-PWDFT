!INPUT/OUTPUT PARAMETERS:
!
!   ngp    : number of G+p-vectors (in,integer)
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   gpc    : length of G+p-vectors (in,real(ngkmax))
!   sfacgp : structure factors of G+p-vectors (in,complex(ngkmax,natmtot))
!   apwalm : APW matching coefficients
!            (out,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
SUBROUTINE match(ngp,vgpc,gpc,sfacgp,apwalm)
  USE m_constants, ONLY: zil, fourpi
  USE m_gkvectors, ONLY: ngkmax
  USE m_atoms, ONLY: natoms, natmtot, idxas, nspecies, rsp
  USE m_muffin_tins, ONLY: lmmaxapw, idxlm, rmt, nrmt, lmaxapw
  USE m_lattice, ONLY: omega
  USE m_apwlo, ONLY: apwordmax, apword, apwfr, npapw
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp
  REAL(8), intent(in) :: vgpc(3,ngkmax),gpc(ngkmax)
  COMPLEX(8), intent(in) :: sfacgp(ngkmax,natmtot)
  COMPLEX(8), intent(out) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  ! local variables
  INTEGER :: is,ia,ias,omax
  INTEGER :: l,m,lm,io,jo,i
  INTEGER :: nr,ir,igp,info
  REAL(8) :: t0,t1
  COMPLEX(8) :: z1,z2,z3
  ! automatic arrays
  INTEGER :: ipiv(apwordmax)
  COMPLEX(8) :: a(apwordmax,apwordmax)
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: djl(:,:,:)
  COMPLEX(8), ALLOCATABLE :: ylmgp(:,:),b(:,:)
  ! external functions
  REAL(8) polynm
  external polynm

  ALLOCATE(djl(0:lmaxapw,apwordmax,ngp),ylmgp(lmmaxapw,ngp))
  IF(apwordmax > 1) ALLOCATE(b(apwordmax,ngp*(2*lmaxapw+1)))
  
  ! compute the spherical harmonics of the G+p-vectors
  DO igp=1,ngp
    CALL genylmv(lmaxapw,vgpc(:,igp),ylmgp(:,igp))
  ENDDO 
  t0=fourpi/sqrt(omega)
  
  ! loop over species
  DO is=1,nspecies
    nr=nrmt(is)
    ! maximum APW order for this species
    omax=maxval(apword(1:lmaxapw,is))
    ! special case of omax=1
    IF(omax==1) THEN 
      DO igp=1,ngp
        t1=gpc(igp)*rmt(is)
        CALL sbessel(lmaxapw,t1,djl(:,1,igp))
      ENDDO 
      DO ia=1,natoms(is)
        ias=idxas(ia,is)
        DO l=0,lmaxapw
          z1=(t0/apwfr(nr,1,1,l,ias))*zil(l)
          DO igp=1,ngp
            z2=djl(l,1,igp)*z1*sfacgp(igp,ias)
            DO m=-l,l
              lm=idxlm(l,m)
              apwalm(igp,1,lm,ias)=z2*conjg(ylmgp(lm,igp))
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
      CYCLE 
    ENDIF 
    ! starting point on radial mesh for fitting polynomial of order npapw
    ir=nr-npapw+1
    ! evaluate the spherical Bessel function derivatives for all G+p-vectors
    DO igp=1,ngp
      t1=gpc(igp)*rmt(is)
      DO io=1,omax
        CALL sbesseldm(io-1,lmaxapw,t1,djl(:,io,igp))
      ENDDO 
      t1=1.d0
      DO io=2,omax
        t1=t1*gpc(igp)
        djl(:,io,igp)=t1*djl(:,io,igp)
      ENDDO 
    ENDDO 
    ! loop over atoms
    DO ia=1,natoms(is)
      ias=idxas(ia,is)
      ! begin loop over l
      DO l=0,lmaxapw
        z1=t0*zil(l)
        ! set up matrix of derivatives
        DO jo=1,apword(l,is)
          DO io=1,apword(l,is)
            a(io,jo)=polynm(io-1,npapw,rsp(ir,is),apwfr(ir,1,jo,l,ias),rmt(is))
          ENDDO 
        ENDDO 
        ! set up target vectors
        i=0
        DO igp=1,ngp
          z2=z1*sfacgp(igp,ias)
          DO m=-l,l
            lm=idxlm(l,m)
            i=i+1
            z3=z2*conjg(ylmgp(lm,igp))
            DO io=1,apword(l,is)
              b(io,i)=djl(l,io,igp)*z3
            ENDDO 
          ENDDO 
        ENDDO 
        ! solve the general complex linear systems
        CALL zgesv(apword(l,is),i,a,apwordmax,ipiv,b,apwordmax,info)
        IF(info.ne.0) THEN 
          WRITE(*,*)
          WRITE(*,'("Error(match): could not find APW matching coefficients")')
          WRITE(*,'(" for species ",I4," and atom ",I4)') is,ia
          WRITE(*,'(" ZGESV RETURN ed INFO = ",I8)') info
          WRITE(*,*)
          stop
        ENDIF 
        i=0
        DO igp=1,ngp
          DO m=-l,l
            lm=idxlm(l,m)
            i=i+1
            DO io=1,apword(l,is)
              apwalm(igp,io,lm,ias)=b(io,i)
            ENDDO 
          ENDDO 
        ENDDO 
        ! end loop over l
      ENDDO 
      ! end loops over atoms and species
    ENDDO 
  ENDDO 
  DEALLOCATE(djl,ylmgp)

  IF(apwordmax > 1) DEALLOCATE(b)
  
  RETURN 
END SUBROUTINE 
