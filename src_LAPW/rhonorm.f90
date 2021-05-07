SUBROUTINE rhonorm()
! !USES:
use modmain
! !DESCRIPTION:
!   Loss of precision of the calculated total charge can result because the
!   muffin-tin density is computed on a set of $(\theta,\phi)$ points and THEN 
!   transformed to a spherical harmonic representation. This routine adds a
!   constant to the density so that the total charge is correct. If the error in
!   total charge exceeds a certain tolerance THEN  a warning is issued.
IMPLICIT NONE 
! local variables
INTEGER is,ia,ias
INTEGER nr,nri,ir,i
REAL(8) t1,t2
IF(.not.trhonorm) RETURN 
! check error in total charge
t1=chgcalc/chgtot-1.d0
IF(abs(t1).gt.epschg) THEN 
  WRITE(*,*)
  WRITE(*,'("Warning(rhonorm): total charge density incorrect for s.c. &
   &loop ",I5)') iscl
  WRITE(*,'(" Calculated : ",G18.10)') chgcalc
  WRITE(*,'(" Required   : ",G18.10)') chgtot
ENDIF 
! error in average density
t1=(chgtot-chgcalc)/omega
! add the constant difference to the density
t2=t1/y00
DO ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  DO ir=1,nri
    rhomt(i,ias)=rhomt(i,ias)+t2
    i=i+lmmaxi
  ENDDO 
  DO ir=nri+1,nr
    rhomt(i,ias)=rhomt(i,ias)+t2
    i=i+lmmaxo
  ENDDO 
ENDDO 
rhoir(:)=rhoir(:)+t1
! add the difference to the charges
DO is=1,nspecies
  t2=t1*(4.d0*pi/3.d0)*rmt(is)**3
  DO ia=1,natoms(is)
    ias=idxas(ia,is)
    chgmt(ias)=chgmt(ias)+t2
    chgmttot=chgmttot+t2
  ENDDO 
ENDDO 
chgir=chgtot-chgmttot
RETURN 
END SUBROUTINE 
