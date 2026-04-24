subroutine my_gradzfmt(nr,nri,r,ri,zfmt,ld,gzfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr),ri(nr)
complex(8), intent(in) :: zfmt(*)
integer, intent(in) :: ld
complex(8), intent(out) :: gzfmt(ld,3)
! local variables
integer nro,iro,ir,i,j
integer np,npi,i1,i2
integer l,m,lm1,lm2
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
real(8) t1,t2,t3
complex(8) z1
! automatic arrays
real(8) f1(nr),f2(nr),g1(nr),g2(nr)
complex(8) drmt(ld)
! external functions
real(8), external :: clebgor

nro = nr-nri
iro = nri + 1
npi = lmmaxi*nri
np = npi + lmmaxo*nro

!----------------------------------------!
!     compute the radial derivatives     !
!----------------------------------------!
do lm1 = 1,lmmaxi
  i1 = lm1
  do ir = 1,nri
    f1(ir) = dble(zfmt(i1))
    f2(ir) = aimag(zfmt(i1))
    i1 = i1 + lmmaxi
  enddo
  do ir = iro,nr
    f1(ir) = dble(zfmt(i1))
    f2(ir) = aimag(zfmt(i1))
    i1 = i1 + lmmaxo
  enddo
  call fderiv(1, nr, r, f1, g1)
  call fderiv(1, nr, r, f2, g2)
  i1 = lm1
  do ir = 1,nri
    drmt(i1) = cmplx(g1(ir),g2(ir),8)
    i1 = i1 + lmmaxi
  enddo
  do ir = iro,nr
    drmt(i1) = cmplx(g1(ir),g2(ir),8)
    i1 = i1 + lmmaxo
  enddo
enddo

do lm1 = lmmaxi+1,lmmaxo
  i1 = npi+lm1
  do ir = iro,nr
    f1(ir) = dble(zfmt(i1))
    f2(ir) = aimag(zfmt(i1))
    i1 = i1 + lmmaxo
  enddo
  call fderiv(1, nro, r(iro), f1(iro), g1(iro))
  call fderiv(1, nro, r(iro), f2(iro), g2(iro))
  i1 = npi + lm1
  do ir = iro,nr
    drmt(i1) = cmplx(g1(ir),g2(ir),8)
    i1 = i1 + lmmaxo
  enddo
enddo

!-------------------------------------------------------!
!     compute the gradient in spherical coordinates     !
!-------------------------------------------------------!
! zero the gradient array
gzfmt(1:np,:) = 0.d0
! inner part of muffin-tin
lm1 = 0
do l = 0,lmaxi
  t1 = sqrt(dble(l+1)/dble(2*l+3))
  if (l > 0) then
    t2 = -sqrt(dble(l)/dble(2*l-1))
  else
    t2 = 0.d0
  endif
  do m = -l,l
    lm1 = lm1 + 1
    j = 1
    do i = -1,1
      if( i == 0 ) then
        j = 3
      endif
      if( i == 1 ) then
        j = 2
      endif
      if( (l+1 <= lmaxi ) .and. ( abs(m+i) <= l+1) ) then
        ! index to (l,m) is l*(l+1)+m+1, therefore index to (l+1,m+i) is
        lm2 = (l+1)*(l+2) + (m+i)+1
        t3 = t1*clebgor(l, 1, l+1, m, i, m+i)
        i1 = lm1
        i2 = lm2
        do ir = 1,nri
          gzfmt(i2,j) = gzfmt(i2,j) + t3*(drmt(i1) - dble(l)*ri(ir)*zfmt(i1))
          i1 = i1 + lmmaxi
          i2 = i2 + lmmaxi
        enddo
      endif
      if (abs(m+i) <= l-1) then
        ! index to (l-1,m+i)
        lm2 = (l-1)*l + (m+i) + 1
        t3 = t2*clebgor(l,1,l-1,m,i,m+i)
        i1 = lm1
        i2 = lm2
        do ir=1,nri
          gzfmt(i2,j) = gzfmt(i2,j) + t3*(drmt(i1) + dble(l+1)*ri(ir)*zfmt(i1))
          i1 = i1 + lmmaxi
          i2 = i2 + lmmaxi
        enddo
      endif
    enddo
  enddo
enddo
! outer part of muffin-tin
lm1 = 0
do l = 0,lmaxo
  t1 = sqrt(dble(l+1)/dble(2*l+3))
  if (l > 0) then
    t2 = -sqrt(dble(l)/dble(2*l-1))
  else
    t2 = 0.d0
  endif
  do m = -l,l
    lm1 = lm1 + 1
    j = 1
    do i = -1,1
      if (i == 0) then
        j = 3
      endif
      if (i == 1) then
        j = 2
      endif
      if( (l+1 <= lmaxo) .and. (abs(m+i) <= l+1) ) then
        lm2 = (l+1)*(l+2) + (m+i) + 1
        t3 = t1*clebgor(l, 1, l+1, m, i, m+i)
        i1 = npi + lm1
        i2 = npi + lm2
        do ir = iro,nr
          gzfmt(i2,j) = gzfmt(i2,j) + t3*(drmt(i1) - dble(l)*ri(ir)*zfmt(i1))
          i1 = i1 + lmmaxo
          i2 = i2 + lmmaxo
        enddo
      endif
      if( abs(m+i) <= l-1 ) then
        lm2 = (l-1)*l+(m+i)+1
        t3 = t2*clebgor(l,1,l-1,m,i,m+i)
        i1 = npi+lm1
        i2 = npi+lm2
        do ir = iro,nr
          gzfmt(i2,j) = gzfmt(i2,j) + t3*(drmt(i1)+dble(l+1)*ri(ir)*zfmt(i1))
          i1 = i1 + lmmaxo
          i2 = i2 + lmmaxo
        enddo
      endif
    enddo
  enddo
enddo

!--------------------------------------------------------!
!     convert from spherical components to Cartesian     !
!--------------------------------------------------------!
i1 = 0
do ir = 1,nri
  do lm1 = 1,lmmaxi
    i1 = i1 + 1
    z1 = gzfmt(i1,1)
    gzfmt(i1,1) = c1*(z1-gzfmt(i1,2))
    z1 = c1*(z1+gzfmt(i1,2))
    gzfmt(i1,2) = cmplx(-aimag(z1),dble(z1),8)
  enddo
enddo

do ir = iro,nr
  do lm1 = 1,lmmaxo
    i1 = i1 + 1
    z1 = gzfmt(i1,1)
    gzfmt(i1,1) = c1*(z1-gzfmt(i1,2))
    z1 = c1*(z1 + gzfmt(i1,2))
    gzfmt(i1,2) = cmplx(-aimag(z1),dble(z1),8)
  enddo
enddo

return
end subroutine
