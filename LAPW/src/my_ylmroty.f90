!----------------------------------------
subroutine my_ylmroty(beta, lmax, ld, dy)
!----------------------------------------
  implicit none
  ! arguments
  real(8), intent(in) :: beta
  integer, intent(in) :: lmax,ld
  real(8), intent(out) :: dy(ld,*)
  ! local variables
  integer j,k,l,m1,m2,lm1,lm2
  real(8) cb,sb,sum,t1,t2
  ! external functions
  real(8) factnm
  external factnm

  cb = cos(beta/2.d0)
  sb = sin(beta/2.d0)
  lm1 = 0
  do l = 0,lmax
    ! generate rotation operator for m-components of current l
    do m1 = -l,l
      lm1 = lm1 + 1
      lm2 = l**2
      do m2 = -l,l
        lm2 = lm2 + 1
        sum = 0.d0
        do k = 0,min(l+m1,l-m2)
          if( ((l+m1-k) >= 0) .and. ((l-m2-k) >= 0) .and. ((m2-m1+k) >= 0)) then
            j = 2*(l-k) + m1 - m2
            if(j == 0) then
              t1 = 1.d0
            else
              t1 = cb**j
            endif
            j = 2*k + m2 - m1
            if( j /= 0 ) then
              t1 = t1*sb**j
            endif
            t2 = t1/(factnm(k,1)*factnm(l+m1-k,1)*factnm(l-m2-k,1) *factnm(m2-m1+k,1))
            if( mod(k,2) /= 0) then
              t2 = -t2
            end
            sum = sum + t2
          endif
        enddo
        t1 = sqrt(factnm(l+m1,1)*factnm(l-m1,1)*factnm(l+m2,1)*factnm(l-m2,1))
        dy(lm1,lm2) = t1*sum
      enddo
    enddo
  enddo
  return
end subroutine


