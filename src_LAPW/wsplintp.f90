subroutine wsplintp(n,x,w)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: w(4,n)
! local variables
integer i
real(8) f(4),t1,t2
! external functions
real(8) polynm
external polynm
if (n.lt.4) then
  write(*,*)
  write(*,'("Error(wsplintp): n < 4 : ",I8)') n
  write(*,*)
  stop
end if
w(:,1)=0.d0
f(:)=0.d0
f(1)=1.d0
w(1,2)=polynm(-1,4,x,f,x(2))
f(1)=0.d0
f(2)=1.d0
w(2,2)=polynm(-1,4,x,f,x(2))
f(2)=0.d0
f(3)=1.d0
w(3,2)=polynm(-1,4,x,f,x(2))
f(3)=0.d0
f(4)=1.d0
w(4,2)=polynm(-1,4,x,f,x(2))
do i=3,n-1
  f(:)=0.d0
  f(1)=1.d0
  t1=polynm(-1,4,x(i-2),f,x(i-1))
  t2=polynm(-1,4,x(i-2),f,x(i))
  w(1,i)=t2-t1
  f(1)=0.d0
  f(2)=1.d0
  t1=polynm(-1,4,x(i-2),f,x(i-1))
  t2=polynm(-1,4,x(i-2),f,x(i))
  w(2,i)=t2-t1
  f(2)=0.d0
  f(3)=1.d0
  t1=polynm(-1,4,x(i-2),f,x(i-1))
  t2=polynm(-1,4,x(i-2),f,x(i))
  w(3,i)=t2-t1
  f(3)=0.d0
  f(4)=1.d0
  t1=polynm(-1,4,x(i-2),f,x(i-1))
  t2=polynm(-1,4,x(i-2),f,x(i))
  w(4,i)=t2-t1
end do
f(:)=0.d0
f(1)=1.d0
t1=polynm(-1,4,x(n-3),f,x(n-1))
t2=polynm(-1,4,x(n-3),f,x(n))
w(1,n)=t2-t1
f(1)=0.d0
f(2)=1.d0
t1=polynm(-1,4,x(n-3),f,x(n-1))
t2=polynm(-1,4,x(n-3),f,x(n))
w(2,n)=t2-t1
f(2)=0.d0
f(3)=1.d0
t1=polynm(-1,4,x(n-3),f,x(n-1))
t2=polynm(-1,4,x(n-3),f,x(n))
w(3,n)=t2-t1
f(3)=0.d0
f(4)=1.d0
t1=polynm(-1,4,x(n-3),f,x(n-1))
t2=polynm(-1,4,x(n-3),f,x(n))
w(4,n)=t2-t1
return
end subroutine

