program test_zfftifc
  implicit none
  integer :: ngridg(3), npoints
  complex(8), allocatable :: zfft(:)  

  ngridg = (/ 48, 48, 48 /)
  npoints = ngridg(1)*ngridg(2)*ngridg(3)

  allocate(zfft(npoints))

  zfft(:) = cmplx(1.d0, 2.1d0, kind=8)
  zfft(1) = cmplx(99.d0, 99.d0, kind=8)

  write(*,*) 'sum zfft before = ', sum(zfft)
  call zfftifc(3,ngridg,1,zfft)
  write(*,*) 'sum zfft after = ', sum(zfft)
  
  deallocate(zfft)

end program