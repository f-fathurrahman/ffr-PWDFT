subroutine oepvcl(vclcv,vclvv)
  use modmain
  implicit none
  ! arguments
  complex(8), intent(out) :: vclcv(ncrmax,natmtot,nstsv,nkpt)
  complex(8), intent(out) :: vclvv(nstsv,nstsv,nkpt)
  ! local variables
  integer ik,ncv,nvv
  integer lp
  do ik=1,nkpt
    write(*,'("Info(oepvcl): ",I6," of ",I6," k-points")') ik,nkpt
    call oepvclk(ik,vclcv(:,:,:,ik),vclvv(:,:,ik))
  end do
  return
end subroutine

