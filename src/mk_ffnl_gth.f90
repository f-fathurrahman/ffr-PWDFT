!
! Adapted from gth.f90 QE-6.1
!
! Copyright (C) 2015 Sebastiano Caravati
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .!
!
! modified by Fadjar Fathurrahman

!-----------------------------------------------------------------------
SUBROUTINE mk_ffnl_gth(isp, ibeta, nq, qg, vq)
!-----------------------------------------------------------------------
  !
  USE m_constants, ONLY: PI, FPI => FOURPI
  USE m_cell, ONLY: omega => CellVolume
  USE m_PsPot, ONLY : Ps_HGH_Params, PsPot_ipr, PsPot_lll

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  !
  ! I/O
  INTEGER,  INTENT(in)  :: isp, ibeta, nq
  REAL(dp), INTENT(in)  :: qg(nq)
  REAL(dp), INTENT(out) :: vq(nq)
  !
  ! Local variables
  INTEGER, PARAMETER :: nprj_max(0:3)=[3, 3, 2, 1]
  INTEGER  :: ii, ll, iproj
  REAL(dp) :: rrl, qr2, fact
  !
  
  iproj = PsPot_ipr(isp,ibeta)
  ll    = PsPot_lll(isp,ibeta)
  rrl   = Ps_HGH_Params(isp)%rc(ll)

  IF( ll < 0 .OR. ll > 3 ) THEN 
    WRITE(*,*) 'ERROR in mk_ffnl_gth: wrong l, l = ', ll
    STOP 
  ENDIF 
  
  IF( iproj > nprj_max(ll) ) THEN 
    WRITE(*,*) 'ERROR in mk_ffnl_gth: projector exceeds max. n. of projectors, iproj = ', iproj
    STOP 
  ENDIF 
  !
  lif: if (ll==0) then     ! s channel
    !
    if(iproj==1)then
       do ii=1,nq
          qr2=(qg(ii)*rrl)**2
          vq(ii)=exp(-0.5_dp*qr2)
       end do
    else if(iproj==2)then
       do ii=1,nq
          qr2=(qg(ii)*rrl)**2
          vq(ii)=2._dp/sqrt(15._dp) * exp(-0.5_dp*qr2) * ( 3._dp-qr2 )
       end do
    else if(iproj==3)then
       do ii=1,nq
          qr2=(qg(ii)*rrl)**2
          vq(ii)=(4._dp/3._dp)/sqrt(105._dp) * exp(-0.5_dp*qr2) * &
               &         (15._dp-10._dp*qr2 + qr2**2)
       end do
    end if
    !
  else if (ll==1) then lif ! p channel
     !
     if(iproj==1)then
        do ii=1,nq
           qr2=(qg(ii)*rrl)**2
           vq(ii)=(1._dp/sqrt(3._dp)) * exp(-0.5_dp*qr2) * qg(ii)
        end do
     else if(iproj==2)then
        do ii=1,nq
           qr2=(qg(ii)*rrl)**2
           vq(ii)=(2._dp/sqrt(105._dp)) * exp(-0.5_dp*qr2) * qg(ii)*(5._dp-qr2)
        end do
     else if(iproj==3)then
        do ii=1,nq
           qr2=(qg(ii)*rrl)**2
           vq(ii)=(4._dp/3._dp)/sqrt(1155._dp) * exp(-0.5_dp*qr2) * &
                &         qg(ii) * (35._dp-14._dp*qr2+qr2**2)
        end do
     end if
     !
  else if (ll==2) then lif ! d channel [ ONLY 2 PROJECTORS!! ]
     !
     if(iproj==1)then
        do ii=1,nq
           qr2=(qg(ii)*rrl)**2
           vq(ii)=(1._dp/sqrt(15._dp)) * exp(-0.5_dp*qr2) * qg(ii)**2
        end do
     else if(iproj==2)then
        do ii=1,nq
           qr2=(qg(ii)*rrl)**2
           vq(ii)=(2._dp/3._dp)/sqrt(105._dp) * exp(-0.5_dp*qr2) * &
                &         qg(ii)**2 * (7._dp-qr2)
        end do
     end if
     !
  else if (ll==3) then lif ! f channel [ ONLY 1 PROJECTOR!! ]
     !
     do ii=1,nq
        qr2=(qg(ii)*rrl)**2
        vq(ii)=qg(ii)**3 * exp(-0.5_dp*qr2)
     end do
     !
  end if lif
  !
  fact = FPI * PI**0.25_dp * sqrt( 2._dp**(ll+1) * rrl**(2*ll+3) / omega )
  vq(:) = fact*vq(:)
  !
END SUBROUTINE 


