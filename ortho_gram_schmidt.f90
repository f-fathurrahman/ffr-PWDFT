! real(8) version
!--------------------------------------------------
SUBROUTINE z_ortho_gram_schmidt( v, ldv, nrow, ncols )
!--------------------------------------------------
  IMPLICIT NONE
  INTEGER :: ldv, nrow, ncols
  COMPLEX(8) :: v(ldv,ncols)
  !
  INTEGER :: ii, jj
  COMPLEX(8) :: zz, puv
  !
  COMPLEX(8) :: zdotc

  DO ii = 1, ncols
    zz = zdotc( nrow, v(1:nrow,ii),1, v(1:nrow,ii),1 )
    v(1:nrow,ii) = v(1:nrow,ii)/sqrt( real(zz,kind=8) )
    DO jj = ii+1, ncols
      puv = prj( nrow, v(1:nrow,ii), v(1:nrow,jj) )
      v(1:nrow,jj) = v(1:nrow,jj) - puv*v(1:nrow,ii)
    ENDDO
  ENDDO

  CONTAINS

    ! compute prj = <v|u>/<u|u>
    FUNCTION prj(N,u,v)
      IMPLICIT NONE
      !
      COMPLEX(8) :: prj
      INTEGER :: N
      COMPLEX(8) :: u(N), v(N)
      !
      COMPLEX(8) :: vu, uu
      COMPLEX(8) :: zdotc
      !
      ! FIXME: I got the vectors to be orthogonal when I reverse the arguments
      ! for zdotc
      vu = zdotc( N, u,1, v,1 )
      uu = zdotc( N, u,1, u,1 )
      prj = vu/uu
    END FUNCTION

END SUBROUTINE

