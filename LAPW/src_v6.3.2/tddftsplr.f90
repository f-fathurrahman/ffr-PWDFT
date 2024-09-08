
! Copyright (C) 2013 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftsplr
use modmain
use modtest
use modmpi
use modomp
implicit none
! local variables
integer ik,isym,iq,iw
integer ig,jg,i,j,n
integer nthd
real(8) v(3)
complex(8) a(4,4),b(4,4),z1
character(256) fname
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: chi(:,:,:,:,:),chit(:),fxc(:,:,:,:)
complex(8), allocatable :: c(:,:),d(:,:,:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(tddftsplr): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
! initialise global variables
call init0
call init1
call init2
call init3
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1).gt.epslat).or.(v(2).gt.epslat).or.(v(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(tddftsplr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! find the equivalent reduced q-point
call findqpt(vecql,isym,iq)
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! generate the G+q-vectors and related quantities
allocate(vgqc(3,ngrf),gqc(ngrf),jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
call gengqrf(vecqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
deallocate(vgqc)
! initialise the OpenMP locks
allocate(lock(nwrf))
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
! compute chi0
allocate(chi(ngrf,4,ngrf,4,nwrf))
chi(:,:,:,:,:)=0.d0
call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(tddftsplr_)
  write(*,'("Info(tddftsplr): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(tddftsplr_)
  call genspchi0(ik,lock,scissor,vecql,jlgqr,ylmgq,sfacgq,chi)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(lock)
! add chi0 from each process and redistribute
if (np_mpi.gt.1) then
  n=ngrf*4*ngrf*4*nwrf
  call mpi_allreduce(mpi_in_place,chi,n,mpi_double_complex,mpi_sum,mpicom, &
   ierror)
end if
! transform chi0 from 2x2 to 1x3 basis
do iw=1,nwrf
  do ig=1,ngrf
    do jg=1,ngrf
      a(:,:)=chi(ig,:,jg,:,iw)
      call tfm2213(a,b)
      chi(ig,:,jg,:,iw)=b(:,:)
    end do
  end do
end do
! generate transverse chi0 for the collinear case
if (.not.ncmag) then
  allocate(chit(nwrf))
  do iw=1,nwrf
    a(:,:)=chi(1,:,1,:,iw)
    call tfm13t(a,b)
    chit(iw)=b(2,2)
  end do
end if
! write chi0 to file
if (mp_mpi) then
! write chi0 to file in 1x3 basis
  do i=1,4
    do j=1,4
      write(fname,'("CHI0_",2I1,".OUT")') i-1,j-1
      open(50,file=trim(fname),form='FORMATTED')
      do iw=1,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),dble(chi(1,i,1,j,iw))
      end do
      write(50,*)
      do iw=1,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),aimag(chi(1,i,1,j,iw))
      end do
      close(50)
    end do
  end do
! write transverse chi0 for collinear case
  if (.not.ncmag) then
    open(50,file='CHI0_T.OUT',form='FORMATTED')
    do iw=1,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(chit(iw))
    end do
    write(50,*)
    do iw=1,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(chit(iw))
      end do
    close(50)
  end if
end if
! compute f_xc in G-space
allocate(fxc(ngrf,4,ngrf,4))
call genspfxcg(fxc)
! generate the Coulomb Green's function in G+q-space regularised for q=0
allocate(gclgq(ngrf))
call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! add the regularised Coulomb interaction to f_xc to give f_Hxc
do ig=1,ngrf
  fxc(ig,1,ig,1)=fxc(ig,1,ig,1)+gclgq(ig)
end do
deallocate(gclgq)
! matrix size
n=4*ngrf
allocate(c(n,n),d(ngrf,4,ngrf,4))
! loop over frequencies
do iw=1,nwrf
! multiply f_Hxc by -chi0 from the left
  z1=-1.d0
  call zgemm('N','N',n,n,n,z1,chi(:,:,:,:,iw),n,fxc,n,zzero,c,n)
! add the identity
  do i=1,n
    c(i,i)=c(i,i)+1.d0
  end do
! invert the matrix
  call zminv(n,c)
! multiply by chi0 on the right and store in chi
  call zgemm('N','N',n,n,n,zone,c,n,chi(:,:,:,:,iw),n,zzero,d,n)
  chi(:,:,:,:,iw)=d(:,:,:,:)
end do
deallocate(c,d)
! generate transverse chi for the collinear case
if (.not.ncmag) then
  do iw=1,nwrf
    a(:,:)=chi(1,:,1,:,iw)
    call tfm13t(a,b)
    chit(iw)=b(2,2)
  end do
end if
if (mp_mpi) then
! write the complete chi matrix if required
  if (task.eq.331) then
    open(50,file='CHI.OUT',form='UNFORMATTED')
    write(50) chi
    close(50)
  end if
! write chi for G = G' = 0 in the 1x3 basis
  do i=1,4
    do j=1,4
      write(fname,'("CHI_",2I1,".OUT")') i-1,j-1
      open(50,file=trim(fname),form='FORMATTED')
      do iw=1,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),dble(chi(1,i,1,j,iw))
      end do
      write(50,*)
      do iw=1,nwrf
        write(50,'(2G18.10)') dble(wrf(iw)),aimag(chi(1,i,1,j,iw))
      end do
      close(50)
    end do
  end do
! write transverse chi for collinear case
  if (.not.ncmag) then
    open(50,file='CHI_T.OUT',form='FORMATTED')
    do iw=1,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(chit(iw))
    end do
    write(50,*)
    do iw=1,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(chit(iw))
    end do
    close(50)
  end if
  write(*,*)
  write(*,'("Info(tddftsplr):")')
  write(*,'(" Spin-dependent response function chi_ij(G,G'',q,w) written to &
   &CHI_ij.OUT")')
  write(*,'(" for i,j = 0-3; G = G'' = 0; and all wplot frequencies")')
  write(*,'(" q-vector (lattice coordinates) : ")')
  write(*,'(3G18.10)') vecql
  write(*,'(" q-vector length : ",G18.10)') gqc(1)
  write(*,*)
  write(*,'(" The elements of chi labeled by (i,j) form the 4x4 matrix :")')
  write(*,*)
  write(*,'("                  (_|_ _ _)")')
  write(*,'("  chi(G,G'',q,w) = ( |     )")')
  write(*,'("                  ( |     )")')
  write(*,'("                  ( |     )")')
  write(*,*)
  write(*,'(" (0,0) is the charge-charge response drho/dv")')
  write(*,'(" (0,1-3) is the charge-magnetisation response drho/dB")')
  write(*,'(" (1-3,0) is the magnetisation-charge response dm/v")')
  write(*,'(" (1-3,1-3) is the magnetisation-magnetisation response dm/dB")')
  write(*,*)
  write(*,'(" Non-interacting Kohn-Sham reponse function written to &
   &CHI0_ij.OUT")')
  if (.not.ncmag) then
    write(*,*)
    write(*,'(" Transverse components corresponding to m_+- = m_x +- im_y")')
    write(*,'(" written to CHI_T.OUT and CHI0_T.OUT")')
  end if
  if (task.eq.331) then
    write(*,*)
    write(*,'(" Complete response function for all G, G'' written to binary &
     &file CHI.OUT")')
    write(*,'(" (array index ordering changed from version 4.5.16 onwards)")')
  end if
end if
! write transverse response to test file
call writetest(330,'transverse response function',nv=nwrf,tol=1.d-2,zva=chit)
deallocate(gqc,ylmgq,sfacgq,chi,fxc)
if (.not.ncmag) deallocate(chit)
return

contains

subroutine tfm2213(a,b)
implicit none
! arguments
complex(8), intent(in) :: a(4,4)
complex(8), intent(out) :: b(4,4)
! local variables
integer i,j
complex(8) c(4,4),z1
do i=1,4
  c(i,1)=a(i,1)+a(i,4)
  c(i,2)=a(i,2)+a(i,3)
  z1=a(i,2)-a(i,3)
  c(i,3)=cmplx(aimag(z1),-dble(z1),8)
  c(i,4)=a(i,1)-a(i,4)
end do
do j=1,4
  b(1,j)=c(1,j)+c(4,j)
  b(2,j)=c(2,j)+c(3,j)
  z1=c(2,j)-c(3,j)
  b(3,j)=cmplx(-aimag(z1),dble(z1),8)
  b(4,j)=c(1,j)-c(4,j)
end do
return
end subroutine

subroutine tfm13t(a,b)
implicit none
! arguments
complex(8), intent(in) :: a(4,4)
complex(8), intent(out) :: b(4,4)
! local variables
integer i,j
complex(8) c(4,4),z1
do i=1,4
  c(i,1)=a(i,1)
  z1=a(i,3)
  c(i,2)=a(i,2)+cmplx(aimag(z1),-dble(z1),8)
  c(i,3)=a(i,2)+cmplx(-aimag(z1),dble(z1),8)
  c(i,4)=a(i,4)
end do
do j=1,4
  b(1,j)=c(1,j)
  z1=c(3,j)
  b(2,j)=c(2,j)+cmplx(-aimag(z1),dble(z1),8)
  b(3,j)=c(2,j)+cmplx(aimag(z1),-dble(z1),8)
  b(4,j)=c(4,j)
end do
return
end subroutine

end subroutine

