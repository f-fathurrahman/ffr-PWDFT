
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkfsm
use modmain
implicit none
! local variables
integer isym,lspn
integer is,ia,ja
real(8) sc(3,3),v(3),t1
if (fsmtype.eq.0) return
do isym=1,nsymcrys
  lspn=lspnsymc(isym)
! proper rotation matrix in Cartesian coordinates
  sc(:,:)=dble(symlatd(lspn))*symlatc(:,:,lspn)
! check invariance of global moment
  if ((abs(fsmtype).eq.1).or.(abs(fsmtype).eq.3)) then
    call r3mv(sc,momfix,v)
    t1=sum(abs(momfix(:)-v(:)))
    if (t1.gt.epslat) then
      write(*,*)
      write(*,'("Error(checkfsm): momfix not invariant under symmetry group")')
      write(*,*)
      stop
    end if
  end if
! check invariance of muffin-tin moments
  if ((abs(fsmtype).eq.2).or.(abs(fsmtype).eq.3)) then
    do is=1,nspecies
      do ia=1,natoms(is)
! if any component is >= 1000 then do not fix the moment
        t1=sum(abs(mommtfix(:,ia,is)))
        if (t1.ge.1000.d0) cycle
! equivalent atom
        ja=ieqatom(ia,is,isym)
        call r3mv(sc,mommtfix(:,ja,is),v)
        t1=sum(abs(mommtfix(:,ia,is)-v(:)))
        if (t1.gt.epslat) then
          write(*,*)
          write(*,'("Error(checkfsm): mommtfix not invariant under symmetry &
           &group")')
          write(*,'(" for species ",I4)') is
          write(*,'(" and equivalent atoms ",2I4)') ia,ja
          write(*,*)
          stop
        end if
      end do
    end do
  end if
end do
return
end subroutine

