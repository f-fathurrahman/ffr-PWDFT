subroutine genidxlo
! !USES:
use modmain, only: nlorb, idxlo, idxlm, lorbl, nlotot, nlomax, idxis, &
                   natmtot, lolmmax
! !DESCRIPTION:
!   Generates an index array which maps the local-orbitals in each atom to their
!   locations in the overlap or Hamiltonian matrices. Also finds the total
!   number of local-orbitals.
implicit none
! local variables
integer is,ias,i,ilo,l,m,lm
! allocate global local-orbital index
if (allocated(idxlo)) deallocate(idxlo)
allocate(idxlo(lolmmax,nlomax,natmtot))
i=0
do ias=1,natmtot
  is=idxis(ias)
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do m=-l,l
      i=i+1
      lm=idxlm(l,m)
      idxlo(lm,ilo,ias)=i
    end do
  end do
end do
nlotot=i
return
end subroutine
!EOC

