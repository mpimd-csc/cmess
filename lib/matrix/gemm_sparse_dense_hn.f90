!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http:!www.gnu.org/licenses/>.
!
! Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
!               2009-2018
!

#define TRANS hn
#if TYPE == DD
      #define TYPE_STR dd
#endif
#if TYPE == DZ
      #define TYPE_STR dz
#endif
#if TYPE == ZD
      #define TYPE_STR zd
#endif
#if TYPE == ZZ
      #define TYPE_STR zz
#endif

#define FNAME _FNAME(TYPE_STR,TRANS)
#define JOIN(X,Y,Z)  X##Y##Z
#define _FNAME(X,Y)  JOIN(X,gemm_sparse_dense,Y)
#if TYPE == DZ || TYPE == ZZ
      #define CONJ(X) dconjg((X))
#else
      #define CONJ(X) (X)
#endif

#if TYPE == ZD || TYPE == ZZ
      #define CONJ2(X) dconjg((X))
#else
      #define CONJ2(X) (X)
#endif

!> SUBROUTINE FNAME(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
!! @param[in] M Number of Rows of A
!! @param[in] N Number of Columns of B
!! @param[in] K Number of Columns of A or Rows of B
!! @param[in] Arowptr rowptr of A
!! @param[in] Acolptr column pointer of A
!! @param[in] Avalues Values of A
!! @param[in] B Matrix B
!! @param[in] LDB leading dimension of B
!! @param[out] C Matrix C
!! @param[out] LDC leading dimension of C
!! @ingroup matrix_op
!!
!! The FNAME subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine FNAME(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer, intent(in) :: M,N,K,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        DATATYPE_A, dimension(*),intent(in) ::  Avalues
        DATATYPE_B, dimension(LDB,N), intent(in) :: B
        DATATYPE_C, dimension(LDC,N), intent(out) :: C

        integer :: i,j,l,pos
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Acolptr(i),Acolptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) +  CONJ2(Avalues(pos))*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine

#undef TRANS
#undef TYPE
#undef FNAME
#undef _FNAME
#undef JOIN
#undef DATATYPE_A
#undef DATATYPE_B
#undef DATATYPE_C
#undef TYPE_STR
#undef CONJ
#undef CONJ2
