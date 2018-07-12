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

#define TRANS nt
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
#define _FNAME(X,Y)  JOIN(X,gemm_dense_sparse,Y)
#if TYPE == DZ || TYPE == ZZ
      #define CONJ(X) dconjg((X))
#else
      #define CONJ(X) (X)
#endif


!> SUBROUTINE FNAME(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
!! @param[in] M Number of Rows of A
!! @param[in] N Number of Columns of B
!! @param[in] K Number of Columns of A or Rows of B
!! @param[in] A Matrix A
!! @param[in] LDA leading dimension of A
!! @param[in] Browptr rowptr of B
!! @param[in] Bcolptr column pointer of B
!! @param[in] Bvalues Values of B
!! @param[out] C Matrix C
!! @param[out] LDC leading dimension of C
!! @ingroup matrix_op
!!
!! The FNAME subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine FNAME(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer, intent(in) :: M,K,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        DATATYPE_A, dimension(LDA,*), intent(in) :: A
        DATATYPE_B, dimension(*),   intent(in) :: Bvalues
        DATATYPE_C, dimension(LDC,*), intent(out)    :: C

        integer :: i,j,l,l2
        DATATYPE_C temp
        do j=1,N
                 ! forall i in col(A)
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                           l2 = l +1
                           temp =temp+Bvalues(l2)*A(i,Bcolptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
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
