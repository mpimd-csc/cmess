! C = A * B
!> SUBROUTINE ddgemm_sparse_densenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densenn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densenn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,*), intent(in) :: B
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(Acolptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densenn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densenn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,*), intent(in) :: B
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(Acolptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densenn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densenn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,*), intent(in) :: B
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(Acolptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densenn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densenn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,*), intent(in) :: B
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(Acolptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A * B^T
!> SUBROUTINE ddgemm_sparse_densent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densent subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densent(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Acolptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densent subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densent(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Acolptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densent subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densent(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Acolptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densent subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densent(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Acolptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A * B^H
!> SUBROUTINE ddgemm_sparse_densenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densenh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densenh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*(B(l,Acolptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densenh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densenh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*dconjg((B(l,Acolptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densenh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densenh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*(B(l,Acolptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densenh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densenh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
        ! forall i in rows(A)
        !$omp parallel private(i,j,l,pos)
        !$omp do
        do i=1,M
                do l=1,N
                        C(i,l) = 0
                end do
                do j=Arowptr(i),Arowptr(i+1)-1
                        pos = j + 1;
                        do l = 1,N
                                C(i,l) = C(i,l) + Avalues(pos)*dconjg((B(l,Acolptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^T * B
!> SUBROUTINE ddgemm_sparse_densetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densetn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densetn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densetn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densetn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densetn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densetn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densetn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densetn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^T * B^T
!> SUBROUTINE ddgemm_sparse_densett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densett subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densett(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densett subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densett(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densett subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densett(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densett subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densett(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^T * B^H
!> SUBROUTINE ddgemm_sparse_denseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_denseth subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_denseth(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*(B(l,Arowptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_denseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_denseth subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_denseth(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*dconjg((B(l,Arowptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_denseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_denseth subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_denseth(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*(B(l,Arowptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_denseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_denseth subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_denseth(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + Avalues(pos)*dconjg((B(l,Arowptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^H * B
!> SUBROUTINE ddgemm_sparse_densehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densehn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densehn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densehn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densehn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densehn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densehn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densehn subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densehn(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*B(Arowptr(pos)+1,l)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^H * B^T
!> SUBROUTINE ddgemm_sparse_denseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_denseht subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_denseht(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_denseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_denseht subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_denseht(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_denseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_denseht subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_denseht(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_denseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_denseht subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_denseht(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*B(l,Arowptr(pos)+1)
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
! C = A^H * B^H
!> SUBROUTINE ddgemm_sparse_densehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_sparse_densehh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine ddgemm_sparse_densehh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*(B(l,Arowptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE dzgemm_sparse_densehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_sparse_densehh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine dzgemm_sparse_densehh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        real*8, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + (Avalues(pos))*dconjg((B(l,Arowptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zdgemm_sparse_densehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_sparse_densehh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zdgemm_sparse_densehh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        real*8, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*(B(l,Arowptr(pos)+1))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
!> SUBROUTINE zzgemm_sparse_densehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_sparse_densehh subroutine computes C = A * B. Where B and C are Dense and
!! A is a CSR sparse matrix.
!!
subroutine zzgemm_sparse_densehh(M,N,K,Arowptr,Acolptr,Avalues,B,LDB,C,LDC)
!subroutine dgemmsparsedensenn(M,N,K,Avalues,Acolptr,Arowptr,Bvalues,Cvalues)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDC,LDB
        integer, dimension(*), intent(in) :: Acolptr, Arowptr
        complex*16, dimension(*),intent(in) :: Avalues
        complex*16, dimension(LDB,N), intent(in) :: B
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,pos
        K=K
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
                                C(i,l) = C(i,l) + dconjg((Avalues(pos)))*dconjg((B(l,Arowptr(pos)+1)))
                        end do
                end do
        end do
        !$omp end parallel
end subroutine
