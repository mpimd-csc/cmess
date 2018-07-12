! C = A * B
!> SUBROUTINE ddgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsenn subroutine computes C = A * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,K), intent(in) :: A
        real*8, dimension(*),intent(in) :: Bvalues
        real*8, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        ! for all columns in C
        K=K
        do j=1,N
                 !$omp parallel private(i,l,l2,temp)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                                l2 = l+1
                                temp =temp+Bvalues(l2)*A(i,Browptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsenn subroutine computes C = A * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,K), intent(in) :: A
        complex*16, dimension(*),intent(in) :: Bvalues
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        ! for all columns in C
        K=K
        do j=1,N
                 !$omp parallel private(i,l,l2,temp)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                                l2 = l+1
                                temp =temp+Bvalues(l2)*A(i,Browptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsenn subroutine computes C = A * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,K), intent(in) :: A
        real*8, dimension(*),intent(in) :: Bvalues
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        ! for all columns in C
        K=K
        do j=1,N
                 !$omp parallel private(i,l,l2,temp)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                                l2 = l+1
                                temp =temp+Bvalues(l2)*A(i,Browptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsenn subroutine computes C = A * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparsenn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,K), intent(in) :: A
        complex*16, dimension(*),intent(in) :: Bvalues
        complex*16, dimension(LDC,N), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        ! for all columns in C
        K=K
        do j=1,N
                 !$omp parallel private(i,l,l2,temp)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                                l2 = l+1
                                temp =temp+Bvalues(l2)*A(i,Browptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A * B^T
!> SUBROUTINE ddgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsent subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine ddgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
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
!> SUBROUTINE dzgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsent subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine dzgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
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
!> SUBROUTINE zdgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsent subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine zdgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
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
!> SUBROUTINE zzgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsent subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine zzgemm_dense_sparsent(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
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
! C = A * B^H
!> SUBROUTINE ddgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsenh subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine ddgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 ! forall i in col(A)
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                           l2 = l +1
                           temp =temp+((Bvalues(l2)))*A(i,Bcolptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsenh subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine dzgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 ! forall i in col(A)
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                           l2 = l +1
                           temp =temp+(dconjg((Bvalues(l2))))*A(i,Bcolptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsenh subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine zdgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 ! forall i in col(A)
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                           l2 = l +1
                           temp =temp+((Bvalues(l2)))*A(i,Bcolptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsenh subroutine computes C = A * B^T. Where A and C are Dense and
!! B is a CSR sparse matrix.
!!
subroutine zzgemm_dense_sparsenh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 ! forall i in col(A)
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                           l2 = l +1
                           temp =temp+(dconjg((Bvalues(l2))))*A(i,Bcolptr(l2)+1)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^T * B
!> SUBROUTINE ddgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsetn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Browptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsetn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Browptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsetn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Browptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsetn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparsetn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Browptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^T * B^T
!> SUBROUTINE ddgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsett subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsett subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsett subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsett subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparsett(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^T * B^H
!> SUBROUTINE ddgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparseth subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+(Bvalues(l2))*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparseth subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+dconjg((Bvalues(l2)))*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparseth subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+(Bvalues(l2))*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparseth subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparseth(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        ! C_ij = A(:,i)'*B(:,j)
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+dconjg((Bvalues(l2)))*A(Bcolptr(l2)+1,i)
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^H * B
!> SUBROUTINE ddgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsehn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*(A(Browptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsehn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*(A(Browptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsehn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*dconjg((A(Browptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsehn subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparsehn(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Bcolptr(j),Bcolptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*dconjg((A(Browptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^H * B^T
!> SUBROUTINE ddgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparseht subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*(A(Bcolptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparseht subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*(A(Bcolptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparseht subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*dconjg((A(Bcolptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparseht subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparseht(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+Bvalues(l2)*dconjg((A(Bcolptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
! C = A^H * B^H
!> SUBROUTINE ddgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The ddgemm_dense_sparsehh subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine ddgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        real*8, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        real*8 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+(Bvalues(l2))*(A(Bcolptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE dzgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The dzgemm_dense_sparsehh subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine dzgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        real*8, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+dconjg((Bvalues(l2)))*(A(Bcolptr(l2)+1,i))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zdgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zdgemm_dense_sparsehh subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zdgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        real*8, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+(Bvalues(l2))*dconjg((A(Bcolptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
!> SUBROUTINE zzgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
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
!! The zzgemm_dense_sparsehh subroutine computes C = A^T * B. Where A and C are Dense and
!! B is a CSC sparse matrix.
!!
subroutine zzgemm_dense_sparsehh(M,N,K,A,LDA,Browptr,Bcolptr,Bvalues,C,LDC)
        implicit none
        integer :: K
        integer, intent(in) :: M,N,LDA,LDC
        integer, dimension(*), intent(in) :: Bcolptr, Browptr
        complex*16, dimension(LDA,*), intent(in) :: A
        complex*16, dimension(*), intent(in) :: Bvalues
        complex*16, dimension(LDC,*), intent(out) :: C
        integer :: i,j,l,l2
        complex*16 temp
        K=K
        do j=1,N
                 !$omp parallel private(i,l,temp,l2)
                 !$omp do
                 do i = 1, M
                        temp = 0
                        do l=Browptr(j),Browptr(j+1)-1
                          l2 = l + 1
                          temp =temp+dconjg((Bvalues(l2)))*dconjg((A(Bcolptr(l2)+1,i)))
                        end do
                        C(i,j) = temp
                end do
                !$omp end do
                !$omp end parallel
        end do
end subroutine
