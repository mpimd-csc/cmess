C>    @brief Add a real matrix to a complex one.
C>    @param[in] M number of rows
C>    @param[in] N number of columns
C>    @param[in] alpha prefix factor for A (COMPLEX*16)
C>    @param[in] A  first matrix of the sum (REAL*8)
C>    @param[in]     LDA leading dimension of A
C>    @param[in]     beta prefix factor for B (COMPLEX*16)
C>    @param[in,out]     B   second matrix of the sum and the output matrix (COMPLEX*16)
C>    @param[in]     LDB leading dimension of the matrix B
C>    @ingroup fortran_blas
C>
C>    The DZGEADD subroutine adds a real matrix A to a complex matrix B.
C>    \f[ B= \alpha *A +\beta *B  \f]
C>    This subroutine is an extend to BLAS because this function is
C>    not available in all BLAS implementations.
      SUBROUTINE DZGEADD(M,N,ALPHA,A,LDA,BETA,B,LDB)
      IMPLICIT NONE
      INTEGER M,N,LDA,LDB
      COMPLEX*16 ALPHA, BETA
      DOUBLE PRECISION A(LDA,N)
      COMPLEX*16 B(LDB,N)
      INTEGER I,J
      IF ( M .LE. 0) THEN
        RETURN
      END IF
      IF ( N .LE. 0) THEN
        RETURN
      END IF
      IF ( LDA .LT. M) THEN
        RETURN
      END IF
      IF ( LDB .LT. M) THEN
        RETURN
      END IF

      IF ( ALPHA .EQ. 0.0 ) THEN
        DO J = 1,N
                DO I = 1,M
                        B(I,J) = BETA * B(I,J)
                END DO
        END DO

      ELSE IF ( BETA .EQ. 0.0 ) THEN
         DO J = 1,N
                DO I = 1,M
                        B(I,J) = ALPHA*A(I,J)
                END DO
        END DO
      ELSE
        DO J = 1,N
                DO I = 1,M
                        B(I,J) = ALPHA*A(I,J) + BETA * B(I,J)
                END DO
        END DO
      END IF
      RETURN
      END

