      SUBROUTINE DGETRP(N, A, LDA)
      IMPLICIT NONE
      INTEGER N,LDA
      DOUBLE PRECISION A(LDA,N)
      INTEGER I,J
      DOUBLE PRECISION T

      DO 100 I = 2,N
        DO 90 J = 1,I-1
                T=A(I,J)
                A(I,J)=A(J,I)
                A(J,I)=T
90      CONTINUE
100   CONTINUE
      RETURN
      END



