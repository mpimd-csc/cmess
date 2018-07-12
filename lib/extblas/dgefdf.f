      DOUBLE PRECISION FUNCTION DGEFDF(N,M,A,LDA,B,LDB)
      INTEGER N,M,LDA,LDB
      DOUBLE PRECISION A(LDA,*), B(LDB,*)
      INTEGER I,J
      DOUBLE PRECISION FNR

      FNR=0d0
      DO 100 J=1,M
       DO 110 I=1,N
        FNR = FNR+(A(I,J)-B(I,J))**2
110    CONTINUE
100   CONTINUE
      DGEFDF = SQRT(FNR)

      END FUNCTION
