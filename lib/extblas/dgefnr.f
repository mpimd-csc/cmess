      DOUBLE PRECISION FUNCTION DGEFNR(N,M,A,LDA)
      INTEGER N,M,LDA
      DOUBLE PRECISION A(LDA,*)
      INTEGER I,J
      DOUBLE PRECISION FNR

      FNR=0d0
      DO 100 J=1,M
       DO 110 I=1,N
        FNR = FNR+A(I,J)**2
110    CONTINUE
100   CONTINUE
      DGEFNR = SQRT(FNR)

      END FUNCTION

