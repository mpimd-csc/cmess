      SUBROUTINE ZGEGSGNNONE(N,A,LDA,B,LDB,MAXIT,TOL,WORK,IWORK,
     $      VERBOSE, INFO)
      IMPLICIT NONE
      INTEGER N,MAXIT,INFO,LDA,LDB
      INTEGER IWORK(*)
      INTEGER VERBOSE
      DOUBLE PRECISION TOL
      COMPLEX*16 A(LDA,N), B(LDB,N), WORK(N,3*N)
      EXTERNAL ZGETRS,ZGETRF,ZGEMM,ZLACPY
      DOUBLE PRECISION ZGEFDF, ZGEFNR
      INTEGER IT,IINFO,ADDIT
      DOUBLE PRECISION DIF
      COMPLEX*16 GA

      GA= 1d0/2d0
      INFO = 0
      IINFO=0
      IF (N .LT. 0 ) THEN
        INFO=-1
      ELSE IF ( LDA .LT. MAX(1,N)) THEN
        INFO=-3
      ELSE IF ( LDB .LT. MAX(1,N)) THEN
        INFO=-5
      ELSE IF (MAXIT .LT. 1 ) THEN
        INFO=-6
      ELSE IF (TOL .LT. 0 ) THEN
        INFO=-7
      ENDIF

      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGEGSGNNONE', -INFO )
        RETURN
      ENDIF
      DIF = 1
      WORK(:,1:3*N) = 0.0

      ADDIT = 0
      DO 100 IT = 0,MAXIT
         CALL ZLACPY("All",N,N,A,LDA,WORK,N)
         CALL ZLACPY("All",N,N,B,LDB,WORK(1,N+1),N)
         CALL ZLACPY("All",N,N,A,LDA,WORK(1,2*N+1),N)
         CALL ZGETRF( N, N, WORK, N, IWORK, IINFO )
         IF (INFO.NE.0) THEN
             CALL XERBLA ( 'ZGETRF', -IINFO)
             RETURN
         ENDIF

         CALL ZGETRS( "N", N, N, WORK, N, IWORK, WORK(1,N+1), N,
     $                IINFO )
         IF (INFO.NE.0) THEN
             CALL XERBLA ( 'ZGETRS', -IINFO)
             RETURN
         ENDIF

         CALL ZGEMM("N","N",N,N,N,GA,B,LDB,WORK(1,N+1),N,
     $              GA, A,LDA)

         DIF = ZGEFDF(N,N,A,LDA,WORK(1,2*N+1),N)/ZGEFNR(N,N,A,LDA)
         IF ( DIF .LT. TOL ) THEN
                ADDIT = ADDIT + 1
         ENDIF
         IF ( VERBOSE .NE. 0 ) THEN
          WRITE ( *, *) IT ,DIF,TOL
         ENDIF
         IF ( ADDIT .GE. 3 )  THEN
                 GOTO 110
         ENDIF
100   CONTINUE
110   CONTINUE
      MAXIT = IT+1
      TOL = DIF
      END SUBROUTINE

