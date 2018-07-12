      SUBROUTINE DGEGSGNDET(N,A,LDA,B,LDB,MAXIT,TOL,WORK,IWORK,
     $      VERBOSE, INFO)
      IMPLICIT NONE
      INTEGER N,MAXIT,INFO,LDA,LDB
      INTEGER IWORK(*)
      INTEGER VERBOSE
      DOUBLE PRECISION A(LDA,N), B(LDB,N), TOL
      DOUBLE PRECISION WORK(N,3*N)
      EXTERNAL DGETRS,DGETRF,DGEMM,DLACPY
      DOUBLE PRECISION DGEFDF, DGEFNR
      INTEGER IT,IINFO,J,ADDIT
      DOUBLE PRECISION DIF
      DOUBLE PRECISION DE,DA,GA


      IINFO=0
      INFO = 0
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
        CALL XERBLA( 'DGEGSGNDET', -INFO )
        RETURN
      ENDIF
      DIF = 1
      WORK(:,1:3*N) = 0.0

      ! Compute det(B)^1/N
      CALL DLACPY("All", N,N, B, LDB, WORK,N)
      CALL DGETRF(N,N,WORK,N,IWORK,IINFO)
      IF ( IINFO .NE. 0 ) THEN
        CALL XERBLA('DGETRF',-IINFO)
        RETURN
      ENDIF

      DE = 1d0
      DO 90 IT = 1,N
        DE = DE * DABS(WORK(IT,IT))**(1d0/N)
90    CONTINUE


      ADDIT = 0
      DO 100 IT = 0,MAXIT
         CALL DLACPY("All",N,N,A,LDA,WORK,N)
         CALL DLACPY("All",N,N,B,LDB,WORK(1,N+1),N)
         CALL DLACPY("All",N,N,A,LDA,WORK(1,2*N+1),N)
         CALL DGETRF( N, N, WORK, N, IWORK, IINFO )
         IF (INFO.NE.0) THEN
             CALL XERBLA ( 'DGETRF', -IINFO)
             RETURN
         ENDIF

         ! det(A)
         GA = 1d0
         IF (DIF .GT. 0.1) THEN
         DA = 1d0
         DO 95, J = 1,N
          DA = DA * DABS(WORK(J,J))**(1d0/N)
95      CONTINUE
         GA = DA/DE
         END IF

         CALL DGETRS( "N", N, N, WORK, N, IWORK, WORK(1,N+1), N,
     $                IINFO )
         IF (INFO.NE.0) THEN
             CALL XERBLA ( 'DGETRS', -IINFO)
             RETURN
         ENDIF

         CALL DGEMM("N","N",N,N,N,GA/2d0,B,LDB,WORK(1,N+1),N,
     $              1d0/(2d0*GA), A,LDA)

         DIF = DGEFDF(N,N,A,LDA,WORK(1,2*N+1),N)/DGEFNR(N,N,A,LDA)
         IF ( DIF .LT. TOL ) THEN
                ADDIT = ADDIT + 1
         ENDIF
         IF ( VERBOSE .NE. 0 ) THEN
          WRITE ( *, *) IT ,DIF,TOL,GA
         ENDIF
         IF ( ADDIT .GE. 3 )  THEN
                 GOTO 110
         ENDIF
100   CONTINUE
110   CONTINUE
      MAXIT = IT+1
      TOL = DIF
      END SUBROUTINE

