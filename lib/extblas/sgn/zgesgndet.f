      SUBROUTINE ZGESGNDET(N,A,LDA,MAXIT,TOL,WORK,IWORK,WORK2,
     $                  LWORK2,VERBOSE,INFO)
      IMPLICIT NONE
      INTEGER N,MAXIT,INFO,LWORK2,LDA
      INTEGER IWORK(*)
      DOUBLE PRECISION TOL
      COMPLEX*16 A(*)
      COMPLEX*16 WORK(*),WORK2(*)
      EXTERNAL ZGETRI,ZGETRF,ZCOPY
      DOUBLE PRECISION ZGEFDF
      INTEGER IT,IINFO,J,VERBOSE,ADDIT
      DOUBLE PRECISION DIF
      DOUBLE PRECISION D


      IINFO=0
      INFO = 0
      IF (LDA .LT. N ) THEN
        INFO=-3
      ELSE IF (N .LT. 0 ) THEN
        INFO=-1
      ELSE IF (MAXIT .LT. 1 ) THEN
        INFO=-4
      ELSE IF (TOL .LT. 0 ) THEN
        INFO=-5
      ENDIF

      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGESGNDET', -INFO )
        RETURN
      ENDIF
      DIF = 1

      ADDIT = 0
      DO 100 IT = 0,MAXIT
        CALL ZCOPY(N*LDA,A,1,WORK,1)
        CALL ZCOPY(N*LDA,A,1,WORK(LDA*N+1),1)
        CALL ZGETRF( N, N, WORK, LDA, IWORK, IINFO )
        IF ( IINFO.NE.0 ) THEN
            CALL XERBLA ( 'ZGETRF', -IINFO)
            RETURN
        ENDIF


        D = 1
        IF ( DIF .GT. 0.1 ) THEN
           ! DETERMINANT SCALING
           ! first compute determinant from LU
           ! after that compute the inverse
           D = 0
           DO J=0,N-1
               D = D+LOG(ABS(WORK(J +LDA*J+1)))
           END DO
           D = D / N
           D = EXP(D)
        END IF

        CALL ZGETRI( N, WORK, LDA, IWORK, WORK2, LWORK2, IINFO )
        IF ( IINFO.NE.0 ) THEN
           CALL XERBLA ( 'ZGETRI', -IINFO)
           RETURN
        ENDIF

        ! Spliten weil das mit LDA nicht ganz korrekt ist.
        A(1:N*LDA) = 0.5*(A(1:N*LDA)/D + D*WORK(1:N*LDA))
        DIF = ZGEFDF(N,N,A,LDA,WORK(LDA*N+1),LDA)
        IF (VERBOSE .NE. 0) THEN
            WRITE ( *, *) DIF, IT
        END IF

        IF ( DIF .LT. TOL )  THEN
            ADDIT = ADDIT +1
        ENDIF
        IF ( ADDIT .GE. 3) THEN
            GOTO 110
        ENDIF

100   CONTINUE
110   CONTINUE
      MAXIT = IT
      TOL = DIF
      END SUBROUTINE

