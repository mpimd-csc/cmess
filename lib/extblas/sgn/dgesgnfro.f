      SUBROUTINE DGESGNFRO(N,A,LDA,MAXIT,TOL,WORK,IWORK,WORK2,
     $                  LWORK2,VERBOSE,INFO)
      IMPLICIT NONE
      INTEGER N,MAXIT,INFO,LWORK2,LDA
      INTEGER IWORK(*)
      DOUBLE PRECISION A(*),TOL
      DOUBLE PRECISION WORK(*),WORK2(*)
      EXTERNAL DGETRI,DGETRF,DCOPY,DSCAL,DAXPY
      DOUBLE PRECISION DGEFDF,DGEFNR
      INTEGER IT,IINFO,VERBOSE,ADDIT
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
        CALL XERBLA( 'DGESGNFRO', -INFO )
        RETURN
      ENDIF
      DIF = 1

      ADDIT = 0
      DO 100 IT = 0,MAXIT
        CALL DCOPY(N*LDA,A,1,WORK,1)
        CALL DCOPY(N*LDA,A,1,WORK(LDA*N+1),1)
        CALL DGETRF( N, N, WORK, LDA, IWORK, IINFO )
        IF ( IINFO.NE.0 ) THEN
            CALL XERBLA ( 'DGETRF', -IINFO)
            RETURN
        ENDIF

        CALL DGETRI( N, WORK, LDA, IWORK, WORK2, LWORK2, IINFO )
        IF ( IINFO.NE.0 ) THEN
             CALL XERBLA ( 'DGETRFI', -IINFO)
             RETURN
         ENDIF

        D = 1
        IF ( DIF .GT. 0.1 ) THEN
            ! FROBENIUS NORM SCALING
            D = SQRT(DGEFNR(N,N,A,LDA)/DGEFNR(N,N,WORK,LDA))
        END IF

        ! Spliten weil das mit LDA nicht ganz korrekt ist.
        A(1:N*LDA) = 0.5*(A(1:N*LDA)/D + D*WORK(1:N*LDA))
        DIF = DGEFDF(N,N,A,LDA,WORK(LDA*N+1),LDA)/DGEFNR(N,N,A,LDA)
        IF (VERBOSE .NE. 0) THEN
            WRITE ( *, *) DIF, IT
        END IF
        IF ( DIF .LT. TOL )  THEN
            ADDIT = ADDIT +1
        ENDIF
        IF ( ADDIT .GE. 3 ) THEN
            GOTO 110
        END IF

100   CONTINUE
110   CONTINUE
      MAXIT = IT
      TOL = DIF
      END SUBROUTINE

