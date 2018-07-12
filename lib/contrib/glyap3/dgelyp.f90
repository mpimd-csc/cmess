 SUBROUTINE DGELYP( DICO, JOB, FACT, TRANA, N, A, LDA, U, LDU, C,      &
     &                  LDC, SCALE, SEP, FERR, WR, WI, IWORK, DWORK,  &
     &                  LDWORK, INFO )
!
!     PURPOSE
!
!     To solve for X either the real continuous-time Lyapunov equation
!
!        op(A)'*X + X*op(A) = scale*C                             (1)
!
!     or the real discrete-time Lyapunov equation (not implemented yet)
!
!        op(A)'*X*op(A) - X = scale*C                             (2)
!
!     and/or estimate an associated condition number, called separation,
!     where op(A) = A or A' (A**T) and C is symmetric (C = C').
!     (A' denotes the transpose of the matrix A.) A is N-by-N, the right
!     hand side C and the solution X are N-by-N, and scale is an output
!     scale factor, set less than or equal to 1 to avoid overflow in X.
!
!     ARGUMENTS
!
!     Mode Parameters
!
!     DICO    CHARACTER*1
!             Specifies the equation from which X is to be determined
!             as follows:
!             = 'C':  Equation (1), continuous-time case;
!             = 'D':  Equation (2), discrete-time case.
!
!     JOB     CHARACTER*1
!             Specifies the computation to be performed, as follows:
!             = 'X':  Compute the solution only;
!             = 'S':  Compute the separation only;
!             = 'B':  Compute both the solution and the separation.
!
!     FACT    CHARACTER*1
!             Specifies whether or not the real Schur factorization
!             of the matrix A is supplied on entry, as follows:
!             = 'F':  On entry, A and U contain the factors from the
!                     real Schur factorization of the matrix A;
!             = 'N':  The Schur factorization of A will be computed
!                     and the factors will be stored in A and U.
!
!     TRANA   CHARACTER*1
!             Specifies the form of op(A) to be used, as follows:
!             = 'N':  op(A) = A    (No transpose);
!             = 'T':  op(A) = A**T (Transpose);
!             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
!
!     Input/Output Parameters
!
!     N       (input) INTEGER
!             The order of the matrices A, X, and C.  N >= 0.
!
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the leading N-by-N part of this array must
!             contain the matrix A. If FACT = 'F', then A contains
!             an upper quasi-triangular matrix in Schur canonical form;
!             the elements below the upper Hessenberg part of the
!             array A are not referenced.
!             On exit, if INFO = 0 or INFO = N+1, the leading N-by-N
!             upper Hessenberg part of this array contains the upper
!             quasi-triangular matrix in Schur canonical form from the
!             Schur factorization of A. The contents of array A is not
!             modified if FACT = 'F'.
!
!     LDA     INTEGER
!             The leading dimension of array A.  LDA >= MAX(1,N).
!
!     U       (input or output) DOUBLE PRECISION array, dimension
!             (LDU,N)
!             If FACT = 'F', then U is an input argument and on entry
!             the leading N-by-N part of this array must contain the
!             orthogonal matrix U of the real Schur factorization of A.
!             If FACT = 'N', then U is an output argument and on exit,
!             if INFO = 0 or INFO = N+1, it contains the orthogonal
!             N-by-N matrix from the real Schur factorization of A.
!
!     LDU     INTEGER
!             The leading dimension of array U.  LDU >= MAX(1,N).
!
!     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!             On entry with JOB = 'X' or 'B', the leading N-by-N part of
!             this array must contain the symmetric matrix C.
!             On exit with JOB = 'X' or 'B', if INFO = 0 or INFO = N+1,
!             the leading N-by-N part of C has been overwritten by the
!             symmetric solution matrix X.
!             If JOB = 'S', C is not referenced.
!
!     LDC     INTEGER
!             The leading dimension of array C.
!             LDC >= 1,        if JOB = 'S';
!             LDC >= MAX(1,N), otherwise.
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor, scale, set less than or equal to 1 to
!             prevent the solution overflowing.
!
!     SEP     (output) DOUBLE PRECISION
!             If JOB = 'S' or JOB = 'B', and INFO = 0 or INFO = N+1, SEP
!             contains the estimated separation of the matrices op(A)
!             and -op(A)', if DICO = 'C' or of op(A) and op(A)', if
!             DICO = 'D'.
!             If JOB = 'X' or N = 0, SEP is not referenced.
!
!     FERR    (output) DOUBLE PRECISION
!             If JOB = 'B', and INFO = 0 or INFO = N+1, FERR contains an
!             estimated forward error bound for the solution X.
!             If XTRUE is the true solution, FERR bounds the relative
!             error in the computed solution, measured in the Frobenius
!             norm:  norm(X - XTRUE)/norm(XTRUE).
!             If JOB = 'X' or JOB = 'S', FERR is not referenced.
!
!     WR      (output) DOUBLE PRECISION array, dimension (N)
!     WI      (output) DOUBLE PRECISION array, dimension (N)
!             If FACT = 'N', and INFO = 0 or INFO = N+1, WR and WI
!             contain the real and imaginary parts, respectively, of
!             the eigenvalues of A.
!             If FACT = 'F', WR and WI are not referenced.
!
!     Workspace
!
!     IWORK   INTEGER array, dimension (N*N)
!             This array is not referenced if JOB = 'X'.
!
!     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
!             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the
!             optimal value of LDWORK.
!
!     LDWORK  INTEGER
!             The length of the array DWORK.  LDWORK >= 1, and
!             If JOB = 'X' then
!                If FACT = 'F', LDWORK >= N*N,           for DICO = 'C';
!                               LDWORK >= MAX(N*N, 2*N), for DICO = 'D';
!                If FACT = 'N', LDWORK >= MAX(N*N, 3*N).
!             If JOB = 'S' or JOB = 'B' then
!                If FACT = 'F', LDWORK >= 2*N*N,       for DICO = 'C';
!                               LDWORK >= 2*N*N + 2*N, for DICO = 'D'.
!                If FACT = 'N', LDWORK >= MAX(2*N*N, 3*N), DICO = 'C';
!                               LDWORK >= 2*N*N + 2*N, for DICO = 'D'.
!             For optimum performance LDWORK should be larger.
!
!     Error Indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             > 0:  if INFO = i, the QR algorithm failed to compute all
!                   the eigenvalues (see LAPACK Library routine DGEES);
!                   elements i+1:n of WR and WI contain eigenvalues
!                   which have converged, and A contains the partially
!                   converged Schur form;
!             = N+1:  if DICO = 'C', and the matrices A and -A' have
!                   common or very close eigenvalues, or
!                   if DICO = 'D', and matrix A has almost reciprocal
!                   eigenvalues (that is, lambda(i) = 1/lambda(j) for
!                   some i and j, where lambda(i) and lambda(j) are
!                   eigenvalues of A and i <> j); perturbed values were
!                   used to solve the equation (but the matrix A is
!                   unchanged).
!
!     METHOD
!     The algorithm uses a block Bartels-Stewart implementation.
!
!     REFERENCES
!
!     [1] Köhler M., Saak J.
!         A level 3 block variant of the Bartel-Stewart Algorithm for
!         the generalized Lyapunov equation.
!
!     This routine is based on SB03MY from SLICOT

!     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
!     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOB, TRANA
      INTEGER           INFO, LDA, LDC, LDU, LDWORK, N
      DOUBLE PRECISION  FERR, SCALE, SEP
!     .. Array Arguments ..
      INTEGER           IWORK( * )
      DOUBLE PRECISION  A( LDA, * ), C( LDC, * ), DWORK( * ), U( LDU, * ), WI( * ), WR( * )
!     .. Local Scalars ..
      LOGICAL           CONT, NOFACT, NOTA, WANTBH, WANTSP, WANTX
      CHARACTER         NOTRA, NTRNST, TRANST, UPLO
      INTEGER           I, IERR, KASE, LWA, MINWRK, NN, NN2, SDIM
      DOUBLE PRECISION  EPS, EST, SCALEF
!     .. Local Arrays ..
      LOGICAL           BWORK( 1 )
!     .. External Functions ..
      LOGICAL           LSAME, SELECT
      DOUBLE PRECISION  DLAMCH, DLANHS
      EXTERNAL          DLAMCH, DLANHS, LSAME
!     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEES, DLACON, DGEMM, DTRLYP, XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX
!     .. Executable Statements ..
!
!     Decode and Test input parameters.
!
      CONT   = LSAME( DICO,  'C' )
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      NOFACT = LSAME( FACT,  'N' )
      NOTA   = LSAME( TRANA, 'N' )
      NN  = N*N
      NN2 = 2*NN
!
      INFO = 0
      IF( .NOT.CONT .AND. .NOT.LSAME( DICO, 'D' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTBH .AND. .NOT.WANTSP .AND. .NOT.WANTX ) THEN
         INFO = -2
      ELSE IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.NOTA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( WANTSP .AND. LDC.LT.1 .OR. .NOT.WANTSP .AND. LDC.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE
         IF ( WANTX ) THEN
            IF ( NOFACT ) THEN
               MINWRK = MAX( NN, 3*N )
            ELSE IF ( CONT ) THEN
               MINWRK = NN
            ELSE
               MINWRK = MAX( NN, 2*N )
            END IF
         ELSE
            IF ( CONT ) THEN
               IF ( NOFACT ) THEN
                  MINWRK = MAX( NN2, 3*N )
               ELSE
                  MINWRK = NN2
               END IF
            ELSE
               MINWRK = NN2 + 2*N
            END IF
         END IF
!       Workspace Query
        IF( LDWORK .EQ. -1 ) THEN
            DWORK(1) = DBLE(MINWRK)
            INFO = 0
            RETURN
        END IF

        IF( LDWORK.LT.MAX( 1, MINWRK ) )  INFO = -19
      END IF
!
      IF ( INFO.NE.0 ) THEN
!        Error return.
         CALL XERBLA( 'DGELYP', -INFO )
         RETURN
      END IF
      IF ( .NOT. CONT ) THEN
          INFO = -1
          CALL XERBLA( 'DGELYP', -INFO )
          RETURN
      END IF

!
!     Quick return if possible.
!
      IF( N.EQ.0 ) THEN
         SCALE = ONE
         IF( WANTBH ) FERR  = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
!
      LWA = 0
!
      IF( NOFACT ) THEN
!
!        Compute the Schur factorization of A.
!        Workspace:  need   3*N;
!                    prefer larger.
!        (Note: Comments in the code beginning "Workspace:" describe the
!        minimal amount of real workspace needed at that point in the
!        code, as well as the preferred amount for good performance.
!        NB refers to the optimal block size for the immediately
!        following subroutine, as returned by ILAENV.)
!
         CALL DGEES( 'Vectors', 'Not ordered', SELECT, N,A, LDA, SDIM,  WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         IF( INFO.GT.0 )  RETURN
         LWA = INT( DWORK( 1 ) )
      END IF
!
      IF( .NOT.WANTSP ) THEN
!
!        Transform the right-hand side.
!        Workspace:  N*N.
!
         NTRNST = 'N'
         TRANST = 'T'
         UPLO   = 'U'

         CALL DGEMM(TRANST, 'N', N, N, N, ONE, U, LDU, C, LDC, ZERO, DWORK, N)
         CALL DGEMM('N', NTRNST, N, N, N, ONE, DWORK, N, U, LDU, ZERO, C, LDC)
!
         DO I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
         END DO
!
         LWA = MAX( LWA, NN )
!
!        Solve the transformed equation.
!        Workspace for DICO = 'D':  2*N.
!
         IF ( CONT ) THEN
            CALL DTRLYP( TRANA, N, A, LDA, C, LDC, SCALE, INFO )
         ELSE
            ! CALL DTRSTN( TRANA, N, A, LDA, C, LDC, SCALE, DWORK, INFO )
         END IF
         IF( INFO.GT.0 )  INFO = N + 1
!
!        Transform back the solution.
!        Workspace:  N*N.
!
         CALL DGEMM(NTRNST, 'N', N, N, N, ONE, U, LDU, C, LDC, ZERO, DWORK, N)
         CALL DGEMM('N', TRANST, N, N, N, ONE, DWORK, N, U, LDU, ZERO, C, LDC)
!
         DO I = 2, N
            CALL DCOPY( I-1, C(1,I), 1, C(I,1), LDC )
         END DO
      END IF
!
      IF( .NOT.WANTX ) THEN
!
!        Estimate the separation.
!        Workspace:  2*N*N       for DICO = 'C';
!                    2*N*N + 2*N for DICO = 'D'.
!
         IF( NOTA ) THEN
            NOTRA = 'T'
         ELSE
            NOTRA = 'N'
         END IF
!
         EST = ZERO
         KASE = 0
!        REPEAT
   30    CONTINUE
         CALL DLACON( NN, DWORK(NN+1), DWORK, IWORK, EST, KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
               IF( CONT ) THEN
                  CALL DTRLYP( TRANA, N, A, LDA, DWORK, N, SCALEF, IERR )
               ELSE
                  ! CALL DTRSTN( TRANA, N, A, LDA, DWORK, N, SCALEF, DWORK(NN2+1), IERR )
               END IF
            ELSE
               IF( CONT ) THEN
                  CALL DTRLYP( NOTRA, N, A, LDA, DWORK, N, SCALEF, IERR )
               ELSE
                  ! CALL DTRSTN( NOTRA, N, A, LDA, DWORK, N, SCALEF, DWORK(NN2+1), IERR )
               END IF
            END IF
            GO TO 30
         END IF
!        UNTIL KASE = 0
!
         SEP = SCALEF / EST
!
         IF( WANTBH ) THEN
!
!           Get the machine precision.
!
            EPS = DLAMCH( 'P' )
!
!           Compute the estimate of the relative error.
!
            IF ( CONT ) THEN
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )/SEP
            ELSE
               FERR = EPS*DLANHS( 'Frobenius', N, A, LDA, DWORK )**2/SEP
            END IF
         END IF
      END IF
!
      DWORK( 1 ) = DBLE( MAX( LWA, MINWRK ) )
      RETURN
END SUBROUTINE
