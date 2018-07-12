SUBROUTINE DGGLYP( DICO, JOB, FACT, TRANS, UPLO, ISOLVE, NB, N, A, LDA, E, &
                 & LDE, Q, LDQ, Z, LDZ, X, LDX, SCALE, SEP, FERR, &
                 &  ALPHAR, ALPHAI, BETA, IWORK, DWORK, LDWORK,&
                 &  INFO )
!
!     PURPOSE
!
!     To solve for X either the generalized continuous-time Lyapunov
!     equation
!
!             T                T
!        op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1)
!
!     or the generalized discrete-time Lyapunov equation
!
!             T                T
!        op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2)
!
!     where op(M) is either M or M**T for M = A, E and the right hand
!     side Y is symmetric. A, E, Y, and the solution X are N-by-N
!     matrices. SCALE is an output scale factor, set to avoid overflow
!     in X.
!
!
!
!     Estimates of the separation and the relative forward error norm
!     are provided.
!
!     The function is based on SG03AD, SG03AX and SG03AY from SLICOT.
!     In contrast to the SLICOT variant this implementation used a blocked
!     LEVEL-3 BLAS implemenatation of the Bartels-Stewart algorithm in order
!     to achieve a higher performance on modern computer architectures.
!
!     ARGUMENTS
!
!     Mode Parameters
!
!     DICO    CHARACTER*1
!             Specifies which type of the equation is considered:
!             = 'C':  Continuous-time equation (1);
!             = 'D':  Discrete-time equation (2).
!
!     JOB     CHARACTER*1
!             Specifies if the solution is to be computed and if the
!             separation is to be estimated:
!             = 'X':  Compute the solution only;
!             = 'S':  Estimate the separation only;
!             = 'B':  Compute the solution and estimate the separation.
!
!     FACT    CHARACTER*1
!             Specifies whether the generalized real Schur
!             factorization of the pencil A - lambda * E is supplied
!             on entry or not:
!             = 'N':  Factorization is not supplied;
!             = 'F':  Factorization is supplied.
!
!     TRANS   CHARACTER*1
!             Specifies whether the transposed equation is to be solved
!             or not:
!             = 'N':  op(A) = A,    op(E) = E;
!             = 'T':  op(A) = A**T, op(E) = E**T.
!
!     UPLO    CHARACTER*1
!             Specifies whether the lower or the upper triangle of the
!             array X is needed on input:
!             = 'L':  Only the lower triangle is needed on input;
!             = 'U':  Only the upper triangle is needed on input.
!     ISOLVE  INTEGER
!             Specifies which solver should be used to solve the small
!             generalized Sylvester equations arising in the blocked
!             algorithm
!             = +/-1 : Use the Gardiner/Laub variant of Bartels-Stewart
!                      employing Forward/Backward Substition with PLUQ
!                      decompositions on the diagonal blocks.
!             = +/-2 : Use the Kagström/Westin generalized Schur algorithm
!                      to solve the small subsystem
!             If ISOLVE is -1 or -2 it uses  the above solvers as well
!             but it solves for the diagonal blocks with Lyapunov solver.
!             Only ISOLVE .EQ. 2 allows to compute the scaling factor
!             SCALE.
!
!    NB       INTEGER
!             Block size in the Bartels-Stewart algorithm. This value
!             must be set an optimal value for the used computer
!             architecture and the used BLAS library. A good value
!             if ISOLVE .EQ. 3 might be 32, 48 or 64 otherwise
!             16, 32 or 48 are good choices. NB must be even and
!             NB must be larger than 1.
!
!
!    Input/Output Parameters
!
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, if FACT = 'F', then the leading N-by-N upper
!             Hessenberg part of this array must contain the
!             generalized Schur factor A_s of the matrix A (see
!             definition (3) in section METHOD). A_s must be an upper
!             quasitriangular matrix. The elements below the upper
!             Hessenberg part of the array A are not referenced.
!             If FACT = 'N', then the leading N-by-N part of this
!             array must contain the matrix A.
!             On exit, the leading N-by-N part of this array contains
!             the generalized Schur factor A_s of the matrix A. (A_s is
!             an upper quasitriangular matrix.)
!
!     LDA     INTEGER
!             The leading dimension of the array A.  LDA >= MAX(1,N).
!
!     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
!             On entry, if FACT = 'F', then the leading N-by-N upper
!             triangular part of this array must contain the
!             generalized Schur factor E_s of the matrix E (see
!             definition (4) in section METHOD). The elements below the
!             upper triangular part of the array E are not referenced.
!             If FACT = 'N', then the leading N-by-N part of this
!             array must contain the coefficient matrix E of the
!             equation.
!             On exit, the leading N-by-N part of this array contains
!             the generalized Schur factor E_s of the matrix E. (E_s is
!             an upper triangular matrix.)
!
!     LDE     INTEGER
!             The leading dimension of the array E.  LDE >= MAX(1,N).
!
!     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!             On entry, if FACT = 'F', then the leading N-by-N part of
!             this array must contain the orthogonal matrix Q from
!             the generalized Schur factorization (see definitions (3)
!             and (4) in section METHOD).
!             If FACT = 'N', Q need not be set on entry.
!             On exit, the leading N-by-N part of this array contains
!             the orthogonal matrix Q from the generalized Schur
!             factorization.
!
!     LDQ     INTEGER
!             The leading dimension of the array Q.  LDQ >= MAX(1,N).
!
!     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!             On entry, if FACT = 'F', then the leading N-by-N part of
!             this array must contain the orthogonal matrix Z from
!             the generalized Schur factorization (see definitions (3)
!             and (4) in section METHOD).
!             If FACT = 'N', Z need not be set on entry.
!             On exit, the leading N-by-N part of this array contains
!             the orthogonal matrix Z from the generalized Schur
!             factorization.
!
!     LDZ     INTEGER
!             The leading dimension of the array Z.  LDZ >= MAX(1,N).
!
!     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
!             On entry, if JOB = 'B' or 'X', then the leading N-by-N
!             part of this array must contain the right hand side matrix
!             Y of the equation. Either the lower or the upper
!             triangular part of this array is needed (see mode
!             parameter UPLO).
!             If JOB = 'S', X is not referenced.
!             On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then
!             the leading N-by-N part of this array contains the
!             solution matrix X of the equation.
!             If JOB = 'S', X is not referenced.
!
!     LDX     INTEGER
!             The leading dimension of the array X.  LDX >= MAX(1,N).
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor set to avoid overflow in X.
!             (0 < SCALE <= 1)
!             It is only computed iff ISOLVE .EQ. 3 otherwise it is set to 1.
!
!     SEP     (output) DOUBLE PRECISION
!             If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then
!             SEP contains an estimate of the separation of the
!             Lyapunov operator.
!
!     FERR    (output) DOUBLE PRECISION
!             If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an
!             estimated forward error bound for the solution X. If XTRUE
!             is the true solution, FERR estimates the relative error
!             in the computed solution, measured in the Frobenius norm:
!             norm(X - XTRUE) / norm(XTRUE)
!
!     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
!     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
!     BETA    (output) DOUBLE PRECISION array, dimension (N)
!             If FACT = 'N' and INFO = 0, 3, or 4, then
!             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the
!             eigenvalues of the matrix pencil A - lambda * E.
!             If FACT = 'F', ALPHAR, ALPHAI, and BETA are not
!             referenced.
!
!     Workspace
!
!     IWORK   INTEGER array, dimension (T)
!             The following table contains the required integer work space T
!             depending on the configured solver:
!
!                  JOB       ISOLVE   |  T
!                 --------------------+---------
!                  'X'       1        | MAX(2*(NB+1),1)
!                  'X'       2        | 1
!                  'B', 'S'  1        | MAX(N**2+2*(NB+1),1)
!                  'B', 'S'  2        | MAX(N**2,1)
!
!     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
!             On exit, if INFO = 0, DWORK(1) returns the optimal value
!             of LDWORK.
!
!     LDWORK  INTEGER
!             The length of the array DWORK. The following table
!             contains the minimal work space requirements depending
!             on the choice of JOB, FACT and ISOLVE.
!
!                    JOB        FACT  ISOLVE  |  LDWORK
!                    -------------------------+-------------------
!                    'X'        'F'   1       |  MAX(1,N**2,(4*(NB+1)*(NB+2)))
!                    'X'        'F'   2       |  MAX(1,N**2,(NB+1)**2 )
!                    'X'        'N'   1       |  MAX(32,N**2,(4*(NB+1)*(NB+2)), 8*N+16 )
!                    'X'        'N'   2       |  MAX(32,N**2,(NB+1)**2, 8*N+16 )
!                    'B','S'    'F'   1       |  MAX(1,2*N**2+(4*(NB+1)*(NB+2)))
!                    'B','S'    'F'   2       |  MAX(1,2*N**2+(NB+1)**2 )
!                    'B','S'    'N'   1       |  MAX(32,2*N**2+(4*(NB+1)*(NB+2)), 8*N+16 )
!                    'B','S'    'N'   2       |  MAX(32,2*N**2+(NB+1)**2, 8*N+16 )
!
!     Error indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  FACT = 'F' and the matrix contained in the upper
!                   Hessenberg part of the array A is not in upper
!                   quasitriangular form;
!             = 2:  FACT = 'N' and the pencil A - lambda * E cannot be
!                   reduced to generalized Schur form: LAPACK routine
!                   DGGES has failed to converge;
!             = 3:  DICO = 'D' and the pencil A - lambda * E has a
!                   pair of reciprocal eigenvalues. That is, lambda_i =
!                   1/lambda_j for some i and j, where lambda_i and
!                   lambda_j are eigenvalues of A - lambda * E. Hence,
!                   equation (2) is singular;  perturbed values were
!                   used to solve the equation (but the matrices A and
!                   E are unchanged);
!             = 4:  DICO = 'C' and the pencil A - lambda * E has a
!                   degenerate pair of eigenvalues. That is, lambda_i =
!                   -lambda_j for some i and j, where lambda_i and
!                   lambda_j are eigenvalues of A - lambda * E. Hence,
!                   equation (1) is singular;  perturbed values were
!                   used to solve the equation (but the matrices A and
!                   E are unchanged).
!
!     METHOD
!
!     A straightforward generalization [3] of the method proposed by
!     Bartels and Stewart [1] is utilized to solve (1) or (2).
!
!     First the pencil A - lambda * E is reduced to real generalized
!     Schur form A_s - lambda * E_s by means of orthogonal
!     transformations (QZ-algorithm):
!
!        A_s = Q**T * A * Z   (upper quasitriangular)                (3)
!
!        E_s = Q**T * E * Z   (upper triangular).                    (4)
!
!     If FACT = 'F', this step is omitted. Assuming SCALE = 1 and
!     defining
!
!              ( Z**T * Y * Z   :   TRANS = 'N'
!        Y_s = <
!              ( Q**T * Y * Q   :   TRANS = 'T'
!
!
!              ( Q**T * X * Q    if TRANS = 'N'
!        X_s = <                                                     (5)
!              ( Z**T * X * Z    if TRANS = 'T'
!
!     leads to the reduced Lyapunov equation
!
!               T                      T
!        op(A_s)  X_s op(E_s) + op(E_s)  X_s op(A_s) = Y_s,          (6)
!
!     or
!               T                      T
!        op(A_s)  X_s op(A_s) - op(E_s)  X_s op(E_s) = Y_s,          (7)
!
!     which are equivalent to (1) or (2), respectively. The solution X_s
!     of (6) or (7) is computed via block back substitution (if TRANS =
!     'N') or block forward substitution (if TRANS = 'T'), where the
!     block order is NB or NB+1 depending if a pair of complex eigenvalues
!     is at the border of a block. (See [1], [4] and [3] for details.)
!     Equation (5) yields the solution matrix X.
!
!     For fast computation the estimates of the separation and the
!     forward error are gained from (6) or (7) rather than (1) or
!     (2), respectively. We consider (6) and (7) as special cases of the
!     generalized Sylvester equation
!
!        R * X * S + U * X * V = Y,                                  (8)
!
!     whose separation is defined as follows
!
!        sep = sep(R,S,U,V) =   min   || R * X * S + U * X * V || .
!                            ||X|| = 1                           F
!                                 F
!
!     Equation (8) is equivalent to the system of linear equations
!
!        K * vec(X) = (kron(S**T,R) + kron(V**T,U)) * vec(X) = vec(Y),
!
!     where kron is the Kronecker product of two matrices and vec
!     is the mapping that stacks the columns of a matrix. If K is
!     nonsingular then
!
!        sep = 1 / ||K**(-1)|| .
!                             2
!
!     We estimate ||K**(-1)|| by a method devised by Higham [2]. Note
!     that this method yields an estimation for the 1-norm but we use it
!     as an approximation for the 2-norm. Estimates for the forward
!     error norm are provided by
!
!        FERR = 2 * EPS * ||A_s||  * ||E_s||  / sep
!                                F          F
!
!     in the continuous-time case (1) and
!
!        FERR = EPS * ( ||A_s|| **2 + ||E_s|| **2 ) / sep
!                              F             F
!
!     in the discrete-time case (2).
!     The reciprocal condition number, RCOND, of the Lyapunov equation
!     can be estimated by FERR/EPS.
!
!     REFERENCES
!
!     [1] Bartels, R.H., Stewart, G.W.
!         Solution of the equation A X + X B = C.
!         Comm. A.C.M., 15, pp. 820-826, 1972.
!
!     [2] Higham, N.J.
!         FORTRAN codes for estimating the one-norm of a real or complex
!         matrix, with applications to condition estimation.
!         A.C.M. Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, 1988.
!
!     [3] Penzl, T.
!         Numerical solution of generalized Lyapunov equations.
!         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
!
!     [4] Köhler M., Saak J.
!         A level 3 block variant of the Bartel-Stewart Algorithm for
!         the generalized Lyapunov equation.
!         to be written
!
!     NUMERICAL ASPECTS
!
!     The number of flops required by the routine is given by the
!     following table. Note that we count a single floating point
!     arithmetic operation as one flop. c is an integer number of modest
!     size (say 4 or 5). The following table contains the flop rate
!     if NB = 1 and ISOLVE = 3.
!
!                   |  FACT = 'F'            FACT = 'N'
!        -----------+------------------------------------------
!        JOB = 'B'  |  (26+8*c)/3 * N**3     (224+8*c)/3 * N**3
!        JOB = 'S'  |  8*c/3 * N**3          (198+8*c)/3 * N**3
!        JOB = 'X'  |  26/3 * N**3           224/3 * N**3
!
!     The algorithm is backward stable if the eigenvalues of the pencil
!     A - lambda * E are real.   Otherwise, linear systems of order at
!     most 4 are involved into the computation. These systems are solved
!     by Gauss elimination with complete pivoting. The loss of stability
!     of the Gauss elimination with complete pivoting is rarely
!     encountered in practice.
!
!     The Lyapunov equation may be very ill-conditioned. In particular,
!     if DICO = 'D' and the pencil A - lambda * E has a pair of almost
!     reciprocal eigenvalues, or DICO = 'C' and the pencil has an almost
!     degenerate pair of eigenvalues, then the Lyapunov equation will be
!     ill-conditioned. Perturbed values were used to solve the equation.
!     Ill-conditioning can be detected by a very small value of the
!     reciprocal condition number RCOND.
!
!     AUTHOR
!
!     Köhler M., MPI Magdeburg
!
!
!     KEYWORDS
!
!     Lyapunov equation
!
!     ******************************************************************
!
!     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         ( ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 )
!     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOB, TRANS, UPLO
      DOUBLE PRECISION  FERR, SCALE, SEP
      INTEGER           INFO, LDA, LDE, LDQ, LDWORK, LDX, LDZ, N, NB, ISOLVE
!     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ALPHAI(*), ALPHAR(*), BETA(*),&
     &                  DWORK(*), E(LDE,*), Q(LDQ,*), X(LDX,*),&
     &                  Z(LDZ,*)
      INTEGER           IWORK(*)
!     .. Local Scalars ..
      CHARACTER         ETRANS
      DOUBLE PRECISION  EST, EPS, NORMA, NORME, SCALE1
      INTEGER           I, INFO1, KASE, MINWRK, SDIM,MINWRKI
      LOGICAL           ISDISC, ISFACT, ISTRAN, ISUPPR, WANTBH, WANTSP, WANTX
      LOGICAL           DUMMY

!     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DNRM2, LSAME
!     .. External Subroutines ..
      EXTERNAL          DCOPY, DGGES, DLACON,DTGLYP, XERBLA,DGEMM,DSYMM,DTGSTN
      EXTERNAL          DLACPY, DAXPY, DTGEXB
!     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, ABS, MOD
!     .. Executable Statements ..
!
!     Decode input parameters.
!
      ISDISC = LSAME( DICO,  'D' )
      WANTX  = LSAME( JOB,   'X' )
      WANTSP = LSAME( JOB,   'S' )
      WANTBH = LSAME( JOB,   'B' )
      ISFACT = LSAME( FACT,  'F' )
      ISTRAN = LSAME( TRANS, 'T' )
      ISUPPR = LSAME( UPLO,  'U' )
      DUMMY = .TRUE.

!
!     Check the scalar input parameters.
!
      IF ( .NOT.( ISDISC .OR. LSAME( DICO,  'C' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( WANTX .OR. WANTSP .OR. WANTBH ) ) THEN
         INFO = -2
      ELSEIF ( .NOT.( ISFACT .OR. LSAME( FACT,  'N' ) ) ) THEN
         INFO = -3
      ELSEIF ( .NOT.( ISTRAN .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -4
      ELSEIF ( .NOT.( ISUPPR .OR. LSAME( UPLO,  'L' ) ) ) THEN
         INFO = -5
      ELSEIF ( ABS(ISOLVE) .LT. 1 .OR. ABS(ISOLVE) .GT. 3) THEN
         INFO = -6
      ELSEIF ( NB .LT. 1 .OR. MOD(NB,2) .EQ. 1 ) THEN
         INFO = -7
      ELSEIF ( N .LT. 0 ) THEN
         INFO = -8
      ELSEIF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -10
      ELSEIF ( LDE .LT. MAX( 1, N ) ) THEN
         INFO = -12
      ELSEIF ( LDQ .LT. MAX( 1, N ) ) THEN
         INFO = -14
      ELSEIF ( LDZ .LT. MAX( 1, N ) ) THEN
         INFO = -16
      ELSEIF ( LDX .LT. MAX( 1, N ) ) THEN
         INFO = -18
      ELSE
         INFO = 0
      END IF
      IF ( INFO .EQ. 0 ) THEN
         IF ( WANTX ) THEN
            IF ( ISFACT ) THEN
                IF ( ABS(ISOLVE).EQ.1 ) THEN
                    MINWRK = MAX(1,N*N,(2*(NB+1))*(2*(NB+2)))
                ELSE
                    MINWRK = MAX(1,N*N,(NB+1)*(NB+1) )
                END IF
            ELSE
                IF ( ABS(ISOLVE).EQ.1) THEN
                    MINWRK = MAX(32,N*N,(2*(NB+1))*(2*(NB+2)),8*N+16)
                ELSE
                    MINWRK = MAX(32,N*N,(NB+1)*(NB+1),8*N+16)
                END IF
            END IF

            IF ( ABS(ISOLVE) .EQ. 1) THEN
               MINWRKI = MAX(1,2*(NB+1))
            ELSE
               MINWRKI = 1
            END IF
         ELSE
            IF ( ISFACT ) THEN
                IF ( ABS(ISOLVE) .EQ. 1) THEN
                    MINWRK = MAX(1,2*N*N+(2*(NB+1))*(2*(NB+2)))
                ELSE
                    MINWRK = MAX(1,2*N*N+(NB+1)*(NB+1) )
                END IF
            ELSE
                IF ( ABS(ISOLVE) .EQ. 1) THEN
                    MINWRK = MAX(32,2*N*N+(2*(NB+1))*(2*(NB+2)),8*N+16)
                ELSE
                    MINWRK = MAX(32,2*N*N+(NB+1)*(NB+1),8*N+16)
                END IF
            END IF

            IF ( ABS(ISOLVE) .EQ. 1) THEN
               MINWRKI = MAX(1,N*N + 2*(NB+1))
            ELSE
               MINWRKI = MAX(1,N*N)
            END IF
         END IF
         IF ( LDWORK .EQ. -1) THEN
                 DWORK(1) = DBLE(MINWRK)
                 IWORK(1) = MINWRKI
                 INFO = 0
                 RETURN
         END IF
         IF ( MINWRK .GT. LDWORK ) THEN
            INFO = -27
         END IF
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DGGLYP', -INFO )
         RETURN
      END IF

!
!     Quick return if possible.
!
      IF ( N .EQ. 0 ) THEN
         SCALE = ONE
         IF ( .NOT.WANTX ) SEP = ZERO
         IF ( WANTBH ) FERR = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
!
      IF ( ISFACT ) THEN
!
!        Make sure the upper Hessenberg part of A is quasitriangular.
!
         DO 20 I = 1, N-2
            IF ( A(I+1,I).NE.ZERO .AND. A(I+2,I+1).NE.ZERO ) THEN
               INFO = 1
               RETURN
            END IF
   20    CONTINUE
      END IF
!
      IF ( .NOT.ISFACT ) THEN
!
!        Reduce A - lambda * E to generalized Schur form.
!
!           A := Q**T * A * Z   (upper quasitriangular)
!           E := Q**T * E * Z   (upper triangular)
!
!        ( Workspace: >= MAX(1,8*N+16,32) )
!
         CALL DGGES( 'Vectors', 'Vectors','NoSort',DUMMY, N, A, LDA, E, LDE, SDIM, &
                   & ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, DWORK, LDWORK,DUMMY, INFO1 )
         IF ( INFO1 .NE. 0 ) THEN
            INFO = 2
            RETURN
         END IF

         IF ( MOD(NB,2) .EQ. 0) THEN
             CALL DTGEXB(N, A, LDA, E, LDE, Q, LDQ, Z, LDZ, NB, DWORK, LDWORK, INFO1)
             IF ( INFO1 .NE. 0 ) THEN
                 INFO = 2
                 RETURN
             END IF
         END IF

      END IF
!
      IF ( WANTBH .OR. WANTX ) THEN
!
!        Transform right hand side.
!
!           X := Z**T * X * Z  or  X := Q**T * X * Q
!
!
!        ( Workspace: >= N )
!! SUBROUTINE DSYORU(TRANS, UPLO, M, X, LDX, Q, LDQ, WORK, LDWORK, INFO)
         IF ( ISTRAN ) THEN
                ! W <- XQ
                CALL DSYMM('Left',UPLO, N,N,ONE,X,LDX,Q,LDQ,ZERO,DWORK,N)
                ! X <- Q^TW
                CALL DGEMM('T','N',N,N,N,ONE, Q, LDQ,DWORK, N, ZERO, X, LDX)
         ELSE
                ! W <- XZ
                CALL DSYMM('Left',UPLO, N,N,ONE,X,LDX,Z,LDZ,ZERO,DWORK,N)
                ! X <- Z^TW
                CALL DGEMM('T','N',N,N,N,ONE, Z, LDZ,DWORK, N, ZERO, X, LDX)
         END IF
!
!        Solve reduced generalized Lyapunov equation.
!
         IF ( ISDISC ) THEN
            CALL DTGSTN( TRANS, ISOLVE, NB, N, A, LDA, E, LDE, X, LDX, SCALE, DWORK, IWORK, INFO1)
            IF ( INFO1 .NE. 0 ) THEN
                    INFO = 3
                    RETURN
            END IF
         ELSE
            CALL DTGLYP( TRANS, ISOLVE, NB, N, A, LDA, E, LDE, X, LDX, SCALE, DWORK, IWORK, INFO1)
            IF ( INFO1 .NE. 0 ) THEN
                    INFO = 4
                    RETURN
            END IF
         END IF
!
!        Transform the solution matrix back.
!
!           X := Q * X * Q**T  or  X := Z * X * Z**T.
!
!        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2.
!
!        ( Workspace: >= N )
!
         IF ( ISTRAN ) THEN
                  ! W <- ZX
                 CALL DSYMM('Right','Upper', N,N,ONE,X,LDX,Z,LDZ,ZERO,DWORK,N)
                 ! X <- XZ^T
                 CALL DGEMM('N','T',N,N,N,ONE, DWORK, N, Z, LDZ, ZERO, X, LDX)
         ELSE
                 ! W <- QX
                 CALL DSYMM('Right','Upper', N,N,ONE,X,LDX,Q,LDQ,ZERO,DWORK,N)
                 ! X <- XQ^T
                 CALL DGEMM('N','T',N,N,N,ONE, DWORK, N, Q, LDQ, ZERO, X, LDX)
         END IF

         ! Not used in Benchmarks
         DO I = 1, N-1
           CALL DCOPY( N-I, X(I,I+1), LDX, X(I+1,I), 1 )
         END DO
      END IF !Both or X
!
      IF ( WANTBH .OR. WANTSP ) THEN
!
!        Estimate the 1-norm of the inverse Kronecker product matrix
!        belonging to the reduced generalized Lyapunov equation.
!
!        ( Workspace: 2*N*N )
!
         EST = ZERO
         KASE = 0
   80    CONTINUE
         CALL DLACON( N*N, DWORK(N*N+1), DWORK, IWORK, EST, KASE )
         IF ( KASE .NE. 0 ) THEN
            IF ( ( KASE.EQ.1 .AND. .NOT.ISTRAN ) .OR. ( KASE.NE.1 .AND. ISTRAN ) ) THEN
               ETRANS = 'N'
            ELSE
               ETRANS = 'T'
            END IF
            IF ( ISDISC ) THEN
               CALL DTGSTN( ETRANS,ISOLVE,NB, N, A, LDA, E, LDE, DWORK, N,SCALE1,DWORK(2*N*N+1),IWORK(N*N+1), INFO1 )
               IF ( INFO1 .NE. 0 ) INFO = 3
            ELSE
               CALL DTGLYP( ETRANS,ISOLVE,NB, N, A, LDA, E, LDE, DWORK, N,SCALE1,DWORK(2*N*N+1),IWORK(N*N+1), INFO1 )
               IF ( INFO1 .NE. 0 ) INFO = 4
            END IF
         GOTO 80
         END IF
         SEP = SCALE1/EST
      END IF
!
!     Estimate the relative forward error.
!
!     ( Workspace: 2*N )
!
      IF ( WANTBH .OR. WANTSP ) THEN
         EPS = DLAMCH( 'Precision' )
         DO 100 I = 1, N
            DWORK(I) = DNRM2( MIN( I+1, N ), A(1,I), 1 )
            DWORK(N+I) = DNRM2( I, E(1,I), 1 )
  100    CONTINUE
         NORMA = DNRM2( N, DWORK, 1 )
         NORME = DNRM2( N, DWORK(N+1), 1 )
         IF ( ISDISC ) THEN
            FERR = ( NORMA**2 + NORME**2 )*EPS/SEP
         ELSE
            FERR = TWO*NORMA*NORME*EPS/SEP
         END IF
      END IF
!
      RETURN
! *** Last line of DGGLYP ***
      END
