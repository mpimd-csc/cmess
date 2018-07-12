SUBROUTINE DTGSYY(TRANS,M,N,A,LDA,B,LDB,C,LDC,D,LDD,X,LDX,WORK,SCALE,INFO)
!
!     PURPOSE
!
!     To solve for X either the reduced generalized Sylvester equation
!
!         T            T
!        A  * X * B + C  * X * D  =  SCALE * Y                       (1)
!
!     or
!
!                 T            T
!        A * X * B  + C * X * D   =  SCALE * Y                       (2)
!
!     where (A,C) and (B,D) are in generalized Schur form. The matrices
!     A and C or respectively B and D are transposed implicitly inside
!     the algorithm. Furthermore the matrices A and B must be upper
!     quasitriangular and the matrices C and D must be upper triangular.
!     SCALE is an output scale factor, set to avoid
!     overflow in X.
!
!     In contrast to DTGSYY2 this variant uses a modification of the
!     Generalized Schur(GS) algorithm by Kagström and Westin. If the
!     matrices A and D are upper quasitriangular  and  the matrices
!     B and C are upper triangular we have to use DTGSYX or DTGSYX2
!     instead.
!
!     This is a level 2 BLAS implementation.
!
!     ARGUMENTS
!
!     Mode Parameters
!
!     TRANS   CHARACTER*1
!             Specifies whether the transposed equation is to be solved
!             or not:
!             = 'N':  Solve equation (1);
!             = 'T':  Solve equation (2).
!
!
!     Input/Output Parameters
!
!     M       (input) INTEGER
!             The order of the matrix A and C.  M >= 0.
!
!     N       (input) INTEGER
!             The order of the matrix B and D.  N >= 0.
!
!     A       (input) DOUBLE PRECISION array, dimension (LDA,M)
!             The leading M-by-M upper Hessenberg part of this array
!             must contain the quasitriangular matrix A.
!
!     LDA     INTEGER
!             The leading dimension of the array A.  LDA >= MAX(1,M).
!
!     B       (input) DOUBLE PRECISION array, dimension (LDB,N)
!             The leading N-by-N upper triangular part of this array
!             must contain the matrix B.
!
!     LDB     INTEGER
!             The leading dimension of the array B.  LDB >= MAX(1,N).
!
!     C       (input) DOUBLE PRECISION array, dimension (LDC,M)
!             The leading M-by-M upper triangular part of this array
!             must contain the matrix C.
!
!     LDC     INTEGER
!             The leading dimension of the array C.  LDB >= MAX(1,M).
!
!     D       (input) DOUBLE PRECISION array, dimension (LDD,N)
!             The leading N-by-N upper Hessenberg part of this array
!             must contain the quasitriangular matrix D.
!
!     LDD     INTEGER
!             The leading dimension of the array A.  LDD >= MAX(1,N).
!
!     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
!             On entry, the leading M-by-N part of this array must
!             contain the right hand side matrix Y of the equation.
!             On exit, the leading M-by-N part of this array contains
!             the solution matrix X of the equation.
!
!     LDX     INTEGER
!             The leading dimension of the array X.  LDX >= MAX(1,M).
!
!     WORK    DOUBLE PRECISION array, dimension (M*N)
!             Double Precision workspace.
!
!     SCALE   (output) DOUBLE PRECISION
!             The scale factor set to avoid overflow in X.
!             (0 < SCALE <= 1)
!             It is only computed iff ISOLVE .EQ. 3 otherwise it is set to 1.
!
!     Error indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  equation is (almost) singular to working precision;
!                   perturbed values were used to solve the equation
!                   (but the matrices A and E are unchanged).
!
!     METHOD
!
!     The solution X of (1) or (2) is computed via block back
!     substitution or block forward substitution, respectively. (See
!     [1] and [2] for details.)
!
!     REFERENCES
!
!     [1] Bartels, R.H., Stewart, G.W.
!         Solution of the equation A X + X B = C.
!         Comm. A.C.M., 15, pp. 820-826, 1972.
!
!     [2] Köhler M., Saak J.
!         A level 3 block variant of the Bartel-Stewart Algorithm for
!         the generalized Lyapunov equation.
!         to be written
!
!
!     Author
!
!      Köhler M. , MPI Magdeburg
!
!     KEYWORDS
!
!     Generalized Sylvester equation
!
!     *********************************
IMPLICIT NONE
! Input Parameters
CHARACTER TRANS
INTEGER  M,N,LDA,LDB,LDC,LDD,LDX
INTEGER  INFO
DOUBLE PRECISION A(LDA,*)
DOUBLE PRECISION B(LDB,*)
DOUBLE PRECISION C(LDC,*)
DOUBLE PRECISION D(LDD,*)
DOUBLE PRECISION X(LDX,*)
DOUBLE PRECISION WORK(M,*)
DOUBLE PRECISION SCALE

! Local Variables
INTEGER I,J,IS,IE,JS,JE,IB,JB
INTEGER DIMAT
DOUBLE PRECISION G(8,8)
DOUBLE PRECISION RHS(8)
INTEGER RPIV(8), CPIV(8)
DOUBLE PRECISION SCAL
INTEGER IINFO,TEQN
DOUBLE PRECISION ZERO, MONE, ONE

EXTERNAL DLASET,DGEMM,DGETC2,DGESC2,DTRSM,DLACPY
EXTERNAL LSAME
LOGICAL LSAME

! Check Input
INFO = 0
IF ( (.NOT. LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'T'))) THEN
    INFO = -1
ELSE IF ( M < 0) THEN
    INFO = -2
ELSE IF ( N < 0) THEN
    INFO = -3
ELSE IF ( LDA < M ) THEN
    INFO = -5
ELSE IF ( LDB < N ) THEN
    INFO = -7
ELSE IF ( LDC < M ) THEN
    INFO = -9
ELSE IF ( LDD < N ) THEN
    INFO = -11
ELSE IF ( LDX < M ) THEN
    INFO = -13
END IF
IF ( INFO .NE. 0 ) THEN
    CALL XERBLA('DTGSYX',-INFO)
    RETURN
END IF
! QUICK RETURN
IF ( ( N == 0 ) .OR. (M==0)) THEN
    RETURN
END IF

IF ( LSAME(TRANS,'N') ) THEN
    TEQN = 1
ELSE
    TEQN = 0
END IF

ZERO = 0d0
ONE = 1d0
MONE = -1d0
SCALE = 1d0

CALL DLASET('All',M,N,ZERO,ZERO,WORK,M)

! Start the algorithm
IF ( TEQN .EQ. 1 ) THEN
    I = 1
    DO WHILE ( I <= M )
        ! Check for BLOCK in A
        IF ( I < M .AND. (A(I+1,I) .NE. 0) ) THEN
            IS = I
            IE = I + 1
            IB = 2
            I = I +2
        ELSE
            IS = I
            IE = I
            IB = 1
            I = I +1
        END IF

        J = 1
        DO WHILE ( J <= N )
            ! CHECK for BLOCK in D
            IF ( J < N .AND. (B(J+1,J) .NE. 0)) THEN
                JS = J
                JE = J+1
                JB = 2
                J = J + 2
            ELSE
                JS = J
                JE = J
                JB = 1
                J = J +1
            END IF
            ! Solve the small inner system
            ! Assemble the Matrix
            IF ( IB == 1 .AND. JB == 1 ) THEN
                G(1,1) =  A (IS,IS)
                G(1,2) =  -D (JS,JS)
                G(2,1) =  C(IS,IS)
                G(2,2) = -B(JS,JS)
                RHS(1) = X(IS,JS)
                RHS(2) = WORK(IS,JS)

                DIMAT = 2

            ELSE IF ( IB == 1 .AND. JB == 2) THEN
                CALL DLASET('All',4,4,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,3) = -D(JS,JS)
                G(2,2) = A(IS,IS)
                G(2,3) = -D(JS,JS+1)
                G(2,4) = -D(JS+1,JS+1)
                G(3,1) = C(IS, IS)
                G(3,3) = -B(JS,JS)
                G(3,4) = -B(JS+1,JS)
                G(4,2) = C(IS,IS)
                G(4,3) = -B(JS,JS+1)
                G(4,4) = -B(JS+1,JS+1)

                RHS(1) = X(IS,JS)
                RHS(2) = X(IS,JS+1)
                RHS(3) = WORK ( IS, JS)
                RHS(4) = WORK ( IS, JS+1)

                DIMAT  =  4

            ELSE IF ( IB == 2 .AND. JB == 1) THEN
                CALL DLASET('All',4,4,ZERO,ZERO,G,8)
                G(1,1) = A(IS, IS )
                G(1,2) = A(IS+1, IS)
                G(1,3) = -D(JS, JS)
                G(2,1) = A(IS, IS+1)
                G(2,2) = A(IS+1,IS+1)
                G(2,4) = -D(JS,JS)
                G(3,1) = C(IS, IS)
                G(3,3) = -B(JS,JS)
                G(4,1) = C(IS,IS+1)
                G(4,2) = C(IS+1,IS+1)
                G(4,4) = -B(JS,JS)

                RHS(1) = X (IS,JS)
                RHS(2) = X (IS+1, JS)
                RHS(3) = WORK ( IS, JS)
                RHS(4) = WORK ( IS+1, JS)

                DIMAT = 4
            ELSE IF ( IB == 2 .AND. JB == 2) THEN
                CALL DLASET('All',8,8,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,2) = A(IS+1,IS)
                G(1,5) = -D(JS,JS)

                G(2,1) = A(IS,IS+1)
                G(2,2) = A(IS+1,IS+1)
                G(2,6) = -D(JS,JS)

                G(3,3) = A(IS,IS)
                G(3,4) = A(IS+1,IS)
                G(3,5) = -D(JS,JS+1)
                G(3,7) = -D(JS+1,JS+1)

                G(4,3) = A(IS,IS+1)
                G(4,4) = A(IS+1,IS+1)
                G(4,6) = -D(JS,JS+1)
                G(4,8) = -D(JS+1,JS+1)

                G(5,1) = C(IS,IS)
                G(5,5) = -B ( JS, JS)
                G(5,7) = -B ( JS+1, JS)

                G(6,1) = C(IS,IS+1)
                G(6,2) = C(IS+1,IS+1)
                G(6,6) = -B(JS, JS )
                G(6,8) = -B(JS+1,JS)

                G(7,3) = C(IS,IS)
                G(7,5) = -B(JS,JS+1)
                G(7,7) = -B(JS+1,JS+1)

                G(8,3) = C(IS,IS+1)
                G(8,4) = C(IS+1,IS+1)
                G(8,6) = -B(JS, JS+1)
                G(8,8) = -B(JS+1,JS+1)

                RHS(1) = X(IS,JS)
                RHS(2) = X(IS+1,JS)
                RHS(3) = X(IS, JS+1)
                RHS(4) = X(IS+1,JS+1)

                RHS(5) = WORK(IS,JS)
                RHS(6) = WORK(IS+1,JS)
                RHS(7) = WORK(IS, JS+1)
                RHS(8) = WORK(IS+1,JS+1)

                DIMAT = 8
            END IF

            ! SOLVE
            CALL DGETC2(DIMAT,G,8,RPIV,CPIV,IINFO)
            IF (IINFO .NE. 0 ) THEN
                INFO = 4
                RETURN
            END IF
            SCAL = ONE
            CALL DGESC2(DIMAT,G,8,RHS,RPIV,CPIV,SCAL)

            ! Scale if necessary
            IF ( SCAL .NE. ONE ) THEN
                DO I = 1, N
                    CALL DSCAL( N, SCAL, X(1,I), 1 )
                    CALL DSCAL( N, SCAL, WORK(1,I),1)
                END DO
                SCALE = SCALE*SCAL
            END IF

            ! Copy the solution back
            IF ( IB == 1 .AND. JB == 1 ) THEN
                X(IS,JS) = RHS(1)
                WORK(IS,JS) = RHS(2)
            ELSE IF ( IB == 1 .AND. JB == 2) THEN
                X(IS,JS)  = RHS(1)
                X(IS,JS+1)  = RHS(2)
                WORK ( IS, JS)  = RHS(3)
                WORK ( IS, JS+1)  = RHS(4)
            ELSE IF ( IB == 2 .AND. JB == 1) THEN
                X (IS,JS) = RHS(1)
                X (IS+1, JS)  = RHS(2)
                WORK ( IS, JS)  = RHS(3)
                WORK ( IS+1, JS)  = RHS(4)
            ELSE IF ( IB == 2 .AND. JB == 2) THEN
                X(IS,JS) = RHS(1)
                X(IS+1,JS) = RHS(2)
                X(IS, JS+1)  = RHS(3)
                X(IS+1,JS+1)  = RHS(4)

                WORK(IS,JS) = RHS(5)
                WORK(IS+1,JS) = RHS(6)
                WORK(IS, JS+1)  = RHS(7)
                WORK(IS+1,JS+1)  = RHS(8)
            END IF

            ! Update the Right hand sides
            IF (IE < M ) THEN
                CALL DGEMM('T','N', M-IE, JB, IB, MONE, A(IS,IE+1), LDA , X(IS,JS), LDX,  &
                    & ONE, X(IE+1,JS),LDX)
                CALL DGEMM('T','N', M-IE, JB, IB, MONE, C(IS,IE+1), LDC , X(IS,JS), LDX,  &
                    & ONE, WORK(IE+1,JS),M)
            END IF

            IF (JE < N ) THEN
                CALL DGEMM('N','N', IB, N-JE, JB, ONE, WORK(IS,JS), M , D(JS,JE+1), LDD, &
                    & ONE, X(IS, JE+1), LDX)
                CALL DGEMM('N','N', IB, N-JE, JB, ONE, WORK(IS,JS), M , B(JS,JE+1), LDB, &
                    & ONE, WORK(IS, JE+1), M)
            END IF
        END DO
    END DO

    CALL DLACPY('A', M,N,WORK,M,X,LDX)
    CALL DTRSM('L','U','T','N',M,N,ONE,C,LDC,X,LDX)
ELSE
    ! SOLVE THE TRANSPOSE EQN
    I = M
    DO WHILE ( I >= 1 )
        IF ( I > 1 .AND. (A(I,I-1) .NE. 0.0d0)) THEN
            IS = I-1
            IE = I
            I = I -2
            IB = 2
        ELSE
            IS = I
            IE = I
            I = I - 1
            IB = 1
        END IF
        J = N

        DO WHILE (J >=1)
            IF ( J > 1 .AND. (B(J,J-1) .NE. 0.0d0)) THEN
                JS = J -1
                JE = J
                J = J - 2
                JB = 2
            ELSE
                JS = J
                JE = J
                J = J - 1
                JB = 1
            END IF

            IF ( IB == 1 .AND. JB == 1 ) THEN
                CALL DLASET('All',2,2,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,2) = -D(JS,JS)
                G(2,1) = C(IS,IS)
                G(2,2) = -B(JS,JS)
                RHS(1) = X(IS,JS)
                RHS(2) = WORK(IS,JS)
                DIMAT = 2
            ELSE IF (IB == 1 .AND. JB == 2) THEN
                CALL DLASET('All',4,4,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,3) = -D(JS,JS)
                G(2,2) = A(IS,IS)
                G(2,3) = -D(JS+1,JS)
                G(2,4) = -D(JS+1, JS+1)

                G(3,1) = C(IS,IS)
                G(3,3) = -B(JS,JS)
                G(3,4) = -B(JS,JS+1)
                G(4,2) = C(IS,IS)
                G(4,3) = -B(JS+1, JS)
                G(4,4) = -B(JS+1, JS+1)

                RHS(1) = X(IS,JS)
                RHS(2) = X(IS,JS+1)
                RHS(3) = WORK ( IS, JS)
                RHS(4) = WORK ( IS, JS+1)
                DIMAT = 4

            ELSE IF (IB == 2 .AND. JB == 1) THEN
                CALL DLASET('All',4,4,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,2) = A(IS,IS+1)
                G(1,3) = -D(JS,JS)
                G(2,1) = A(IS+1,IS)
                G(2,2) = A(IS+1,IS+1)
                G(2,4) = -D(JS,JS)
                G(3,1) = C(IS,IS)
                G(3,2) = C(IS,IS+1)
                G(3,3) = -B(JS,JS)
                G(4,2) = C(IS+1,IS+1)
                G(4,4) = -B(JS,JS)

                RHS(1) = X (IS,JS)
                RHS(2) = X (IS+1, JS)
                RHS(3) = WORK ( IS, JS)
                RHS(4) = WORK ( IS+1, JS)
                DIMAT = 4
            ELSE
                CALL DLASET('All',8,8,ZERO,ZERO,G,8)
                G(1,1) = A(IS,IS)
                G(1,2) = A(IS,IS+1)
                G(1,5) = -D(JS,JS)
                G(1,7) = -D(JS,JS+1)

                G(2,1) = A(IS+1,IS)
                G(2,2) = A(IS+1,IS+1)
                G(2,6) = -D(JS,JS)
                G(2,8) = -D(JS,JS+1)

                G(3,3) = A(IS,IS)
                G(3,4) = A(IS,IS+1)
                G(3,7) = -D(JS+1,JS+1)

                G(4,3) = A(IS+1,IS)
                G(4,4) = A(IS+1,IS+1)
                G(4,8) = -D(JS+1,JS+1)

                G(5,1) = C(IS,IS)
                G(5,2) = C(IS,IS+1)
                G(5,5) = -B ( JS, JS)
                G(5,7) = -B ( JS, JS+1)

                G(6,2) = C(IS+1,IS+1)
                G(6,6) = -B(JS,JS)
                G(6,8) = -B(JS, JS+1)

                G(7,3) = C(IS,IS)
                G(7,4) = C(IS,IS+1)
                G(7,5) = -B(JS+1,JS)
                G(7,7) = -B(JS+1,JS+1)

                G(8,4) = C(IS+1,IS+1)
                G(8,6) = -B(JS+1, JS)
                G(8,8) = -B(JS+1,JS+1)

                RHS(1) = X(IS,JS)
                RHS(2) = X(IS+1,JS)
                RHS(3) = X(IS, JS+1)
                RHS(4) = X(IS+1,JS+1)

                RHS(5) = WORK(IS,JS)
                RHS(6) = WORK(IS+1,JS)
                RHS(7) = WORK(IS, JS+1)
                RHS(8) = WORK(IS+1,JS+1)
                DIMAT = 8
            ENDIF

            CALL DGETC2(DIMAT,G,8,RPIV,CPIV,IINFO)
            IF (IINFO .NE. 0 ) THEN
                INFO = 4
                RETURN
            END IF
            SCAL = ONE
            CALL DGESC2(DIMAT,G,8,RHS,RPIV,CPIV,SCAL)

            ! Scale if necessary
            IF ( SCAL .NE. ONE ) THEN
                DO I = 1, N
                    CALL DSCAL( N, SCAL, X(1,I), 1 )
                    CALL DSCAL( N, SCAL, WORK(1,I),1)
                END DO
                SCALE = SCALE*SCAL
            END IF

            IF ( IB == 1 .AND. JB == 1 ) THEN
                X(IS,JS) = RHS(1)
                WORK(IS,JS) = RHS(2)
            ELSE IF (IB == 1 .AND. JB == 2) THEN
                X(IS,JS)  = RHS(1)
                X(IS,JS+1)  = RHS(2)
                WORK ( IS, JS)  = RHS(3)
                WORK ( IS, JS+1)  = RHS(4)
            ELSE IF (IB == 2 .AND. JB == 1) THEN
                X (IS,JS) = RHS(1)
                X (IS+1, JS)  = RHS(2)
                WORK ( IS, JS)  = RHS(3)
                WORK ( IS+1, JS)  = RHS(4)
            ELSE
                X(IS,JS) = RHS(1)
                X(IS+1,JS) = RHS(2)
                X(IS, JS+1)  = RHS(3)
                X(IS+1,JS+1)  = RHS(4)

                WORK(IS,JS) = RHS(5)
                WORK(IS+1,JS) = RHS(6)
                WORK(IS, JS+1)  = RHS(7)
                WORK(IS+1,JS+1)  = RHS(8)

            ENDIF

            IF ( IS > 1 ) THEN
                CALL DGEMM('N','N',IS-1,JB,IB,MONE,A(1,IS),LDA,X(IS,JS),LDX,ONE,X(1,JS),LDX)
                CALL DGEMM('N','N',IS-1,JB,IB,MONE,C(1,IS),LDC,X(IS,JS),LDX,ONE,WORK(1,JS),M)
            END IF

            IF ( JS >  1 ) THEN
                CALL DGEMM('N','T',IB,JS-1,JB,ONE,WORK(IS,JS),M,D(1,JS),LDD,ONE, X(IS,1),LDX)
                CALL DGEMM('N','T',IB,JS-1,JB,ONE,WORK(IS,JS),M,B(1,JS),LDB,ONE, WORK(IS,1),M)
            END IF
        END DO
    END DO

    CALL DLACPY('A', M,N,WORK,M,X,LDX)
    CALL DTRSM('L','U','N','N',M,N,ONE,C,LDC,X,LDX)
END IF
END SUBROUTINE


