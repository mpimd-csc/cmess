SUBROUTINE DTRLYP(TRANS,N,A,LDA, X, LDX, SCALE, INFO)
!     PURPOSE
!
!     To solve for X either the reduced generalized continuous-time
!     Lyapunov equation
!
!         T
!        A  * X  + X * A  =  SCALE * Y                       (1)
!
!     or
!
!                      T
!        A * X  + X * A   =  SCALE * Y                       (2)
!
!     where the right hand side Y is symmetric. A, Y, and the
!     solution X are N-by-N matrices. The pencil A - lambda * E must be
!     in generalized Schur form (A upper quasitriangular, E upper
!     triangular). SCALE is an output scale factor, set to avoid
!     overflow in X.
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
!     TRANS   CHARACTER*1
!             Specifies whether the transposed equation is to be solved
!             or not:
!             = 'N':  Solve equation (1);
!             = 'T':  Solve equation (2).
!
!     NB      INTEGER
!             Block size in the Bartels-Stewart algorithm. This value
!             must be set an optimal value for the used computer
!             architecture and the used BLAS library. A good value
!             if ISOLVE .EQ. 3 might be 32, 48 or 64 otherwise
!             16, 32 or 48 are good choices.
!
!     Input/Output Parameters
!
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!
!     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!             The leading N-by-N upper Hessenberg part of this array
!             must contain the quasitriangular matrix A.
!
!     LDA     INTEGER
!             The leading dimension of the array A.  LDA >= MAX(1,N).
!
!     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
!             On entry, the leading N-by-N part of this array must
!             contain the right hand side matrix Y of the equation. Only
!             the upper triangular part of this matrix need be given.
!             On exit, the leading N-by-N part of this array contains
!             the solution matrix X of the equation.
!
!     LDX     INTEGER
!             The leading dimension of the array X.  LDX >= MAX(1,N).
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
!     [2] Penzl, T.
!         Numerical solution of generalized Lyapunov equations.
!         Advances in Comp. Math., vol. 8, pp. 33-48, 1998.
!
!     [4] Köhler M., Saak J.
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
!     Lyapunov equation
!
!     *********************************
 IMPLICIT NONE

 ! Input Arguments
 CHARACTER      TRANS
 INTEGER        N,LDA,LDX,INFO
 DOUBLE PRECISION A(LDA,*),X(LDX,*)
 DOUBLE PRECISION SCALE


 ! Local Variables
 LOGICAL NOTRAN
 INTEGER I
 INTEGER KL, KH, LL, LH, KB, LB
 INTEGER IINFO
 DOUBLE PRECISION SCAL
 LOGICAL SYM

 ! Local Constants
 DOUBLE PRECISION ONE,ZERO,MONE,HALF
 PARAMETER (ONE = 1d0, ZERO = 0d0, MONE = -1d0, HALF=0.5d0)
 INTEGER IONE

 INTEGER NB
 PARAMETER (IONE = 1, NB = 64)

 DOUBLE PRECISION WORK(NB, NB)

 ! External Subroutines
 EXTERNAL DCOPY,DGEMM,DLACPY,DTRSYL, DLASET


 ! External Functions
 EXTERNAL LSAME
 LOGICAL LSAME

 ! Intrisics
 INTRINSIC MAX, ABS

 ! Check Input
 INFO = 0
 SYM = .TRUE.
 IF ( (.NOT. LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'T'))) THEN
         INFO = -1
 ELSE IF ( N < 0 ) THEN
         INFO = -2
 ELSE IF ( LDA < MAX(1,N) ) THEN
         INFO = -4
 ELSE IF ( LDX < MAX(1,N) ) THEN
         INFO = -6
 END IF
 IF ( INFO .NE. 0 ) THEN
         CALL XERBLA('DTRLYP',-INFO)
         RETURN
 END IF
 ! QUICK RETURN
 IF ( ( N .EQ. 0 ) ) THEN
         RETURN
 END IF



 NOTRAN = LSAME ( TRANS, 'N')
 SCALE = ONE
 SCAL = ONE


 ! Solve the non transposed case
 IF (NOTRAN) THEN
     KL = 0
     KH = 0
     DO WHILE (KH .NE. N)
           KL = KH + 1
           KH = KH + NB
           IF ( KH < N .AND. A(KH+1,KH) .NE. ZERO ) THEN
               KH = KH + 1
           END IF
           IF ( KH > N ) THEN
               KH = N
           END IF
           KB = KH - KL + 1

           ! Copy Symmetry
           IF ( KL > 1 ) THEN
                DO I = KL,KH
                    CALL DCOPY( KL-1, X(1,I), 1, X(I,1), LDX )
                END DO
           END IF


           LL = 0
           LH = KL-1

           DO WHILE ( LH .NE. N )
                LL = LH + 1
                LH = LH + NB
                IF ( LH < N .AND. A(LH+1,LH) .NE. ZERO ) THEN
                    LH = LH + 1
                END IF
                IF ( LH > N ) THEN
                    LH = N
                END IF
                LB = LH - LL + 1

                ! Update RHS
                IF ( LH > 1 ) THEN
                    CALL DGEMM( 'N', 'N', KB, LB, LL-1, MONE, X(KL,1), LDX, A(1,LL), LDA, ONE, X(KL,LL), LDX )
                END IF


                ! Solve inner System
                CALL DTRSYL('T', 'N', IONE, KB, LB, A(KL,KL), LDA, A(LL,LL), LDA, X(KL,LL), LDX, SCAL, IINFO)
                IF (SCAL .NE. ONE ) THEN
                    CALL DLACPY("All",KB,LB,X(KL,LL),LDX,WORK,KB)
                    DO I = 1,N
                        CALL DSCAL(N,SCAL,X(1,I),IONE)
                    END DO
                    CALL DLACPY("All",KB,KL,WORK,KB,X(KL,LL),LDX)
                    SCALE = SCALE * SCAL
                END IF

                ! Fix Symmetry of the Diagonal block
                IF ( LL .EQ. KL ) THEN
                    DO I = KL,KH
                        CALL DSCAL(KH-I,HALF,X(I+1,I),IONE)
                        CALL DAXPY(KH-I,HALF,X(I,I+1),LDX,X(I+1,I),IONE)
                        CALL DCOPY(KH-I,X(I+1,I),IONE,X(I,I+1),LDX)
                    END DO
                END IF

                ! Update RHS
                IF ( KL < LL ) THEN
                    CALL DGEMM("T","N",LH-KH,LB,KB,MONE,A(KL,KH+1),LDA,X(KL,LL), LDX, ONE,X(KH+1,LL), LDX)
                END IF
           END DO
   END DO
 ELSE
     ! CALL DLASET("All", N,N, ZERO, ONE, E, LDE)
 ! Solve the transposed case
   LH = N+1
   LL = N+1
   DO WHILE ( LL .NE. 1 )
        LH = LL - 1
        LL = LH - NB +1
        IF ( LL .GT. 1 .AND. A(LL, LL-1) .NE. ZERO) THEN
            LL = LL -1
        END IF
        IF ( LL .LT. 1 ) THEN
            LL = 1
        END IF
        LB = LH - LL + 1

        IF ( LH < N ) THEN
                DO I = LL, LH
                        CALL DCOPY( N-LH, X(I,LH+1), LDX, X(LH+1,I), 1 )
                END DO
        END IF

        KL = LH +1
        KH = LH +1
        DO WHILE ( KL .NE. 1)
           KH = KL - 1
           KL = KH - NB +1
           IF ( KL .GT. 1 .AND. A(KL, KL-1) .NE. ZERO) THEN
               KL = KL -1
           END IF

           IF ( KL .LT. 1 ) THEN
               KL = 1
           END IF
           KB = KH - KL +1

           IF ( KH < N ) THEN
                  CALL DGEMM( 'N', 'N', KB, LB, N-KH, MONE, A(KL,KH+1), LDA, X(KH+1,LL), LDX, ONE, X(KL,LL), LDX)
           END IF

           CALL DTRSYL('N', 'T', IONE, KB, LB, A(KL,KL), LDA, A(LL,LL), LDA, X(KL,LL), LDX, SCAL, IINFO)

           IF (SCAL .NE. ONE ) THEN
               CALL DLACPY("All",KB,LB,X(KL,LL),LDX,WORK,KB)
               DO I = 1,N
                   CALL DSCAL(N,SCAL,X(1,I),IONE)
               END DO
               CALL DLACPY("All",KB,KL,WORK,KB,X(KL,LL),LDX)
               SCALE = SCALE * SCAL
           END IF

           ! Fix Symmetry of the Diagonal block
           IF ( LL .EQ. KL ) THEN
               DO I = KL,KH
                   CALL DSCAL(KH-I,HALF,X(I+1,I),IONE)
                   CALL DAXPY(KH-I,HALF,X(I,I+1),LDX,X(I+1,I),IONE)
                   CALL DCOPY(KH-I,X(I+1,I),IONE,X(I,I+1),LDX)
               END DO
           END IF


           IF ( KH < LH ) THEN
               CALL DGEMM( 'N', 'T', KB, LL-KL, LB, MONE, X(KL,LL), LDX,A(KL,LL), LDA, ONE, X(KL,KL), LDX )
           END IF
        END DO
   END DO
 END IF
END SUBROUTINE

