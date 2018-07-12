SUBROUTINE DTGSYX2(TRANS,M,N,A,LDA,B,LDB,C,LDC,D,LDD,X,LDX,WORK,IWORK,INFO)
!
!     PURPOSE
!
!     To solve for X either the reduced generalized Sylvester equation
!
!         T            T
!        A  * X * B + C  * X * D  =  Y                       (1)
!
!     or
!
!                 T            T
!        A * X * B  + C * X * D   =  Y                       (2)
!
!     where (A,C) and (D,B) are in generalized Schur form. The matrices
!     A and C or respectively B and D are transposed implicitly inside
!     the algorithm. Furthermore the matrices A and D must be upper
!     quasitriangular and the matrices B and C must be upper triangular.
!
!     In contrast to DTGSYX this variant uses a modification of the
!     Gardiner and Laub method for the generalized Sylvester equation.
!     Due the internal design of the algorithm this variant DOES NOT
!     provide a scaling factor for the right hand side Y as it is done
!     by DTGSYX.
!     If the matrices A and B are upper quasitriangular  and  the matrices
!     C and D are upper triangular we have to use DTGSYY or DTGSYY2
!     instead.
!
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
!     WORK    DOUBLE PRECISION array, dimension (4*M*(M+1)+4*(M+1))
!             Double Precision workspace.
!
!     IWORK   INTEGER array, dimension (2*M)
!             Integer workspace
!
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
!     The solution X of (1) or (2) is computed via a modified Bartels-Stewart
!     algorithm. Depending on the ISOLVE parameter the internal solver is
!     only a PLU decompositions or a forward/backward subsitution scheme
!     (See  [1] and [2] for details.)
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
! Arguments
CHARACTER TRANS
INTEGER  M,N,LDA,LDB,LDC,LDD,LDX
INTEGER  INFO
DOUBLE PRECISION A(LDA,*)
DOUBLE PRECISION B(LDB,*)
DOUBLE PRECISION C(LDC,*)
DOUBLE PRECISION D(LDD,*)
DOUBLE PRECISION X(LDX,*)
DOUBLE PRECISION WORK(2*M,2*M+1)
INTEGER IWORK(2*M)


! Local Variables
DOUBLE PRECISION ZERO, MONE, ONE
PARAMETER(ONE = 1d0, ZERO = 0d0, MONE = -1d0)
DOUBLE PRECISION AA,BB,CC,SS
INTEGER IINFO,TEQN,IONE,ITWO
PARAMETER(IONE = 1, ITWO = 2)
INTEGER K,L,NGIVENS, I
INTEGER IPIV(4), JPIV(4)

! Additional Functions
EXTERNAL DLASET,DGEMM,DGETC2,DLACPY,DGESC2X
EXTERNAL DGEMV,DCOPY
EXTERNAL DROT, DROTG,DTRSV,DSCAL,DAXPY
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
    CALL XERBLA('DTGSYX2',-INFO)
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

IINFO = 0

IF (TEQN .EQ. 1 ) THEN
    ! SOLVE A^TXB+C^TXD = X
    K = 1
    DO WHILE ( K <= N)
        IF (((K .LT. N) .AND. (D(K+1,K) .EQ. ZERO)) .OR. (K .EQ. N)) THEN
            ! X(:,k)=(B(k,k)*A+D(k,k)*C)'\X(:,k);
            ! WORK(1:2*M,1:2*M) = ZERO
            DO I = 1, M
                CALL DCOPY(M, A(1,I), IONE, WORK(1,I), IONE)
                CALL DSCAL(M, B(K,K), WORK(1,I), IONE)
                CALL DAXPY(M, D(K,K), C(1,I), IONE, WORK(1,I), IONE)
            END DO

            L = 1
            NGIVENS = 0
            DO WHILE ( L < M )
                IF ( WORK(L+1,L) .NE. ZERO ) THEN
                    NGIVENS = NGIVENS + 1
                    IWORK(NGIVENS) = L
                    AA = WORK(L,L)
                    BB = WORK(L+1,L)
                    CALL DROTG(AA,BB,WORK(2*NGIVENS-1,2*M+1),WORK(2*NGIVENS,2*M+1))
                    CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,WORK(2*NGIVENS-1,2*M+1),WORK(2*NGIVENS,2*M+1))
                    L = L +2
                ELSE
                    L = L +1
                END IF
            END DO
            CALL DTRSV("U","T","N",M,WORK,2*M,X(1,K),IONE)

            DO L = 1,NGIVENS
                CALL DROT(IONE,X(IWORK(L),K),IONE,X(IWORK(L)+1,K),IONE,WORK(2*L-1,2*M+1),-WORK(2*L,2*M+1))
            END DO

            !Update Rest
            !Precompute A^T X_k and C^T X_k
            CALL DGEMV("T",M,M,ONE, A, LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
            CALL DGEMV("T",M,M,ONE, C, LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)
            DO L = K+1,N
                CALL DAXPY(M,-B(K,L),WORK(1,1),IONE,X(1,L),IONE)
                CALL DAXPY(M,-D(K,L),WORK(1,2),IONE,X(1,L),IONE)
            END DO
            K = K + 1
        ELSE
            WORK(1:2*M, 1:2*M+1) = ZERO
            ! % Transponieren der matrizen entfernen dafür transponiert lösen
            ! TMP = [ (B(k,k)*A+D(k,k)*C) ,(B(k,k+1)*A+D(k,k+1)*C) ;
            !         D(k+1,k)*C, (B(k+1,k+1)*A+D(k+1,k+1)*C)]'\reshape(X(:,k:k+1),2*m,1);
            DO L = 1, M
                CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                CALL DSCAL(M,B(K,K),WORK(1,2*L-1),ITWO)
                CALL DAXPY(M,D(K,K),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                CALL DSCAL(M,B(K,K+1),WORK(1,2*L),ITWO)
                CALL DAXPY(M,D(K,K+1),C(1,L),IONE,WORK(1,2*L),ITWO)

                CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                CALL DSCAL(M,D(K+1,K),WORK(2,2*L-1),ITWO)

                CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                CALL DSCAL(M,B(K+1,K+1),WORK(2,2*L),ITWO)
                CALL DAXPY(M,D(K+1,K+1),C(1,L),IONE,WORK(2,2*L),ITWO)

            END DO

            CALL DCOPY(M, X(1,K), IONE, WORK(1,2*M+1),ITWO)
            CALL DCOPY(M, X(1,K+1),IONE, WORK(2,2*M+1),ITWO)

            L = 1
            DO WHILE ( L <= 2*M  )
                IF ( (L+3 <= 2*M) .AND. ( WORK(L+2,L) .NE. ZERO) ) THEN
                    CALL DGETC2(4, WORK(L,L), 2*M , IPIV, JPIV, IINFO)
                    CALL DGESC2X("T",4, WORK(L,L), 2* M , WORK(L,2*M+1), IPIV, JPIV,AA)
                    CALL DGEMV("T",4, 2*M-L-3, MONE, WORK(L,L+4),2*M,WORK(L,2*M+1), &
                        & IONE,ONE,WORK(L+4,2*M+1),IONE)
                    L = L + 4
                ELSE
                    CALL DGETC2(2, WORK(L,L), 2*M , IPIV, JPIV, IINFO)
                    CALL DGESC2X("T",2, WORK(L,L), 2* M , WORK(L,2*M+1), IPIV, JPIV,AA)

                    CALL DGEMV("T",2, 2*M-L-1, MONE, WORK(L,L+2),2*M,WORK(L,2*M+1), &
                        & IONE,ONE,WORK(L+2,2*M+1),IONE)
                    L = L + 2
                END IF
            END DO

            CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K), IONE)
            CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K+1),IONE)

            ! Precompute the Update
            CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
            CALL DGEMV("T",M,M,ONE,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
            CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K+1),IONE, ZERO, WORK(1,2),IONE)
            CALL DGEMV("T",M,M,ONE,C,LDA, X(1,K+1),IONE, ZERO, WORK(M+1,2),IONE)

            DO L = K+2,N
                CALL DAXPY(M,-B(K,L),WORK(1,1),IONE, X(1,L),IONE)
                CALL DAXPY(M,-D(K,L),WORK(M+1,1),IONE, X(1,L),IONE)

                CALL DAXPY(M,-B(K+1,L),WORK(1,2),IONE, X(1,L),IONE)
                CALL DAXPY(M,-D(K+1,L),WORK(M+1,2),IONE, X(1,L),IONE)
            END DO
            K = K + 2
        END IF
    END DO
ELSE
    K = N
    DO WHILE ( K >= 1 )
        IF ( K .EQ. 1 .OR. (K > 1 .AND. D(K,K-1) .EQ. ZERO)) THEN
            ! WORK(1:M,1:M) = ZERO
            DO I = 1, M
                CALL DCOPY(M, A(1,I), IONE, WORK(1,I), IONE)
                CALL DSCAL(M, B(K,K), WORK(1,I), IONE)
                CALL DAXPY(M, D(K,K), C(1,I), IONE, WORK(1,I),IONE)
            END DO

            ! New one using givens rotations
            L = 1
            DO WHILE ( L < M )
                IF ( WORK(L+1,L) .NE. ZERO ) THEN
                    AA = WORK(L,L)
                    BB = WORK(L+1,L)
                    CALL DROTG(AA,BB,CC,SS)
                    CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,CC,SS)
                    CALL DROT(IONE,X(L,K),IONE,X(L+1,K),IONE,CC,SS)
                    L = L +2
                ELSE
                    L = L +1
                END IF
            END DO
            CALL DTRSV("U","N","N",M,WORK,2*M,X(1,K),IONE)

            !Update Rest
            CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
            CALL DGEMV("N",M,M,ONE,C,LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)

            DO L = 1,K-1
                CALL DAXPY(M,-B(L,K),WORK(1,1),IONE,X(1,L),IONE)
                CALL DAXPY(M,-D(L,K),WORK(1,2),IONE,X(1,L),IONE)
            END DO
            K = K - 1
        ELSE
            ! WORK(1:2*M, 1:2*M+1) = ZERO

            ! New Idea
            DO L = 1, M
                CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                CALL DSCAL(M,B(K-1,K-1),WORK(1,2*L-1),ITWO)
                CALL DAXPY(M,D(K-1,K-1),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                CALL DSCAL(M,B(K-1,K),WORK(1,2*L),ITWO)
                CALL DAXPY(M,D(K-1,K),C(1,L),IONE,WORK(1,2*L),ITWO)

                CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                CALL DSCAL(M,D(K,K-1),WORK(2,2*L-1),ITWO)

                CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                CALL DSCAL(M,B(K,K),WORK(2,2*L),ITWO)
                CALL DAXPY(M,D(K,K),C(1,L),IONE,WORK(2,2*L),ITWO)

            END DO

            CALL DCOPY(M, X(1,K-1), IONE, WORK(1,2*M+1),ITWO)
            CALL DCOPY(M, X(1,K),IONE, WORK(2,2*M+1),ITWO)

            L = 2*M
            DO WHILE ( L >= 1  )
                IF ( (L >= 4) .AND. ( WORK(L,L-2) .NE. ZERO) ) THEN
                    CALL DGETC2(4, WORK(L-3,L-3), 2*M , IPIV, JPIV, IINFO)
                    CALL DGESC2X("N",4, WORK(L-3,L-3), 2* M , WORK(L-3,2*M+1), IPIV, JPIV,AA)

                    CALL DGEMV("N",L-4, 4, MONE, WORK(1,L-3),2*M,WORK(L-3,2*M+1), &
                        & IONE,ONE,WORK(1,2*M+1),IONE)
                    L = L - 4
                ELSE
                    CALL DGETC2(2, WORK(L-1,L-1), 2*M , IPIV, JPIV, IINFO)
                    CALL DGESC2X("N",2, WORK(L-1,L-1), 2* M , WORK(L-1,2*M+1), IPIV, JPIV,AA)
                    CALL DGEMV("N",L-2,2, MONE, WORK(1,L-1),2*M,WORK(L-1,2*M+1), &
                        & IONE,ONE,WORK(1,2*M+1),IONE)
                    L = L - 2
                END IF
            END DO

            CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K-1), IONE)
            CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K),IONE)


            CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
            CALL DGEMV("N",M,M,ONE,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
            CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K-1),IONE, ZERO, WORK(1,2),IONE)
            CALL DGEMV("N",M,M,ONE,C,LDA, X(1,K-1),IONE, ZERO, WORK(M+1,2),IONE)

            DO L = 1,K-2
                CALL DAXPY(M,-B(L,K),  WORK(1,1),IONE,X(1,L),IONE)
                CALL DAXPY(M,-D(L,K),  WORK(M+1,1),IONE, X(1,L),IONE)

                CALL DAXPY(M,-B(L,K-1),WORK(1,2),IONE, X(1,L),IONE)
                CALL DAXPY(M,-D(L,K-1),WORK(M+1,2),IONE, X(1,L),IONE)
            END DO
            K = K - 2
        END IF
    END DO
END IF


END SUBROUTINE
