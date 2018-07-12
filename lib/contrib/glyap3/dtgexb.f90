SUBROUTINE DTGEXB( N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NB, WORK, LWORK, INFO )
!     PURPOSE
!
!     Reorder the eigenvalues of a matrix pencil (A,B) such that
!     the complex eigenvalues pairs will not cause a change of the block
!     size in blocked linear matrix equation solver codes.
!
!     The function uses DTGEX2 from LAPACK for moving the eigenvalues
!     on the diagonal.
!
!     ARGUMENTS
!
!
!
!    Input/Output Parameters
!
!     N       (input) INTEGER
!             The order of the matrix A.  N >= 0.
!
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the leading N-by-N upper
!             Hessenberg part of this array must contain the
!             generalized Schur factor A_s of the matrix A.
!             On exit, the leading N-by-N part of this array contains
!             the generalized Schur factor A_s of the matrix A with
!             no complex eigenvalues pairs an A(k*NB,k*NB), where k is
!             a non negative integer.
!
!     LDA     INTEGER
!             The leading dimension of the array A.  LDA >= MAX(1,N).
!
!     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
!             On entry, the leading N-by-N upper
!             triangular part of this array must contain the
!             generalized Schur factor E_s of the matrix E,
!             On exit, the leading N-by-N part of this array contains
!             the generalized Schur factor E_s of the matrix E
!             with reordered eigenvalues.
!
!     LDE     INTEGER
!             The leading dimension of the array E.  LDE >= MAX(1,N).
!
!     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!             On entry, the leading N-by-N part of
!             this array must contain the orthogonal matrix Q from
!             the generalized Schur factorization.
!             On exit, the leading N-by-N part of this array contains
!             the updated orthogonal matrix Q from the generalized Schur
!             factorization with integrated eigenvalue reordering.
!
!     LDQ     INTEGER
!             The leading dimension of the array Q.  LDQ >= MAX(1,N).
!
!     Z       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
!             On entry, the leading N-by-N part of
!             this array must contain the orthogonal matrix Z from
!             the generalized Schur factorization.
!             On exit, the leading N-by-N part of this array contains
!             the updated orthogonal matrix Z from the generalized Schur
!             factorization with integrated eigenvalue reordering.
!
!     LDZ     INTEGER
!             The leading dimension of the array Z.  LDZ >= MAX(1,N).
!
!
!     NB      INTEGER
!             Block size of the solver planed to use. Typical values
!             are 16, 32, or 64. The block size must be an even postive
!             integer otherwise an error is returned.
!
!     Workspace
!
!     WORK    DOUBLE PRECISION array, dimension (LWORK)
!             Workspace
!     LWORK   INTEGER
!             Size of the workspace.
!             LWORK >= MAX(1, 4*N, 32)
!
!     Error indicator
!
!     INFO    INTEGER
!             = 0:  successful exit;
!             < 0:  if INFO = -i, the i-th argument had an illegal
!                   value;
!             = 1:  The internal DTGEX2 caused an error. ]
!
!
!
    IMPLICIT NONE

    ! .. Arguments ..
    INTEGER INFO, LDA, LDB, LDQ, LDZ, LWORK, N, NB
    DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )

    ! .. Local Variables ..
    INTEGER K, POS, SJ, J, INFO2

    DOUBLE PRECISION ZERO
    PARAMETER (ZERO = 0.0d0)
    LOGICAL WANT
    PARAMETER (WANT = .TRUE.)

    ! .. External Functions
    INTRINSIC MAX, MIN, MOD
    EXTERNAL XERBLA, DTGEX2


    ! .. Check Parameters
    INFO = 0
    IF ( N .LT. 0) THEN
        INFO = -1
    ELSE IF ( LDA .LT. MAX(1, N))  THEN
        INFO = -3
    ELSE IF ( LDB .LT. MAX(1, N))  THEN
        INFO = -5
    ELSE IF ( LDQ .LT. MAX(1, N))  THEN
        INFO = -7
    ELSE IF ( LDZ .LT. MAX(1, N))  THEN
        INFO = -9
    ELSE IF (( MOD(NB,2) .EQ. 1 ) .OR. (NB .LT. 1)) THEN
        INFO = -10
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -10
    ELSE IF ( LWORK .LT. MAX(1, N*4, 32)) THEN
        INFO = -12
    END IF

    IF ( INFO .NE. 0 ) THEN
        CALL XERBLA('DTGEXB', -INFO)
    END IF


    ! .. Sort the Eigenvalues on the Block Boundary..
    DO K = 1, N, NB
        POS = K+NB-1
        IF ( POS .GE. N ) THEN
            EXIT
        END IF

        IF ( A(POS+1, POS) .NE. ZERO) THEN
            SJ = K
            ! .. Search for the last previous real eigenvalue ..
            DO J = POS-2, K, -2
                IF( J .GT. 0 .AND. A(J+1,J) .EQ. ZERO) THEN
                    SJ = J + 1
                    EXIT
                END IF
            END DO

            ! .. Move the real eigenvalue to the blocksize boundary ..
            DO WHILE ( SJ .LT. POS )
                CALL DTGEX2(WANT, WANT, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, SJ, 1, 2, WORK, LWORK, INFO2)
                SJ = SJ + 2
                IF ( INFO2 .NE. 0 ) THEN
                    INFO = 1
                    RETURN
                END IF
            END DO
        END IF
    END DO


END SUBROUTINE
