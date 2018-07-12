!> \brief \b DGESC2X solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESC2X( TRANS, N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDA, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESC2X solves a system of linear equations
!>
!>           op(A) * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by DGETC2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Decides whether op(A) = A or op(A) = A^T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by DGETC2:  A = P * L * U * Q
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is DOUBLE PRECISION array, dimension (N).
!>          On entry, the right hand side vector b.
!>          On exit, the solution vector X.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE contains the scale factor. SCALE is chosen
!>          0 <= SCALE <= 1 to prevent owerflow in the solution.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Peter Benner, Martin Koehler, MPI Magdeburg
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date January 2014
!
!> \ingroup doubleGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!>     Extention to solve with op(A), Peter Benner, Martin Koehler MPI Magdeburg
!
!  =====================================================================
      SUBROUTINE DGESC2X( TRANS, N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDA, N
      DOUBLE PRECISION   SCALE
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), RHS( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, EPS, SMLNUM, TEMP
      LOGICAL NOTRAN

!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASWP, DSCAL
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      LOGICAL            LSAME
      EXTERNAL           IDAMAX, DLAMCH, LSAME

!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!      Set constant to control owerflow
!
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      NOTRAN = LSAME(TRANS,'N')

      IF ( NOTRAN ) THEN
        !
        !     Apply permutations IPIV to RHS
        !
              CALL DLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
        !
        !     Solve for L part
        !
              DO 20 I = 1, N - 1
                 DO 10 J = I + 1, N
                    RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
           10    CONTINUE
           20 CONTINUE
        !
        !     Solve for U part
        !
              SCALE = ONE
        !
        !     Check for scaling
        !
              I = IDAMAX( N, RHS, 1 )
              IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
                 TEMP = ( ONE / TWO ) / ABS( RHS( I ) )
                 CALL DSCAL( N, TEMP, RHS( 1 ), 1 )
                 SCALE = SCALE*TEMP
              END IF
        !
              DO 40 I = N, 1, -1
                 TEMP = ONE / A( I, I )
                 RHS( I ) = RHS( I )*TEMP
                 DO 30 J = I + 1, N
                    RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
           30    CONTINUE
           40 CONTINUE
        !
        !     Apply permutations JPIV to the solution (RHS)
        !
              CALL DLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      ELSE
        !
        !     Apply permutations JPIV to RHS
        !
              CALL DLASWP( 1, RHS, LDA, 1, N-1, JPIV, 1 )
        !
        !     Solve for U part
        !
              SCALE = ONE
        !
        !     Check for scaling
        !
              I = IDAMAX( N, RHS, 1 )
              IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
                 TEMP = ( ONE / TWO ) / ABS( RHS( I ) )
                 CALL DSCAL( N, TEMP, RHS( 1 ), 1 )
                 SCALE = SCALE*TEMP
              END IF
        !
              DO 140 I = 1, N
                 TEMP = ONE / A( I, I )
                 RHS( I ) = RHS( I )*TEMP
                 DO 130 J = I +1 , N
                    RHS( J ) = RHS( J ) - RHS( I  )* A( I, J )
           130    CONTINUE
           140 CONTINUE
        !
        !     Solve for L part
        !
              DO 120 I = N,1,-1
                 DO 110 J = I-1, 1, -1
                    RHS( J ) = RHS( J ) - A( I , J )*RHS( I )
           110    CONTINUE
           120 CONTINUE
        !
        !     Apply permutations JPIV to the solution (RHS)
        !
              CALL DLASWP( 1, RHS, LDA, 1, N-1, IPIV, -1 )
      END IF

      RETURN
!
!     End of DGESC2X
!
      END
