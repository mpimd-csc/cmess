*> \brief \b ZDDOTC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       COMPLEX*16 FUNCTION ZDDOTC(N,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       DOUBLE PRECISION ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZDDOTC forms the dot product of a complex and real vector
*>      ZDOTC = X^H * Y
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in] ZY
*> \verbatim
*>          ZY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*> \author Maximilian Behr
*
*> \date October 2017
*
*
*  =====================================================================
      COMPLEX*16 FUNCTION ZDDOTC(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
      DOUBLE PRECISION ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG
*     ..
      ZTEMP = (0.0d0,0.0d0)
      ZDDOTC = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDDOTC = ZTEMP
      RETURN
      END
