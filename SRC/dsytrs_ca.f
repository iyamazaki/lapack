*> \brief \b DSYTRS_CA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DSYTRS_CA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs_ca.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs_ca.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs_ca.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*      SUBROUTINE DSYTRS_CA( UPLO, N, NB, NRHS, A, LDA, TB, LDTB, IPIV, 
*                            IPIV2, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N, NB, NRHS, LDA, LDTB, LDB, INFO
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), IPIV2( * )
*       DOUBLE PRECISION   A( LDA, * ), TB( LDTB, * ), B( LDB, * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYTRS_CA solves a system of linear equations A*X = B with a real
*> symmetric matrix A using the factorization A = U*T*U**T or
*> A = L*T*L**T computed by DSYTRF_CA.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*T*U**T;
*>          = 'L':  Lower triangular, form is A = L*T*L**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The bandwidth of the matrix T.  NB > 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          Details of factors computed by DSYTRF_CA.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TB
*> \verbatim
*>          TB is DOUBLE PRECISION array, dimension (LDTB, N)
*>          Details of factors computed by DSYTRF_CA.
*> \endverbatim
*>
*> \param[in] LDTB
*> \verbatim
*>          The leading dimension of the array TB. LDTB >= 3*NB+1.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges as computed by DSYTRF_CA.
*> \endverbatim
*>
*> \param[in] IPIV2
*> \verbatim
*>          IPIV2 is INTEGER array, dimension (N)
*>          Details of the interchanges as computed by DSYTRF_CA.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup doubleSYcomputational
*
*  =====================================================================
      SUBROUTINE DSYTRS_CA( UPLO, N, NB, NRHS, A, LDA, TB, LDTB, IPIV, 
     $                      IPIV2, B, LDB, INFO )
*
*  -- LAPACK computational routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, NB, NRHS, LDA, LDTB, LDB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IPIV2( * )
      DOUBLE PRECISION   A( LDA, * ), TB( LDTB, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGBTRS, DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRS_CA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U*T*U**T.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B
*
            CALL DLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (U**T \P**T * B) -> B    [ (U**T \P**T * B) ]
*
            CALL DTRSM( 'L', 'U', 'T', 'U', N-NB, NRHS, ONE, A(1, NB+1),
     $                 LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
*
         CALL DGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB,
     $               INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]
*
            CALL DTRSM( 'L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1),
     $                  LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
*
            CALL DLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
*
      ELSE
*
*        Solve A*X = B, where A = L*T*L**T.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B
*
            CALL DLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
*
            CALL DTRSM( 'L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1),
     $                 LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL DGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB,
     $               INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
*
            CALL DTRSM( 'L', 'L', 'T', 'U', N-NB, NRHS, ONE, A(NB+1, 1),
     $                  LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
*
            CALL DLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
      END IF
*
      RETURN
*
*     End of DSYTRS_CA
*
      END
