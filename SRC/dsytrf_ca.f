*> \brief \b DSYTRF_CA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DSYTRF_CA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrf_ca.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrf_ca.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrf_ca.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYTRF_CA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N, LDA, LWORK, INFO
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYTRF_CA computes the factorization of a real symmetric matrix A
*> using the Aasen's algorithm.  The form of the factorization is
*>
*>    A = U*T*U**T  or  A = L*T*L**T
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, and T is a symmetric tridiagonal matrix.
*>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the tridiagonal matrix is stored in the diagonals
*>          and the subdiagonals of A just below (or above) the diagonals,
*>          and L is stored below (or above) the subdiaonals, when UPLO
*>          is 'L' (or 'U').
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          On exit, it contains the details of the interchanges, i.e.,
*>          the row and column k of A were interchanged with the
*>          row and column IPIV(k).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.  LWORK >= MAX(1,2*N). For optimum performance
*>          LWORK >= N*(1+NB), where NB is the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
*> \date December 2016
*
*> \ingroup doubleSYcomputational
*
*  =====================================================================
      SUBROUTINE DSYTRF_CA( UPLO, N, NB, A, LDA, TB, LDTB, H, LDH,
     $                      IPIV, INFO)
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, LDA, LDTB, LDH, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), TB( LDTB, *), H( LDH, * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 5.0D-1 )
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, J, K, I1, I2, TD
      INTEGER            NB, KB, NT
      DOUBLE PRECISION   PIV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
*     Determine the block size
*
c      NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
c      NB = 5
      NT = (N+NB-1)/NB
      TD = 2*NB
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF_CA', -INFO )
         RETURN
      END IF
*
*     Quick return
*
      IF ( N.EQ.0 ) THEN
          RETURN
      ENDIF
      KB = MIN(NB, N)
      DO J = 1, KB
         IPIV( J ) = J
      END DO
      IF ( N.LE.KB ) THEN
         CALL DLACPY( 'Full', KB, KB, A( 1, 1 ), LDA,
     $                TB( TD+1, 1 ), LDTB-1 )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        .....................................................
*        Factorize A as L*D*L**T using the upper triangle of A
*        .....................................................
*
      ELSE
*
*        .....................................................
*        Factorize A as L*D*L**T using the lower triangle of A
*        .....................................................
*
         DO J = 0, NT-1
*         
*           Generate Jth column of W and H
* 
            KB = MIN(NB, N-J*NB)
            DO I = 1, J-1
*              H(I,J) = T(I,I)*L(J,I)'
               CALL DGEMM( 'NoTranspose', 'Transpose',
     $                      NB, KB, NB,
     $                      ONE,  TB( TD+1, I*NB+1 ), LDTB-1,
     $                            A( J*NB+1, (I-1)*NB+1 ), LDA,
     $                      ZERO, H( I*NB+1, 1 ), LDH )
*              H(I,J) += Z, where Z = T(I+1,I)'*L(J,I+1)'
               CALL DGEMM( 'Transpose', 'Transpose',
     $                      NB, KB, NB,
     $                      ONE, TB( TD+NB+1, I*NB+1 ), LDTB-1,
     $                           A( J*NB+1, I*NB+1 ), LDA,
     $                      ONE, H( I*NB+1, 1 ), LDH )
               IF( I.GT.1 ) THEN
*                 H(I,J) += X where X = T(I,I-1)*L(J,I-1)'
                  CALL DGEMM( 'NoTranspose', 'Transpose',
     $                         NB, KB, NB,
     $                         ONE, TB( TD+NB+1, (I-1)*NB+1 ), LDTB-1,
     $                              A( J*NB+1, (I-2)*NB+1 ), LDA,
     $                         ONE, H( I*NB+1, 1 ), LDH )
               END IF
            END DO
*         
*           Compute T(J,J)
*     
            CALL DLACPY( 'Full', KB, KB, A( J*NB+1, J*NB+1 ), LDA,
     $                   TB( TD+1, J*NB+1 ), LDTB-1 ) 
            IF( J.GT.1 ) THEN
*              T(J,J) = L(J,1:J)*H(1:J)             
               CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                      KB, KB, (J-1)*NB,
     $                     -ONE, A( J*NB+1, 1 ), LDA,
     $                           H( NB+1, 1 ), LDH,
     $                      ONE, TB( TD+1, J*NB+1 ), LDTB-1 )
*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                      KB, NB, KB,
     $                      ONE,  A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                            TB( TD+NB+1, (J-1)*NB+1 ), LDTB-1,
     $                      ZERO, H( 1, 1 ), LDH )
               CALL DGEMM( 'NoTranspose', 'Transpose',
     $                      KB, KB, NB,
     $                     -ONE, H( 1, 1 ), LDH,
     $                           A( J*NB+1, (J-2)*NB+1 ), LDA,
     $                      ONE, TB( TD+1, J*NB+1 ), LDTB-1 )
            END IF
            IF( J.GT.0 ) THEN 
               CALL DSYGST( 1, 'Lower', KB, TB( TD+1, J*NB+1 ), LDTB-1, 
     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, INFO )
            END IF
*
*           Expand T(J,J) into full format
*
            DO I = 1, KB
               DO K = I+1, KB
                  TB( TD-(K-(I+1)), J*NB+K ) = TB( TD+(K-I)+1, J*NB+I )
               END DO
            END DO

            IF( J.LT.NT-1 ) THEN
               IF( J.GT.0 ) THEN
*
*                 Compute H(J,J)
*
                  CALL DGEMM( 'NoTranspose', 'Transpose',
     $                         KB, KB, KB,
     $                         ONE,  TB( TD+1, J*NB+1 ), LDTB-1,
     $                               A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                         ZERO, H( J*NB+1, 1 ), LDH )
                  IF( J.GT.1 ) THEN 
                     CALL DGEMM( 'NoTranspose', 'Transpose',
     $                            KB, KB, KB,
     $                           ONE, TB( TD+NB+1, (J-1)*NB+1 ), LDTB-1,
     $                                 A( J*NB+1, (J-2)*NB+1 ), LDA,
     $                            ONE, H( J*NB+1, 1 ), LDH )
                  END IF
*
*                 Update with the previous column
*
                  CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                         N-(J+1)*NB, NB, J*NB,
     $                        -ONE, A( (J+1)*NB+1, 1 ), LDA,
     $                              H( NB+1, 1 ), LDH,
     $                         ONE, A( (J+1)*NB+1, J*NB+1 ), LDA )
               END IF
               CALL DGETRF( N-(J+1)*NB, NB, 
     $                      A( (J+1)*NB+1, J*NB+1 ), LDA,
     $                      IPIV( (J+1)*NB+1 ), INFO )
               IF (INFO.NE.0) THEN
                  WRITE(*,*) 'DGETRF returned INFO=',INFO,' at J=',J
               END IF
*         
*              Compute T(J+1, J)     
*     
               KB = MIN(NB, N-(J+1)*NB)
               CALL DLACPY( 'Upper', KB, NB,
     $                      A( (J+1)*NB+1, J*NB+1 ), LDA,
     $                      TB( TD+NB+1, J*NB+1 ), LDTB-1 )
               IF( J.GT.0 ) THEN 
                  CALL DTRSM( 'R', 'L', 'T', 'U', KB, NB, ONE,
     $                        A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                        TB( TD+NB+1, J*NB+1 ), LDTB-1 )
               END IF
*
*              Copy T(J+1,J) into T(J, J+1)
*
               DO K = 1, NB
                  DO I = 1, MIN(K, KB)
                     TB( TD-NB+K-I+1, J*NB+NB+I ) =
     $                  TB( TD+NB+I-K+1, J*NB+K )
                  END DO
               END DO
               CALL DLASET( 'Upper', KB, NB, ZERO, ONE, 
     $                      A( (J+1)*NB+1, J*NB+1), LDA )
*              
*              Apply pivots to trailing submatrix of A
*     
               DO K = 1, KB
*                 > Adjust ipiv               
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
*                  
                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
                  IF( I1.NE.I2 ) THEN 
*                    > Apply pivots to previous columns of L
                     CALL DSWAP( K-1, A( I1, (J+1)*NB+1 ), LDA, 
     $                                A( I2, (J+1)*NB+1 ), LDA )
*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)               
                     CALL DSWAP( I2-I1-1, A( I1+1, I1 ), 1,
     $                                    A( I2, I1+1 ), LDA )
*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     CALL DSWAP( N-I2, A( I2+1, I1 ), 1,
     $                                 A( I2+1, I2 ), 1 ) 
*                    > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 )
                     A( I1, I1 ) = A( I2, I2 )
                     A( I2, I2 ) = PIV
                  ENDIF   
               END DO   
*         
*              Apply pivots to previous columns of L
*         
               CALL DLASWP( J*NB, A( 1, 1 ), LDA, 
     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            END IF
         END DO
      END IF
*
*     End of DSYTRF_CA
*
      END
