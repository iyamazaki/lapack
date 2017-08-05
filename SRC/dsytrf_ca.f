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
      SUBROUTINE DSYTRF_CA( UPLO, N, NB, A, LDA, T, LDT, H, LDH,
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
      INTEGER            N, LDA, LDT, LDH, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), T( LDT, *), H( LDH, * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 5.0D-1 )
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, J, K, I1, I2
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
c      WRITE(*,*) 'NB=',NB,'NT=',NT
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
      IPIV( 1 ) = 1
      IF ( N.EQ.1 ) THEN
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
c            WRITE(11,*) '--',J,'--'
            DO I = 1, J-1
*              H(I,J) = T(I,I)*L(J,I)'
               CALL DGEMM( 'NoTranspose', 'Transpose',
     $                      NB, KB, NB,
     $                      ONE,  T( I*NB+1, I*NB+1 ), LDT,
     $                            A( J*NB+1, (I-1)*NB+1 ), LDA,
     $                      ZERO, H( I*NB+1, 1 ), LDH )
c               WRITE(11,*) 'H',I*NB+1
c               DO K = 1, NB
c                  WRITE(11,*) H(I*NB+K, 1:NB)
c               END DO 
c               WRITE(11,*)
*              H(I,J) += Z, where Z = T(I+1,I)'*L(J,I+1)'
               CALL DGEMM( 'Transpose', 'Transpose',
     $                      NB, KB, NB,
     $                      ONE, T( (I+1)*NB+1, I*NB+1 ), LDT,
     $                           A( J*NB+1, I*NB+1 ), LDA,
     $                      ONE, H( I*NB+1, 1 ), LDH )
c               DO K = 1, NB
c                  WRITE(11,*) H(I*NB+K, 1:NB)
c               END DO 
c               WRITE(11,*)
               IF( I.GT.1 ) THEN
*                 H(I,J) += X where X = T(I,I-1)*L(J,I-1)'
c                  WRITE(11,*) ' xxxxxxxxxxxxxxx'
c                  DO K = 1, NB
c                     WRITE(11,*) T(I*NB+K, (I-1)*NB+1:(I-1)*NB+NB)
c                  END DO 
c                  WRITE(11,*)
c                  DO K = 1, NB
c                     WRITE(11,*) A(J*NB+K, (I-2)*NB+1:(I-2)*NB+NB)
c                  END DO 
c                  WRITE(11,*)
                  CALL DGEMM( 'NoTranspose', 'Transpose',
     $                         NB, KB, NB,
     $                         ONE, T( I*NB+1, (I-1)*NB+1 ), LDT,
     $                              A( J*NB+1, (I-2)*NB+1 ), LDA,
     $                         ONE, H( I*NB+1, 1 ), LDH )
               END IF
            END DO
*         
*           Compute T(J,J)
*     
c            WRITE(11,*) ' >> T(',J,',',J,')',NB
            CALL DLACPY( 'Full', KB, KB, A( J*NB+1, J*NB+1 ), LDA,
     $                   T( J*NB+1, J*NB+1 ), LDT ) 
            IF( J.GT.1 ) THEN
c               WRITE(11,*) ' A(',J*NB+1,',1)'
c               DO K = 1, KB
c                  WRITE(11,*) A(J*NB+K, 1:J*NB)
c               END DO 
c               WRITE(11,*)
c               WRITE(11,*) ' H(',NB+1,',1)'
c               DO K = 1, J*NB
c                  WRITE(11,*) H(NB+K, 1:KB)
c               END DO 
c               WRITE(11,*)
c               WRITE(11,*) ' T(',J*NB+1,',',J*NB+1,')'
c               DO K = 1, NB
c                  WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c               END DO 
c               WRITE(11,*)
*              T(J,J) = L(J,1:J)*H(1:J)             
               CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                      KB, KB, (J-1)*NB,
     $                     -ONE, A( J*NB+1, 1 ), LDA,
     $                           H( NB+1, 1 ), LDH,
     $                      ONE, T( J*NB+1, J*NB+1 ), LDT )
c               DO K = 1, NB
c                  WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c               END DO 
*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
c               WRITE(11,*) 'A(',J,',',J-1,'), T(',J,',',J-1,
c     $                     '), A(',J,',',J-2,')'
               CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                      KB, NB, KB,
     $                      ONE,  A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                            T( J*NB+1, (J-1)*NB+1 ), LDT,
     $                      ZERO, H( 1, 1 ), LDH )
c               DO K = 1, NB
c                  WRITE(11,*) H(K, 1:NB)
c               END DO 
c               WRITE(11,*)
c               DO K = 1, NB
c                  WRITE(11,*) A(J*NB+K, (J-2)*NB+1:(J-1)*NB)
c               END DO 
               CALL DGEMM( 'NoTranspose', 'Transpose',
     $                      KB, KB, NB,
     $                     -ONE, H( 1, 1 ), LDH,
     $                           A( J*NB+1, (J-2)*NB+1 ), LDA,
     $                      ONE, T( J*NB+1, J*NB+1 ), LDT )
c               DO K = 1, NB
c                  WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c               END DO 
            END IF
            IF( J.GT.0 ) THEN 
c               WRITE(11,*) 'T(',J,',',J,') / A(',J,',',J-1,')'
c               DO K = 1, NB
c                  WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c               END DO 
c               WRITE(11,*)
c               DO K = 1, NB
c                  WRITE(11,*) A(J*NB+K, (J-1)*NB+1:(J-1)*NB+NB)
c               END DO 
c               WRITE(11,*)
               CALL DSYGST( 1, 'Lower', KB, T( J*NB+1, J*NB+1 ), LDT, 
     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, INFO )
c               DO K = 1, NB
c                  WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c               END DO 
            END IF
*
*           Expand T(J,J) into full format
*
            DO I = 1, KB
               DO K = I+1, KB
                  T( J*NB+I, J*NB+K ) = T( J*NB+K, J*NB+I )
               END DO
            END DO

            IF( J.LT.NT-1 ) THEN
               IF( J.GT.0 ) THEN
*
*                 Compute H(J,J)
*
                  CALL DGEMM( 'NoTranspose', 'Transpose',
     $                         KB, KB, KB,
     $                         ONE,  T( J*NB+1, J*NB+1 ), LDT,
     $                               A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                         ZERO, H( J*NB+1, 1 ), LDH )
                  IF( J.GT.1 ) THEN 
                     CALL DGEMM( 'NoTranspose', 'Transpose',
     $                            KB, KB, KB,
     $                            ONE, T( J*NB+1, (J-1)*NB+1 ), LDT,
     $                                 A( J*NB+1, (J-2)*NB+1 ), LDA,
     $                            ONE, H( J*NB+1, 1 ), LDH )
                  END IF
c                  WRITE(11,*) 'H',J
c                  DO K = 1, NB
c                     WRITE(11,*) T(J*NB+K, J*NB+1:J*NB+NB)
c                  END DO 
c                  WRITE(11,*)
c                  DO K = 1, NB
c                     WRITE(11,*) A(J*NB+K, (J-1)*NB+1:(J-1)*NB+NB)
c                  END DO 
c                  WRITE(11,*)
c                  DO K = 1, NB
c                     WRITE(11,*) H(J*NB+K, 1:NB)
c                  END DO 
*
*                 Update with the previous column
*
c                  WRITE(11,*) 'UPDATE',N-(J+1)*NB,J*NB
c                  DO K=(J+1)*NB+1, N
c                     WRITE(11,*) A(K,1:(J+1)*NB)
c                  END DO    
                  CALL DGEMM( 'NoTranspose', 'NoTranspose',
     $                         N-(J+1)*NB, NB, J*NB,
     $                        -ONE, A( (J+1)*NB+1, 1 ), LDA,
     $                              H( NB+1, 1 ), LDH,
     $                         ONE, A( (J+1)*NB+1, J*NB+1 ), LDA )
               END IF
c               WRITE(11,*) J,'DGETRF',N-(J+1)*NB, NB
c               DO K = (J+1)*NB+1, N
c                  WRITE(11,*) A(K, J*NB+1:(J+1)*NB)
c               END DO 
c               WRITE(11,*)
               CALL DGETRF( N-(J+1)*NB, NB, 
     $                      A( (J+1)*NB+1, J*NB+1 ), LDA,
     $                      IPIV( (J+1)*NB+1 ), INFO )
c               WRITE(11,*) 'INFO=',INFO
               IF (INFO.NE.0) THEN
                  WRITE(*,*) 'DGETRF returned INFO=',INFO,' at J=',J
               END IF
c               DO K = (J+1)*NB+1, N
c                  WRITE(*,*) K,' ',A(K, J*NB+1:(J+1)*NB)
c               END DO 
*         
*              Compute T(J+1, J)     
*     
c               WRITE(*,*) 'T(',J+1,',',J,')'
               KB = MIN(NB, N-(J+1)*NB)
               CALL DLACPY( 'Upper', KB, NB,
     $                      A( (J+1)*NB+1, J*NB+1 ), LDA,
     $                      T( (J+1)*NB+1, J*NB+1 ), LDT )
               IF( J.GT.0 ) THEN 
                  CALL DTRSM( 'R', 'L', 'T', 'U', KB, NB, ONE,
     $                        A( J*NB+1, (J-1)*NB+1 ), LDA,
     $                        T( (J+1)*NB+1, J*NB+1 ), LDT )
               END IF
               CALL DLASET( 'Upper', KB, NB, ZERO, ONE, 
     $                      A( (J+1)*NB+1, J*NB+1), LDA )
*              
*              Apply pivots to trailing submatrix of A
*     
c               WRITE(11,*) 'PIVOT' 
c               DO K = 1, N
c                  WRITE(11,*) A(K, 1:N)
c               END DO 
c               WRITE(11,*)
               DO K = 1, KB
*                 > Adjust ipiv               
c                  WRITE(11,*) 'IPIV',(J+1)*NB+K
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
*                  
                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
c                  WRITE(11,*) I1,I2
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
c                     DO I1 = 1, N
c                        WRITE(11,*) A(I1, 1:N)
c                     END DO 
                  ENDIF   
c                  WRITE(11,*) 'A(',(J+1)*NB+1,')'
c                  DO I1 = 1, NB
c                     WRITE(11,*) A((J+1)*NB+I1,
c     $                             (J+1)*NB+1:(J+2)*NB)
c                  END DO
               END DO   
*         
*              Apply pivots to previous columns of L
*         
               CALL DLASWP( J*NB, A( 1, 1 ), LDA, 
     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
c               DO K = 1, N
c                  WRITE(11,*) A(K, 1:N)
c               END DO 
            END IF
         END DO
      END IF
*
*     End of DSYTRF_CA
*
      END
