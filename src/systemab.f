C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE SYSTEMBF (NUMNP, MBAND, B, LB, A, LA, MXBND, MXNODE)
      IMPLICIT NONE
      INTEGER MXBND, MXNODE
      INTEGER NUMNP, MBAND, LA, LB
      DOUBLE PRECISION A, B
      DIMENSION A(LA), B(LB)
      DOUBLE PRECISION ALOC, BLOC
      DIMENSION ALOC(LA), BLOC(LB)
      INTEGER I, N, K, L
C----------------------------------------------------------------------
      DO I = 1, LA
         ALOC(I) = A(I)
      END DO
      DO I = 1, LB
         BLOC(I) = B(I)
      END DO
C----------- FORWARD SUBSTITUTION -------------------------------------
      DO 31 N = 1 , NUMNP-1
         DO 21 L = 2 , MBAND
            I = N + L - 1
            IF ( I .GT. NUMNP ) GO TO 21
            BLOC(I) = BLOC(I) - ALOC(N+(L-1)*MXNODE) * BLOC(N)
   21    CONTINUE
         BLOC(N) = BLOC(N) / ALOC(N)
   31 CONTINUE
C---------- BACKSUBSTITUTION -------------
      N=NUMNP-1
      BLOC(NUMNP)=0.D0
   40 DO 50 K = 2, MBAND
         L = N + K - 1
         IF ( L .GT. NUMNP ) GO TO 60
         BLOC(N) = BLOC(N) - ALOC(N+(K-1)*MXNODE) * BLOC(L)
   50 CONTINUE
   60 N = N - 1
      IF ( N .GT. 0 ) GO TO 40
C----------------------------------------------------------------------
      DO I = 1, LB
         B(I) = BLOC(I)
      END DO
      RETURN
      END SUBROUTINE SYSTEMBF
C
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE SYSTEMAF (NUMNP, MBAND, A, LA, MXBND, MXNODE)
      IMPLICIT NONE
      INTEGER MXBND, MXNODE, LA
      INTEGER NUMNP, MBAND
      DOUBLE PRECISION A
      DIMENSION A(LA)
      DOUBLE PRECISION ALOC
      DIMENSION ALOC(LA)
      DOUBLE PRECISION C
      INTEGER I, J, K, L, N
C-----------------------------------------------------------------------
      DO I = 1, LA
         ALOC(I) = A(I)
      END DO
      DO 30 N = 1, NUMNP-1
         DO 20 L = 2, MBAND
            C = ALOC(N+(L-1)*MXNODE) / ALOC(N)
            I = N + L - 1
            IF ( I .GT. NUMNP ) GO TO 20
            J = 0
            DO 10 K = L , MBAND
               J = J + 1
               ALOC(I+(J-1)*MXNODE) = ALOC(I+(J-1)*MXNODE)
     &                                -C*ALOC(N+(K-1)*MXNODE)
   10       CONTINUE
            ALOC(N+(L-1)*MXNODE) = C
   20    CONTINUE
   30 CONTINUE
C-----------------------------------------------------------------------
      DO I =1, LA
         A(I) = ALOC(I)
      END DO
      RETURN
      END SUBROUTINE SYSTEMAF
C
C
