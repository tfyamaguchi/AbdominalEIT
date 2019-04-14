C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE FIELDF (NEX, NEZ, LCEX, CEX, LBB, BB, NNODEZ,
     & LNODEX, NODEX, LNODEZ, NODEZ, LEXYZ, EXYZ, MXE, MZE, MXNODE)
      IMPLICIT NONE
      INTEGER MXE, MZE, MXNODE
      INTEGER  NEX,NEZ,LCEX,LBB,NNODEZ,LNODEX,NODEX,LNODEZ,NODEZ,LEXYZ
      INTEGER  NEXL,NEZL,LCEXL,LBBL,NNODEZL,LNODEXL,LNODEZL,LEXYZL
      INTEGER  NODEZLOC,NODEXLOC
      DOUBLE PRECISION CEX,CEXLOC,BB,BBLOC,EXYZ,EXYZLOC,SUM
      DIMENSION BB(MXNODE),BBLOC(MXNODE),
     &                 NODEX(4*MZE*MXE), NODEXLOC(4,MZE,MXE),
     &                 NODEZ(4*MZE*MXE), NODEZLOC(4,MZE,MXE),
     &                 EXYZ(3*MZE*MXE), EXYZLOC(3,MZE,MXE),
     &                 CEX(3*4*MZE*MXE), CEXLOC(3,4,MZE,MXE)
      INTEGER I, J, K, L, N
C----------------------------------------------------------------------
      NEXL = NEX
      NEZL = NEZ
      LCEXL = LCEX
      LBBL = LBB
      NNODEZL = NNODEZ
      LNODEXL = LNODEX
      LNODEZL = LNODEZ
      LEXYZL = LEXYZ
      DO I = 1, MXNODE
         BBLOC(I) = BB(I)
      END DO
      DO L = 1, MXE
         DO K = 1, MZE
            DO J= 1, 4
               DO I = 1, 3
                  CEXLOC(I,J,K,L)=CEX((L-1)*3*4*MZE + (K-1)*3*4
     &                            + (J-1)*3 + I)
               END DO
            END DO
         END DO
      END DO
      DO K = 1, MXE
         DO J= 1, MZE
            DO I = 1, 3
               EXYZLOC(I,J,K)= EXYZ((K-1)*3*MZE + (J-1)*3 + I)
            END DO
         END DO
      END DO
      DO K = 1, MXE
         DO J= 1, MZE
            DO I = 1, 4
               NODEXLOC(I,J,K) = NODEX((K-1)*4*MZE + (J-1)*4 + I)
               NODEZLOC(I,J,K) = NODEZ((K-1)*4*MZE + (J-1)*4 + I)
            END DO
         END DO
      END DO
C----------- CALCUALTION  ----------------------------------------
      DO N = 1, NEXL
         DO L= 1,NEZL
            DO J= 1,3
               SUM=0.D0
               DO K=1,4
               SUM=SUM-CEXLOC(J,K,L,N)*BBLOC(NNODEZL*(NODEXLOC(K,L,N)-1)
     &               +NODEZLOC(K,L,N))
               END DO
               EXYZLOC(J,L,N) = SUM
            END DO
         END DO
      END DO
C------------ WRITE BACK -------------------------------------------------
      DO K = 1, MXE
         DO J= 1, MZE
            DO I = 1, 3
               EXYZ((K-1)*3*MZE + (J-1)*3 + I)= EXYZLOC(I,J,K)
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE FIELDF
C============================================================================
