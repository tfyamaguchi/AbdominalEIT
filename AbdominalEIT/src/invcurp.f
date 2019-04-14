C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE INVCURP (NCUR, NEX, NEZ, VOL, EXYZ, EXYZ0,
     & SIGMA, POD, SENSEMAT, MXE, MZE, MXNODE, JCUR, KVOL, GW)
      IMPLICIT NONE
      INTEGER MXE,MZE,MXNODE,JCUR,KVOL
      INTEGER NCUR,NEX,NEZ
      INTEGER NCURL,NEXL,NEZL
      DOUBLE PRECISION VOL,EXYZ,EXYZ0,SIGMA,POD,SENSEMAT,SENSEMATL
      DOUBLE PRECISION VOLL,EXYZL,EXYZ0L,SIGMAL,PODL,SUM,GW,GWL
      DIMENSION VOL(MZE*MXE),VOLL(MZE,MXE),GW(KVOL*JCUR),GWL(KVOL,JCUR),
     & EXYZ(3*MZE*MXE), EXYZL(3,MZE,MXE),
     & EXYZ0(3*MZE*MXE*JCUR),EXYZ0L(3,MZE,MXE,JCUR),
     & SIGMA(MXE),SIGMAL(MXE),POD(KVOL*JCUR),PODL(KVOL,JCUR),
     & SENSEMAT(MXE*JCUR),SENSEMATL(MXE,JCUR)
      INTEGER I, J, K, L, N
C----------------------------------------------------------------------
      NCURL = NCUR
      NEXL = NEX
      NEZL = NEZ
      DO J = 1, MXE
         DO I = 1, MZE
            VOLL(I,J) = VOL((J-1)*MZE + I)
         END DO
      END DO
      DO I = 1, MXE
         SIGMAL(I) = SIGMA(I)
      END DO
      DO K = 1, MXE
         DO J= 1, MZE
            DO I = 1, 3
               EXYZL(I,J,K)= EXYZ((K-1)*3*MZE + (J-1)*3 + I)
            END DO
         END DO
      END DO
      DO L = 1, JCUR
         DO K = 1, MXE
            DO J= 1, MZE
               DO I = 1, 3
                  EXYZ0L(I,J,K,L)=EXYZ0((L-1)*3*MZE*MXE + (K-1)*3*MZE
     &                            + (J-1)*3 + I)
               END DO
            END DO
         END DO
      END DO
      DO J = 1, JCUR
         DO I = 1, KVOL
            PODL(I,J) = POD((J-1)*KVOL + I)
            GWL(I,J) = GW((J-1)*KVOL + I)
         END DO
      END DO
      DO J = 1, JCUR
         DO N = 1, MXE
            SENSEMATL(N,J) = SENSEMAT((J-1)*MXE + N)
         END DO
      END DO
C----------- CALCUALTION  -------------------------------------
      DO J = 1, JCUR
         DO N = 1, NEXL
            SUM=0.D0
            DO L = 1, NEZL
               SUM=SUM+2.D0*VOLL(L,N)*(EXYZL(1,L,N)*EXYZ0L(1,L,N,J)
     &         + EXYZL(2,L,N)*EXYZ0L(2,L,N,J)
     &         + EXYZL(3,L,N)*EXYZ0L(3,L,N,J))
            END DO
            SENSEMATL(N,J) = -SUM*SIGMAL(N)/PODL(NCURL,J)*GWL(NCURL,J)
         END DO
      END DO
      DO J = 1, JCUR
         DO N = 1, MXE
            SENSEMAT((J-1)*MXE + N) = SENSEMATL(N,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE INVCURP
C============================================================================
