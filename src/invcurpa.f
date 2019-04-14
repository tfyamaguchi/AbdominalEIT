C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE INVCURPA (NCUR, NEX, NEZ, VOL, EXYZ, EXYZ0,
     & SIGXY, SIGZ, POD, SENXY,SENSZ, MXE, MZE, MXNODE, JCUR, KVOL,GW)
      IMPLICIT NONE
      INTEGER MXE,MZE,MXNODE,JCUR,KVOL
      INTEGER NCUR,NEX,NEZ
      INTEGER NCURL,NEXL,NEZL
      DOUBLE PRECISION VOL,EZFAC,EXYZ,EXYZ0,SIGXY,SIGZ,POD,SENXY,SENXYL
      DOUBLE PRECISION SENSZ,SENSZL,GW,GWL
      DOUBLE PRECISION VOLL,EZFACL,EXYZL,EXYZ0L,SIGXYL,SIGZL,PODL,SUM
      DIMENSION VOL(MZE*MXE),VOLL(MZE,MXE),
     & EXYZ(3*MZE*MXE), EXYZL(3,MZE,MXE),
     & EXYZ0(3*MZE*MXE*JCUR),EXYZ0L(3,MZE,MXE,JCUR),
     & SIGXY(MXE),SIGXYL(MXE),SIGZ(MXE),SIGZL(MXE),POD(KVOL*JCUR),
     & PODL(KVOL,JCUR),
     & SENXY(MXE*JCUR),SENXYL(MXE,JCUR),SENSZ(MXE*JCUR),
     & SENSZL(MXE,JCUR),GWL(KVOL,JCUR),GW(KVOL*JCUR)
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
         SIGXYL(I) = SIGXY(I)
         SIGZL(I) = SIGZ(I)
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
            SENXYL(N,J) = SENXY((J-1)*MXE + N)
            SENSZL(N,J) = SENSZ((J-1)*MXE + N)         
         END DO
      END DO
C----------- CALCUALTION  -------------------------------------
      DO J = 1, JCUR
         DO N = 1, NEXL
            SUM=0.D0
            DO L = 1, NEZL
               SUM=SUM+2.D0*VOLL(L,N)*(EXYZL(1,L,N)*EXYZ0L(1,L,N,J)
     &         + EXYZL(2,L,N)*EXYZ0L(2,L,N,J))
            END DO
            SENXYL(N,J) = -SUM*SIGXYL(N)/PODL(NCURL,J)*GWL(NCURL,J)
            SUM=0.D0
            DO L = 1, NEZL
               SUM=SUM+2.D0*VOLL(L,N)*(EXYZL(3,L,N)*EXYZ0L(3,L,N,J))
            END DO
            SENSZL(N,J) = -SUM*SIGZL(N)/PODL(NCURL,J)*GWL(NCURL,J)
         END DO
      END DO
      DO J = 1, JCUR
         DO N = 1, MXE
            SENXY((J-1)*MXE + N) = SENXYL(N,J)
            SENSZ((J-1)*MXE + N) = SENSZL(N,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE INVCURPA
C============================================================================
