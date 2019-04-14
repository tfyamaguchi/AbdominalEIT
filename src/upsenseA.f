C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE UPSENSEA (NEE,NUMPOL,SSETXYL,SSETZL,
     & SMATXY,SMATZ,BB,NCONV,NCONV2,BX0L,R0XYL,R0ZL,DIFP)
      IMPLICIT NONE
      INTEGER MXE,JCUR,KVOL,MAXCUR,NEX0,NEX2,NEXX
      PARAMETER (MXE=708,JCUR=32,KVOL=32,MAXCUR=1024)
      PARAMETER (NEX0=164,NEX2=356,NEXX=420)
      INTEGER NEE,NUMPOL,NCONV,NCONV2
      INTEGER NCUR,NVV,N1,N2,NN,NUMPOLL,NCONVL,NCONV2L
      DOUBLE PRECISION SMATXY,SMATZ,BB, DIFP
      DOUBLE PRECISION SMATXYL,SMATZL,BBL,DIFPL
      DOUBLE PRECISION SSETXYL, SSETZL, BX0L
      DOUBLE PRECISION R0XYL,R0ZL
      DIMENSION NUMPOLL(KVOL,JCUR),NUMPOL(*)
      DIMENSION SSETXYL(MXE*MAXCUR)
      DIMENSION SSETZL(MXE*MAXCUR)
      DIMENSION SMATXY(*),SMATXYL(MXE*KVOL*JCUR)
      DIMENSION SMATZ(*),SMATZL(MXE*KVOL*JCUR)
      DIMENSION BB(*),BBL(NEXX*NEXX),NCONVL(2,MXE),NCONV(*)
      DIMENSION NCONV2L(6,64),NCONV2(*),BX0L(2*MXE)
      DIMENSION R0XYL(MXE*MXE),R0ZL(MXE*MXE)
      DIMENSION DIFP(*),DIFPL(KVOL*JCUR)
      INTEGER I,J,K,L,N
C----------------------------------------------------------------------
      DO J=1,JCUR
         DO I=1,KVOL
            NUMPOLL(I,J) = NUMPOL((J-1)*KVOL+I)
         END DO
      END DO
      DO J=1,MXE
         DO I=1,2
            NCONVL(I,J) = NCONV((J-1)*2+I)
         END DO
      END DO
      DO J=1,64
         DO I=1,6
            NCONV2L(I,J) = NCONV2((J-1)*6+I)
         END DO
      END DO
      DO I=1,MXE*KVOL*JCUR
         SMATXYL(I)=SMATXY(I)
      END DO
      DO I=1,MXE*KVOL*JCUR
         SMATZL(I)=SMATZ(I)
      END DO
      DO I=1,NEXX*NEXX
         BBL(I)=BB(I)
      END DO
      DO I=1,KVOL*JCUR
         DIFPL(I)=DIFP(I)
      END DO
C      DO I=1,2*MXE
C         BX0L(I)=BX0(I)
C      END DO
C      DO I=1,4*MXE*MXE
C         R0L(I) = R0(I)
C      END DO
C      DO I=1,MXE*MAXCUR
C         SSETXYL(I)=SSETXY(I)
C      END DO
C      DO I=1,MXE*MAXCUR
C         SSETZL(I) =SSETZ(I)
C      END DO
C----------- CALCUALTION  -------------------------------------
      DO K=1,NEE
         BX0L(K)=0.D0
         BX0L(MXE+K)=0.D0
         DO NCUR=1,JCUR
            DO NVV=1,KVOL
               L=NUMPOLL(NVV,NCUR)
               IF(L.GT.0) THEN
                  SSETXYL((L-1)*MXE+K)=0.D0
                  SSETZL((L-1)*MXE+K)=0.D0
                  DO N=1,NEX0
                     SSETXYL((L-1)*MXE+K)=SSETXYL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                    *SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N)
                     SSETZL((L-1)*MXE+K)=SSETZL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                    *SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N)
                  END DO
                  DO N=NEX0+1,NEX2
                     N1=NCONVL(1,N)
                     N2=NCONVL(2,N)
                     IF(N2.EQ.0) THEN
                        SSETXYL((L-1)*MXE+K)=SSETXYL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                 *SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
                        SSETZL((L-1)*MXE+K)=SSETZL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                 *SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
                     ELSE
                        SSETXYL((L-1)*MXE+K)=SSETXYL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &             *(SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
     &              +SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N2))
                        SSETZL((L-1)*MXE+K)=SSETZL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &             *(SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
     &              +SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N2))
                     ENDIF
                  END DO
                  DO N=NEX2+1,NEXX
                     NN=N-NEX2
                     SSETXYL((L-1)*MXE+K)=SSETXYL((L-1)*MXE+K)
     &                +BBL((K-1)*NEXX+N)*(
     &            SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(1,NN))
     &          +(SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(2,NN))
     &           +SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(3,NN))
     &           +SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(4,NN))
     &    +SMATXYL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(5,NN)))/2.D0)
                     SSETZL((L-1)*MXE+K)=SSETZL((L-1)*MXE+K)
     &                +BBL((K-1)*NEXX+N)*(
     &            SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(1,NN))
     &          +(SMATZL((NCUR-1)*MXE+KVOL+(NVV-1)*MXE+NCONV2L(2,NN))
     &           +SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(3,NN))
     &           +SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(4,NN))
     &    +SMATZL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(5,NN)))/2.D0)
                  END DO
                  BX0L(K)=BX0L(K)-SSETXYL((L-1)*MXE+K)
     &                 *DIFPL((NCUR-1)*KVOL+NVV)
                  BX0L(MXE+K)=BX0L(MXE+K)-SSETZL((L-1)*MXE+K)
     &                 *DIFPL((NCUR-1)*KVOL+NVV)
               ENDIF
            END DO
         END DO
         DO J=K,1,-1
            R0XYL((K-1)*MXE+J)=0.D0
            R0ZL((K-1)*MXE+J)=0.D0
            DO L=1,MAXCUR
               R0XYL((K-1)*MXE+J)=R0XYL((K-1)*MXE+J)
     &              +SSETXYL((L-1)*MXE+J)*SSETXYL((L-1)*MXE+K)
               R0ZL((K-1)*MXE+J)=R0ZL((K-1)*MXE+J)  
     &              +SSETZL((L-1)*MXE+J)*SSETZL((L-1)*MXE+K)
            END DO
         END DO
      END DO
C      DO I=1,2*MXE
C         BX0(I)=BX0L(I)
C      END DO
C      DO I=1,4*MXE*MXE
C         R0(I)=R0L(I)
C      END DO
C      DO I=1,MXE*MAXCUR
C         SSETXY(I)=SSETXYL(I)
C      END DO
C      DO I=1,MXE*MAXCUR
C         SSETZ(I)=SSETZL(I)
C      END DO
      RETURN
      END SUBROUTINE UPSENSEA
C============================================================================
