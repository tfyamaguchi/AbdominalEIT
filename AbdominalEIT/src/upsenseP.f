C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE UPSENSEP (NEE,SSETL,SMAT,BB,NCONV,NCONV2,NEX0,
     & NEX2,NEXX,KVOL,CDIM)
      IMPLICIT NONE
      INTEGER MXE,JCUR,KVOL,MAXCUR,NEX0,NEX2,NEXX,CDIM
      PARAMETER (MXE=708,JCUR=32,MAXCUR=1024)
      INTEGER NEE,NCONV,NCONV2
      INTEGER NCUR,NVV,N1,N2,NN,NCONVL,NCONV2L
      DOUBLE PRECISION SMAT,SMATL,BB,BBL,DIFP,DIFPL
      DOUBLE PRECISION SSETL,BX0L
      DIMENSION SSETL(MXE*MAXCUR)
      DIMENSION SMAT(MXE*KVOL*JCUR),SMATL(MXE*KVOL*JCUR)
      DIMENSION BB(NEXX*NEXX),BBL(NEXX*NEXX),NCONVL(2,MXE),NCONV(2*MXE)
      DIMENSION NCONV2L(6,64),NCONV2(6*64)
      INTEGER I,J,K,L,N
C----------------------------------------------------------------------
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
         SMATL(I)=SMAT(I)
      END DO
      DO I=1,NEXX*NEXX
         BBL(I)=BB(I)
      END DO
C----------- CALCUALTION  -------------------------------------
      DO K=1,NEE
         DO NCUR=1,CDIM
            DO NVV=1,KVOL
               L=NVV + (NCUR-1)*KVOL
               IF(L.GT.0) THEN
                  SSETL((L-1)*MXE+K)=0.D0
                  DO N=1,NEX0
                     SSETL((L-1)*MXE+K)=SSETL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                    *SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N)
                  END DO
                  DO N=NEX0+1,NEX2
                     N1=NCONVL(1,N)
                     N2=NCONVL(2,N)
                     IF(N2.EQ.0) THEN
                        SSETL((L-1)*MXE+K)=SSETL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &                 *SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
                     ELSE
                        SSETL((L-1)*MXE+K)=SSETL((L-1)*MXE+K)
     &                    +BBL((K-1)*NEXX+N)
     &             *(SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N1)
     &              +SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+N2))
                     ENDIF
                  END DO
                  DO N=NEX2+1,NEXX
                     NN=N-NEX2
                     SSETL((L-1)*MXE+K)=SSETL((L-1)*MXE+K)
     &                +BBL((K-1)*NEXX+N)*(
     &            SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(1,NN))
     &          +(SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(2,NN))
     &           +SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(3,NN))
     &           +SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(4,NN))
     &    +SMATL((NCUR-1)*MXE*KVOL+(NVV-1)*MXE+NCONV2L(5,NN)))/2.D0)
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE UPSENSEP
C============================================================================
