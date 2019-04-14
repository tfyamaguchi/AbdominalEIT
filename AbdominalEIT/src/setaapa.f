CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCSCCCC
C
C     ±ßÅû·Á¾õ¤Î3¼¡¸µFEM
C
C     X: ÅÅÎ®ÅÅ¶Ë¤ÎÊÒÂ¦ 2Í×ÁÇÃ±°Ì´ÖÅÅ°Ìº¹¤Ç¡¢¥Ç¡¼¥¿¿ô¡á32*29¡£
C     FA: ºÇ³°¼þ¤ÎÅÁÆ³ÅÙÊ¬ÉÛ¤ËÀ©¸Â¤òÀß¤±¤ë¡£NEWTON¤¬Í×ÁÇÊ¬³ä¤Ë°ÍÂ¸¤¹¤ë¤Î¤ÇÃí°Õ¡£
C        8.30: STEP UP ¤ò¸ÇÍ­ÃÍ¸Â³¦Êý¼°¤«¤é¸ÇÍ­¥Ù¥¯¥È¥ë¿ôÊý¼°¤ËÊÑ¤¨¤¿¡£CNEE¡£
C
C       MXE,        MZE,         MXN,        MZN,         MXB=128
C  ºÇÂç¡§ÌÌÆâÍ×ÁÇ¿ô¡¢ZÊý¸þÍ×ÁÇ¿ô¡¢ÌÌÆâÀáÅÀ¿ô¡¢ZÊý¸þÀáÅÀ¿ô¡¢¼þ°ÏÀáÅÀ¿ô
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     AA0
C
      SUBROUTINE SETAAPA(AA0,LAAL,SIGXY,LSIGXY,SIGZ,LSIGZ,NODEX,LNODEX,
     &     NODEZ,LNODEZ,ST,LST,STZ,LSTZ,NNODEZ,NBW,NEZ,NEX,MXNODE)
      IMPLICIT NONE
      DOUBLE PRECISION AA0,SIGXY,SIGZ,ST,STZ
      DOUBLE PRECISION SIGXYL,SIGZL,STL,STZL
      INTEGER LAAL,LSIGXY,LSIGZ,LNODEX,LNODEZ,LST,LSTZ
      INTEGER NNODEZ,NODEX,NODEZ,NEZ,NEX,NBW
      INTEGER NODEXL,NODEZL
      INTEGER MXNODE
      DIMENSION AA0(LAAL),
     &     SIGXY(LSIGXY),SIGXYL(LSIGXY),SIGZ(LSIGZ),SIGZL(LSIGZ),
     &     NODEX(LNODEX),NODEXL(LNODEX),
     &     NODEZ(LNODEZ),NODEZL(LNODEZ),
     &     ST(LST),STL(LST),
     &     STZ(LSTZ),STZL(LSTZ)
      INTEGER I,J,K,L,NX,NZ,JXZ,KXZ,JXZS
      DOUBLE PRECISION SIGMXY,SIGMZ
C        DO I=1,LAAL
C           AA0(I)=AAL(I)
C        END DO
        DO I=1,LSIGXY
           SIGXYL(I) = SIGXY(I)
        END DO
        DO I=1,LSIGZ
           SIGZL(I) = SIGZ(I)
        END DO
        DO I=1,LNODEX
           NODEXL(I)=NODEX(I)
        END DO
        DO I=1,LNODEZ
           NODEZL(I)=NODEZ(I)
        END DO
        DO I=1,LST
           STL(I)=ST(I)
        END DO
        DO I=1,LSTZ
           STZL(I)=STZ(I)
        END DO
C----------------------------------------------------------------------
        DO NZ=1,NEZ
           DO NX=1,NEX
              SIGMXY=SIGXYL(NX)
              SIGMZ =SIGZL(NX)
              DO J=1,4
                 JXZ=NNODEZ*(NODEXL((NX-1)*4*NEZ+(NZ-1)*4+J)-1)
     &                +NODEZL((NX-1)*4*NEZ+(NZ-1)*4+J)
                 DO K=1,4
                    KXZ=NNODEZ*(NODEXL((NX-1)*4*NEZ+(NZ-1)*4+K)-1)
     &                   +NODEZL((NX-1)*4*NEZ+(NZ-1)*4+K)
                    JXZS = JXZ - KXZ + 1
                    IF(JXZS.GT.NBW) THEN
                       STOP
                    ENDIF
                    IF (JXZS.GE.1) AA0((JXZS-1)*MXNODE+KXZ)
     &                  =AA0((JXZS-1)*MXNODE+KXZ)
     &                  +SIGMXY*STL((NX-1)*4*4*NEZ+(NZ-1)*4*4+(J-1)*4+K)
     &                  +SIGMZ*STZL((NX-1)*4*4*NEZ+(NZ-1)*4*4+(J-1)*4+K)
                 END DO
              END DO
           END DO
        END DO
C-------------------------------------------------------------------------
C        DO I=1,LAAL
C           AAL(I)=AA0(I)
C        END DO
        RETURN
        END SUBROUTINE
C===========================================================================
