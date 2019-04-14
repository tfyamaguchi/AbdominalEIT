CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ??????????Ǥ??оݤˤ????ե???E?????͡?
C
C     ?Ȥ??ǡ????ϳ???????????ɽ??4?��Ǥ???E??
C     ??????Ⱦ??1?αߤˤʡ?E褦?ˡ???Ķ?????????��?????��ֳ֤???¸?????褦?ˡ??
C     Ⱦ???????????��˼̡?E??Ƥ??????
C     ???????Ǥ?Ⱦ???????ˤʡ?E٤???ֳ֤ˤʤ??褦??Ⱦ??ʡ???��???��?��?��???????
C     ???δؿ??��? R??(STRT*R**3+R)/(1.+STRT) ?Ǥ???E??
C     ???ˤ???Ⱦ??1?αߤ?Ⱦ??1/(??*RATMAP)?ε??̾??˼ͱƤ???E??
C     ?ߤ??濴?ϵ??̤??̶ˤˡ??߼????̶ˤ??????̤˱??ä???Υ1?ε??̾??α߼??ˡ??ͱƤ???E????
C     ???ξ??ǡ??ǡ????��ޤǤε??̾???Υ??R?˱????ƽŤߤ??դ?????Ƴ?٤??п???ʿ?Ѥ???E??
C     ?Ťߴؿ???3??E??Ǥ??????
C     KFILT=0??Ⱦ?¦?R ??ERFILT ?ǽŤ?1??¾??0
C     KFILT=1??Ⱦ?¦?R ??ERFILT*??E ?? ?Ťߤ? 1-(??R/RFILT)**2/2??¾??0
C     KFILT=2???Ťߤ? exp{-(??R/RFILT)**2}
C     ??????E⡢?Ťߤ???ʡ?ͤ????RFILT**2??Ϳ???顦E????
C     ????E??Υѥ??᡼?????ILTER.PRN
C     ???ϥǡ????ե?????E??GMOUT0.PRN
C     ???ǥǡ????Ȥ???FEM5DATA.QQQ???????ǡ????Ȥ???BOUNDRY.QQQ??ɬ?פǤ???E??
C     Ʊ???��??ν??ϥե?????E??GMOUT.PRN?ȥ????ס?E??F?֡?EPRN
C     ?ʻҥѥ͡?E??ν??ϥǡ???????Ƴ?٤ξ????п??Ǥ??????
C     ?ʥե?????E??ϥ????ץ????G?֡?EDAT??
C     ???????ȹԤ˥ѥ??᡼???ȡ?
C     ???κ?ɸ?ϤǤΤ????��γ???????Ⱦ?¤?ñ?̤Ȥ???????E١ɱFILT?ˡ?E??ˤ????????Ƥ??????
C     ?濴????????Ū?Ǥ???E????????Ǥβ?��?٤?Ⱦ??ʡ?????⤡????ʡ?????㤤???
C2002.06.30 SIGB ?η׻??????ǿ??????ä??ɽ?????E褦?˽?��?
C2002.06.30 COMMENT ???Ϥ?Y?????ǡ????????ɲáʹԿ????ѡ?E餺?ˡ??
C2002.06.30 SUBROUTINE INPUT ?Ǥ?????????ARAS()?????̾??????ѤȤ?????
C2002.06.30 ʿ???ͤ??ե???E???��???Ǥɤ??Ѥ??ä????򼨤??????
C           0:????ʬ?ۡ?1:???????Ƿ��??Τޤޤǥե???E??????:��???ʻҤǤΥե???E????
C           ?ѡ?E?ʡ?ۤ??Ѱ??ˤ??ơ???????ʿ?Ѥ??ɤ???ʡ?ۤǤ???¸?????褦?ˤϽ????ʤ????
C           ???äơ??????̡?E2: ?򤽤Τޤ޽??Ϥ??Ƥ???E????: ?˹硦E??Ƥ??餹???Ȥ????褡?
C2002.06.30 ��???��?Υ?????????????????Ѥ??Ѳ?Ψ????��????E????Υ?????????��?򵭤??????
C
      SUBROUTINE FILTA06(FILENM,SIGMA,NEXD,NEZ,NODEXD,NODEZD,NNODEXD
     & ,NNODEZ,XCODD,YCODD,XYD,NBD,IBNDD
     & ,BOUND,IBOUND,NBUN,NBUNY,KFILT,RFILT,STRT,RATMAP,ZPLXY,IRET)
C      IMPLICIT REAL*8 ( A-H , O-Z )
      IMPLICIT NONE
      INTEGER MXE,MXN,MXB
      PARAMETER (MXE=708,MXN=419,MXB=128)
      DOUBLE PRECISION XCODD,YCODD,BOUND,RFILT,STRT,RATMAP,BND,ARAR
     & ,ARAS,ARSUM,AVESIGL0,AVESIGL1,AVESIGL2,AVSIG6,AVSIG7,BICHOSEI
     & ,DCS,DDX,DLG10,DRR,DSN,DTT,DXX,DYY,FRR,PI,R0MP,RAV,RMP,RRB
     & ,RRMIN,RRNG,RRNGLM,RRP,RRR,RST0,RST1,SCALE,SGMX,SGMXL,SIGB
     & ,SIGBB,SIGMA,SIGMAL,SIGMAV,SIGMLAV,SIGMLMN,SIGLMX,SIGMMN
     & ,SIGMMX,SIGMLMX,SIGOUT,SIZE,SLIM,STRC,SUMAREA,SUMAREA1,SUMR
     & ,SUMSIG,SUMSIGL0,SUMSIGL1,SUMSIGL2,SUMSIGT,TAKASA,TH0,THAV
     & ,THET,THETA,THETT,THST1,WEIGHT,WGHT,X0,XAV,XB,XCOD,XMAR,XMAX
     & ,XMIN,XXX,XY,Y0,YAV,YB,YCOD,YMAX,YMIN,YYY,ZPL,ZPLXY,XYD
      INTEGER NEXD,NEZ,NODEXD,NODEZD,NNODEXD,NNODEZ,NBD,IBNDD,NBUN
     & ,IBOUND,KFILT,IRET,HABA,I,J,K,L,LL,N,IBND,NB,NBUNY,NEL,NEX
     & ,NFO,NFOUND,NMIN,NN,NN11,NN22,NNODEX,NNSUM,NODEX,NODEZ,NX,NY
      COMMON NEX,NNODEX,NB
     & ,NODEX(3,MXE),NODEZ(3)
     & ,XCOD(MXN),YCOD(MXN),IBND(MXB+2)
     & ,XB(MXB+1),YB(MXB+1),SUMR(MXB+1)
     & ,ARAS(MXE),XAV(4,MXE),YAV(4,MXE)
     & ,RAV(4,MXE),THAV(4,MXE),THETA(MXB+1)
     & ,RRB(MXB+1),XMAX,XMIN,YMAX,YMIN,BND(264)
C      COMMON /MESH/F9M5(260520),IF9M5
      DIMENSION SIGMA(MXE),XY(6,MXE),SIGMAL(MXE),ARAR(MXE)
      DIMENSION ZPL(501),ZPLXY(501,501),BOUND(264)
      DIMENSION NODEXD(3,MXE),NODEZD(3),XCODD(MXN),YCODD(MXN)
     & ,IBNDD(MXB+2),XYD(6,MXE)
      CHARACTER*8 FILENM
      CHARACTER*8 WINDW(3)
      EXTERNAL STRC,WEIGHT
      DATA WINDW/' SQUARE ','PARABORA','GAUSSIAN'/
      NEX=NEXD
      NNODEX = NNODEXD
      NB = NBD
C      DO J=1,MXE
C       SIGMAL(J)=SIGMA(J)
C      END DO
      DO I=1,3
       DO J=1,MXE
         NODEX(I,J) = NODEXD(I,J)
       END DO
       NODEZ(I) = NODEZD(I)
      END DO
      DO I=1, MXN
        XCOD(I)=XCODD(I)
        YCOD(I)=YCODD(I)
      END DO
      DO I=1,MXB+2
       IBND(I)=IBNDD(I)
      END DO
      DO I=1, IBOUND
         BND(I)=BOUND(I)
      END DO
      PI=DATAN(1.D0)*4.D0
      SLIM=1.D-10
      DLG10=DLOG(1.D+01)
      OPEN(65,FILE='FILT06.LOG')
C      WRITE(65,*) ' FILTER'
C      WRITE(65,*) 'PREPARE SGMOUT0.PRN.'
C      WRITE(65,*) 'NOTE SGMOUT1.PRN BE OVERWRITTEN!'
C      OPEN(77,FILE='FILTER.PRN')
C      READ(77,*) NBUN
C      IF9M5=IF9M5D
C      DO I=1, 260520
C         F9M5(I)=F9M5D(I)
C      END DO
      IF(NBUN.LT.100) THEN
        WRITE(65,*) 'NBUN MUST BE LARGER THAN 99'
C        PAUSE
C        STOP
        IRET=2
        RETURN
      END IF
      IF(NBUN.GT.500) THEN
        WRITE(65,*) 'NBUN MUST BE LESS THAN 501'
C        PAUSE
C        STOP
        IRET=3
        RETURN
      END IF
C      READ(77,*) KFILT
      WRITE(65,'(A8,7H FILTER)') WINDW(KFILT+1)
C      READ(77,*) RFILT
      WRITE(65, '(7HRANGE= ,F10.5)') RFILT
      RRNGLM=RFILT**2
      IF(KFILT.EQ.1) RRNGLM=2.D0*RRNGLM
C      READ(77,*) STRT
      WRITE(65,'(5HSTRT=,F10.3)') STRT
C      READ(77,*) RATMAP
      WRITE(65,'(7HRATMAP=,F10.5,4H  PI)')  RATMAP
      RMP=RATMAP*PI
C      CLOSE(77)
      RST0=RFILT*(1.D0+STRT)
      RST1=RFILT*(1.D0+STRT)/(3.D0+STRT)
      IF(RMP.LT.SLIM) THEN
        THST1=RFILT
      ELSE IF(RMP.GT.(PI-SLIM)) THEN
        THST1=RFILT*(PI-SLIM)/DSIN(PI-SLIM)
      ELSE
        THST1=RFILT*RMP/DSIN(RMP)
      ENDIF
      WRITE(65,'(16HRANGE AT CENTER=,E15.5)') RST0
      WRITE(65,'(16H AT BOUNDRY DR=,E15.5,8H DTHETA=,E15.5)') RST1,THST1
      WRITE(65,'(10H NEX NEZ= ,2I10)') NEX,NEZ
      WRITE(65,'(10H NODEZ=   ,3I10)') NODEZ(1),NODEZ(2),NODEZ(3)
      DO I=1,MXE
       WRITE(65,'(3I8)') NODEX(1,I),NODEX(2,I),NODEX(3,I)
      END DO
      DO I=1,MXN
       WRITE(65,'(2E15.5)') XCOD(I),YCOD(I)
      END DO
      WRITE(65,'(20H NNODEX NNODEZ NB=  ,3I10)') NNODEX, NNODEZ, NB
      DO I=1,NB+2
         WRITE(65,*) IBND(I)
      END DO
C      OPEN(8,FILE='SGMOUT0.PRN')
C      READ(8,*) NEL
      NEL=NEX
      SGMXL=-1.D10
C      OPEN(88,FILE='XYTST.TXT')
      SUMAREA=0.D0
      SUMSIGL0=0.D0
      DO N=1,NEL
C        READ(8,*) SIGMA(N),(XY(J,N),J=1,6)
CCC SLIM
        DO J=1,6
         XY(J,N)=XYD(J,N)
        END DO
        SIGMAL(N)=DLOG(SIGMA(N))/DLG10
        IF(SGMXL.LT.SIGMAL(N)) SGMXL=SIGMAL(N)
        ARAR(N)=(XY(1,N)*XY(4,N)+XY(3,N)*XY(6,N)+XY(5,N)*XY(2,N)
     &          -XY(2,N)*XY(3,N)-XY(4,N)*XY(5,N)-XY(6,N)*XY(1,N))/2.D0
        SUMAREA=SUMAREA+ARAR(N)
        SUMSIGL0=SUMSIGL0+ARAR(N)*SIGMAL(N)
CCC NOTE
      END DO
      AVESIGL0=SUMSIGL0/SUMAREA
      SGMX=DEXP(SGMXL*DLG10)
      WRITE(65,*) 'SIGMA WAS READ'
      DO I=1,MXE
        WRITE(65,'(2E15.8)') SIGMA(I),SIGMAL(I)
      END DO
C      READ(8,*) FILENM
C      FILENM(7:7)='F'
C      CLOSE(8)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ???ޤǶ????Ǥ?ʿ????Ƴ?٤????ᡦE??
C     ????ʬ???˰?¸????E??
C     BACKGROUND LEVEL ?ˤ????ƶ��??ʤ???
C     BACKGROUND LEVEL??SIGOUT??SIGB-0.5*(SIGMAXL-SIGB)
C
      NN11=NEL-320
      NN22=NEL-192
      AVSIG6=0.D0
      DO N=NN11+1,NN22
        AVSIG6=AVSIG6+SIGMAL(N)
      END DO
      AVSIG6=AVSIG6/128.D0
      AVSIG7=0.D0
      DO N=1,64
        AVSIG7=AVSIG7+SIGMAL(3*N+NN22-1)*2.D0
     &               +SIGMAL(3*N+NN22-2)+SIGMAL(3*N+NN22)
      END DO
      AVSIG7=AVSIG7/256.D0
C      SIGB=(3.5D0*AVSIG7-AVSIG6)/2.5D0
      SIGB=(3.5D0*AVSIG7-AVSIG6)*10.0D0
C
CCCCCCCCCCCCCCCCCCCC
C
      CALL INPUTS(STRT,RMP,NEL)
C
CCCCCCCCCCCCCCCCCCCCCC
C
C
C      OPEN(8,FILE='SGMOUT1.PRN')
      SUMAREA1=0
      SUMSIGL1=0.D0
C      WRITE(8,*) NEX
      DO NN=1,NEX
        SUMSIGT=0.D0
        DO LL=1,4
          X0=XAV(LL,NN)
          Y0=YAV(LL,NN)
C              R0MP: SURFACE DISTANCE FROM NORTH POLE
          R0MP=RAV(LL,NN)*RMP
          TH0=THAV(LL,NN)
          IF(R0MP.GE.SLIM) THEN
            FRR=(DSIN(R0MP)/R0MP)**2
            DCS=DCOS(TH0)
            DSN=DSIN(TH0)
          ENDIF
          SUMSIG=0.D0
          ARSUM=0.D0
          DO N=1,NEL
            DO L=1,4
                DXX=(X0-XAV(L,N))
                DYY=(Y0-YAV(L,N))
                IF(R0MP.GE.SLIM) THEN
                  DRR=DXX*DCS+DYY*DSN
                  DTT=-DXX*DSN+DYY*DCS
                  RRNG=DRR**2+FRR*DTT**2
                ELSE
                  RRNG=DXX**2+DYY**2
                ENDIF
                WGHT=ARAS(N)*WEIGHT(RRNGLM,RRNG,KFILT)
                ARSUM=ARSUM+WGHT
                SUMSIG=SUMSIG+SIGMAL(N)*WGHT
            END DO
          END DO
          SUMSIGT=SUMSIGT+SUMSIG/ARSUM
        END DO
        SUMSIGT=SUMSIGT/4.D0
        SUMAREA1=SUMAREA1+ARAR(NN)
        SUMSIGL1=SUMSIGL1+SUMSIGT*ARAR(NN)
        SIGMA(NN)=DEXP(SUMSIGT*DLG10)
C        WRITE(8,'(7E15.5)') SIGMA(NN),(XY(J,NN),J=1,6)
      END DO
      AVESIGL1=SUMSIGL1/SUMAREA1
C      WRITE(8,'(1H ,A8)') FILENM
C      WRITE(8,'(E10.2)') SGMX
C      CLOSE(8)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      SIGOUT=SIGB-0.5D0*(SGMXL-SIGB)
C      SIGOUT=0.D0
C      HABA=384
      HABA=XMAX-XMIN
C      TAKASA=288
      TAKASA=YMAX-YMIN
C      XMAR=0.01D0*HABA
      XMAR=0
      SIZE=HABA+XMAR
      IF(HABA.LT.TAKASA) SIZE=TAKASA+XMAR
      DDX=SIZE/NBUN
      NBUNY=(TAKASA+XMAR)/DDX+1.D0
      DO J=NB/2,NB+1
        IF(THETA(J).LT.THETA(1)) THETA(J)=THETA(J)+2.D0*PI
      END DO
      THETA(NB+1)=THETA(NB+1)+SLIM
C      OPEN(77,FILE='THETEST.TXT')
      DO J=1,NB+1
C      WRITE(77,*) J, THETA(J)
      END DO
C      CLOSE(77)
C      FILENM(7:7)='G'
C      OPEN(77,FILE=FILENM//'.DAT')
C      WRITE(77,'(6H"NBUN=,2I10)') NBUN, NBUNY
C      WRITE(77,'(1H",A8,7H FILTER)') WINDW(KFILT+1)
C      WRITE(77,'(8H"RANGE= ,F10.5)') RFILT
C      WRITE(77,'(6H"STRT=,F10.3)') STRT
C      WRITE(77,'(8H"RATMAP=,F10.5,4H  PI)')  RATMAP
C      WRITE(77,'(17H"RANGE AT CENTER=,E15.5)') RST0
C      WRITE(77,'(17H"  AT BOUNDRY DR=,E15.5)') RST1
C      WRITE(77,'(17H"         DTHETA=,E15.5)') THST1
      SIGMLMX=-1.D+10
      SIGMLMN=1.D0+10
      SUMSIGL2=0.D0
      NNSUM=0
      DO NY=1,NBUNY
        YYY=YMAX+0.5D0*XMAR-(NY-1)*DDX
        DO NX=1,NBUN
C     CC NBUN+1
          XXX=XMIN-0.5D0*XMAR+(NX-1)*DDX
          THET=DATAN2(YYY,XXX)
          IF(THET.LT.THETA(1)) THET=THET+2.D0*PI
          RRP=DSQRT(XXX**2+YYY**2)
          NFOUND=0
          DO J=1,NB
            IF((THET.GE.THETA(J)).AND.(THET.LT.THETA(J+1))) THEN
              NFOUND=1
              THETT
     &        =((THET-THETA(J))*SUMR(J+1)+(THETA(J+1)-THET)*SUMR(J))
     &         /(THETA(J+1)-THETA(J))
              THETT=THETT-PI*0.5D0
              RRR
     &        =((THET-THETA(J))*RRB(J+1)+(THETA(J+1)-THET)*RRB(J))
     &         /(THETA(J+1)-THETA(J))
            ENDIF
          END DO
          IF(NFOUND.EQ.0) THEN
            WRITE(65,*) 'NOT FOUND',NY,NX,THET
C            PAUSE
C            STOP
            IRET=4
            RETURN
          ENDIF
          IF(RRP.GT.RRR) THEN
            ZPLXY(NX,NY)=SIGOUT
            X0=RRP*DCOS(THETT)/RRR
            Y0=RRP*DSIN(THETT)/RRR
          ELSE
            RRP=STRC(RRP/RRR,STRT)
            X0=RRP*DCOS(THETT)
            Y0=RRP*DSIN(THETT)
C              R0MP: SURFACE DISTANCE FROM NORTH POLE
            R0MP=RRP*RMP
            IF(R0MP.GE.SLIM) THEN
              FRR=(DSIN(R0MP)/R0MP)**2
              DCS=DCOS(THETT)
              DSN=DSIN(THETT)
            ENDIF
            ZPLXY(NX,NY)=0.D0
            SUMSIG=0.D0
            ARSUM=0.D0
            RRMIN=1.D10
            NFO=0
            NMIN=0
            DO N=1,NEL
              DO L=1,4
                DXX=X0-XAV(L,N)
                DYY=Y0-YAV(L,N)
                IF(R0MP.GE.SLIM) THEN
                  DRR=DXX*DCS+DYY*DSN
                  DTT=-DXX*DSN+DYY*DCS
                  RRNG=DRR**2+FRR*DTT**2
                ELSE
                  RRNG=DXX**2+DYY**2
                ENDIF
                IF(RRMIN.GT.RRNG) THEN
                  RRMIN=RRNG
                  NMIN=N
                ENDIF
                  WGHT=ARAS(N)*WEIGHT(RRNGLM,RRNG,KFILT)
                  ARSUM=ARSUM+WGHT
                  SUMSIG=SUMSIG+SIGMAL(N)*WGHT
              END DO
            END DO
            IF(ARSUM.GT.SLIM) THEN
              ZPLXY(NX,NY)=SUMSIG/ARSUM
            ELSE
              ZPLXY(NX,NY)=SIGMAL(NMIN)
            ENDIF
            NNSUM=NNSUM+1
            SUMSIGL2=SUMSIGL2+ZPLXY(NX,NY)
            IF(SIGMLMX.LT.ZPL(NX)) SIGMLMX=ZPLXY(NX,NY)
            IF(SIGMLMN.GT.ZPL(NX)) SIGMLMN=ZPLXY(NX,NY)
          ENDIF
        END DO
      END DO
      AVESIGL2=SUMSIGL2/NNSUM
      SIGMLAV=AVESIGL2
C
C     ????ʬ?ۤ?ʿ???ͤ??硦E????٤ˤϰʲ???��???????
C
C      SIGMLAV=AVESIGL0
C      SIGMLMX=SIGMLMX+AVESIGL0-AVESIGL2
C      SIGMLMN=SIGMLMN+AVESIGL0-AVESIGL2
C      SIGOUT=SIGOUT+AVESIGL0-AVESIGL2
C      DO NY=1,NBUNY
C      DO NX=1,NBUN+1
C        ZPLXY(NX,NY)=ZPLXY(NX,NY)+AVESIGL0-AVESIGL2
C      END DO
C      END DO
C
      DO NY=1,NBUNY
C        WRITE(77,'(E15.5,501(1H,,E15.5))')(ZPLXY(NX,NY),NX=1,NBUN)
C     CC NBUN+1
      END DO
      DO J=1,2
        DO K=1,NBUN
C     CC NBUN+1
          ZPL(K)=SIGOUT
        END DO
C        WRITE(77,'(E15.5,501(1H,,E15.5))') (ZPL(NX),NX=1,NBUN)
C     CC NBUN+1
      END DO
      DO J=1,2
        DO K=26,75
          ZPL(K)=SIGMLMX
        END DO
C        WRITE(77,'(E15.5,501(1H,,E15.5))') (ZPL(NX),NX=1,NBUN)
C     CC NBUN+1
      END DO
      DO J=1,2
        DO K=1,NBUN
C     CC
          ZPL(K)=SIGOUT
        END DO
C        WRITE(77,'(E15.5,501(1H,,E15.5))') (ZPL(NX),NX=1,NBUN)
C     CC
      END DO
      SCALE=50.D0*DDX
      BICHOSEI=DSQRT(SUMAREA/NNSUM)/DDX
      SCALE=50.D0*DDX*BICHOSEI
      SIGMAV=DEXP(SIGMLAV*DLG10)
      SIGMMX=DEXP(SIGMLMX*DLG10)
      SIGMMN=DEXP(SIGMLMN*DLG10)
      SIGBB=DEXP(SIGB*DLG10)
C      WRITE(77,'(17H"LOGAVE SIGMAL  =,E15.5)') SIGMLAV
C      WRITE(77,'(17H"MAX SIGMAL     =,E15.5)') SIGMLMX
C      WRITE(77,'(17H"MIN SIGMAL     =,E15.5)') SIGMLMN
C      WRITE(77,'(17H"BOUNDRY SIGMAL =,E15.5)') SIGB
C      WRITE(77,'(17H"LOGAVE SIGMA   =,E15.5)') SIGMAV
C      WRITE(77,'(17H"MAX SIGMA      =,E15.5)') SIGMMX
C      WRITE(77,'(17H"MIN SIGMA      =,E15.5)') SIGMMN
C      WRITE(77,'(17H"BOUNDRY SIGMA  =,E15.5)') SIGBB
C      WRITE(77,'(17H"SCALE          =,E15.5)') SCALE
C      WRITE(77,'(17H"LGV SIGMAL0    =,E15.5)') AVESIGL0
C      WRITE(77,'(17H"LGV SIGMAL1    =,E15.5)') AVESIGL1
C      WRITE(77,'(17H"LGV SIGMAL2    =,E15.5)') AVESIGL2
C      WRITE(77,'(17H"BICHOSEI       =,E15.5)') BICHOSEI
C      CLOSE(77)
      CLOSE(65)
C     STOP
      IRET=0
      END SUBROUTINE FILTA06
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION STRC(R,STRT)
      IMPLICIT REAL*8 ( A-H , O-Z )
        STRC=(STRT*R*R*R+R)/(STRT+1.D0)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION WEIGHT(R0,R,KF)
      IMPLICIT REAL*8 ( A-H , O-Z )
                IF(KF.GE.2) THEN
                  WEIGHT=DEXP(-R/R0)
                  RETURN
                ELSE IF(R.LT.R0) THEN
                  IF(KF.EQ.1) THEN
                    WEIGHT=(R0-R)/R0
                    RETURN
                  ELSE IF(KF.LE.0) THEN
                    WEIGHT=1.D0
                    RETURN
                  ENDIF
                ELSE
                  WEIGHT=0.D0
                  RETURN
                ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE INPUTS(STRT,RMP,NEL)
      IMPLICIT REAL*8 ( A-H , O-Z )
      PARAMETER (MXE=708,MXN=419,MXB=128)
      COMMON NEX,NNODEX,NB
     &  ,NODEX(3,MXE),NODEZ(3)
     &  ,XCOD(MXN),YCOD(MXN),IBND(MXB+2)
     &  ,XB(MXB+1),YB(MXB+1),SUMR(MXB+1)
     &  ,ARAS(MXE),XAV(4,MXE),YAV(4,MXE)
     &  ,RAV(4,MXE),THAV(4,MXE),THETA(MXB+1)
     &  ,RRB(MXB+1),XMAX,XMIN,YMAX,YMIN,BND(264)
C      COMMON /MESH/F9M5(260520),IF9M5
      DIMENSION X(3),Y(3)
      EXTERNAL STRC
      SLIM=1.D-10
      PI=4.D0*DATAN(1.D0)
	IF(NEL.EQ.516) THEN
C        OPEN (8,FILE ='FEM5DATA.QQQ')
	ELSE IF(NEL.EQ.644) THEN
C	  OPEN(8,FILE='F8M5DATA.QQQ')
	ELSE IF(NEL.EQ.708) THEN
C	  OPEN(8,FILE='F9M5DATA.QQQ')
	ELSE
	  WRITE(65,*) 'BUNKATUDATA NOT AVAILABLE',NEL
C	  PAUSE
C	  STOP
          IRET=5
          RETURN
	ENDIF
C      READ (8,*) NEX,NEZ
C      NEX=INT(F9M5(1))
C      NEZ=INT(F9M5(2))
      L=1
      IF(NEX.NE.NEL) THEN
	  WRITE(65,*) 'INCONSISTENT BUNKATUDATA,NEX=',NEL,NEX
C	  PAUSE
C	  STOP
	IRET=6
        RETURN
       ENDIF
      DO K=1,NEX
      DO J=1,NEZ
        IF(J.EQ.1) THEN
C          READ (8,*) IEX,IEZ,(NODEX(I,K),NODEZ(I),I=1,3)
C          IEX=INT(F9M5(L*10+1))
C          IEZ=INT(F9M5(L*10+2))
          DO I=1,3
C            NODEX(I,K)=INT(F9M5(L*10+1+2*I))
C            NODEZ(I)=INT(F9M5(L*10+2+2*I))
          END DO
        ELSE
C          READ(8,*) IDUM
C           IDUM=INT(F9M5(L*10+1))
        ENDIF
        L=L+1
      END DO
      END DO
          WRITE(65,*) 'START READ'
C      READ (8,*) NNODEX,NNODEZ
C      NNODEX=INT(F9M5(L*10+1))
C      NNODEZ=INT(F9M5(L*10+2))
      L=L+1
      DO J=1,NNODEX
C        READ(8,*) NODE,XCOD(J),YCOD(J)
C        NODE=INT(F9M5(L*10+1))
C        XCOD(J)=F9M5(L*10+2)
C        YCOD(J)=F9M5(L*10+3)
        L=L+1
        RTT=DSQRT(XCOD(J)**2+YCOD(J)**2)
        IF(RTT.GT.SLIM) THEN
          RTT=STRC(RTT,STRT)/RTT
          XCOD(J)=RTT*XCOD(J)
          YCOD(J)=RTT*YCOD(J)
        ENDIF
C        IF(J.NE.NODE) THEN
C          WRITE(65,*) 'X-NODE NUMBER INCONSISTENT'
C         PAUSE
C          STOP
C          IRET=7
C          RETURN
C        ENDIF
      END DO
          WRITE(65,*) 'START READ2'
      DO J=1,NNODEZ
C        READ(8,*) NODE
C        NODE=INT(F9M5(L*10+1))
        L=L+1
        IF(J.NE.NODE) THEN
          WRITE(65,*) 'Z-NODE NUMBER INCONSISTENT'
C          PAUSE
C          STOP
          IRET=8
          RETURN
        ENDIF
      END DO
C      READ(8,*) NB
C      NB=INT(F9M5(L*10+1))
      L=L+1
      DO J=1,NB
C        READ(8,*) IBND(J)
C        IBND(J)=INT(F9M5(L*10+1))
        L=L+1
      END DO
      IBND(NB+1)=IBND(1)
          WRITE(65,*) 'START READ3'
      CLOSE(8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      OPEN(8,FILE='BOUNDRY.QQQ',FORM='FORMATTED')
C      READ(8,*) NBOUND
      NBOUND=INT(BND(1))
      IF(NB.NE.2*NBOUND) THEN
        WRITE(6,*) 'INCONSISTENT BOUNDARY SIZE',NB,NBOUND
C        PAUSE
C        STOP
        IRET=9
        RETURN
      ENDIF
      SUMXX=0.D0
      SUMYY=0.D0
      XMAX=0.D0
      XMIN=0.D0
      YMAX=0.D0
      YMIN=0.D0
      DO J=1,NBOUND
C        READ(8,*) NUMB,XB(2*J-1),YB(2*J-1)
        NUMB=BND(J*4+1)
        XB(2*J-1)=BND(J*4+2)
        YB(2*J-1)=BND(J*4+3)
        IF(XMAX.LT.XB(2*J-1)) XMAX=XB(2*J-1)
        IF(XMIN.GT.XB(2*J-1)) XMIN=XB(2*J-1)
        IF(YMAX.LT.YB(2*J-1)) YMAX=YB(2*J-1)
        IF(YMIN.GT.YB(2*J-1)) YMIN=YB(2*J-1)
        SUMXX=SUMXX+XB(2*J-1)
        SUMYY=SUMYY+YB(2*J-1)
      END DO
C      CLOSE(8)
C
C
      SUMXX=SUMXX/NBOUND
      SUMYY=SUMYY/NBOUND
      DO J=1,NBOUND
        XB(2*J-1)=XB(2*J-1)-SUMXX
        YB(2*J-1)=YB(2*J-1)-SUMYY
      END DO
        XB(2*NBOUND+1)=XB(1)
        YB(2*NBOUND+1)=YB(1)
        XMAX=XMAX-SUMXX
        XMIN=XMIN-SUMXX
        YMAX=YMAX-SUMYY
        YMIN=YMIN-SUMYY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ?¸??ǡ????????֤?Ϳ???顦E??????��ϴ񿡦ֹ????б????????
C     ?????ֹ??ˤĤ??Ƥϡ????��??Ρ?E??
C
      DO J=1,NBOUND
        XB(2*J)=(XB(2*J-1)+XB(2*J+1))/2.D0
        YB(2*J)=(YB(2*J-1)+YB(2*J+1))/2.D0
      END DO
        SUMR(1)=0.D0
C
      DO J=1,NBOUND
          RR=(XB(2*J+1)-XB(2*J-1))**2+(YB(2*J+1)-YB(2*J-1))**2
          RR=DSQRT(RR)
          SUMR(2*J+1)=SUMR(2*J-1)+RR
      END DO
        DSIZE=SUMR(2*NBOUND+1)/NB
        F360=2.D0*PI/SUMR(2*NBOUND+1)
        C360=3.6D+02/SUMR(2*NBOUND+1)
      DO J=1,NBOUND
        SUMR(2*J)=(SUMR(2*J-1)+SUMR(2*J+1))/2.D0
      END DO
      DO J=1,2*NBOUND+1
        SUMR(J)=SUMR(J)*F360
        THETA(J)=DATAN2(YB(J),XB(J))
        RRB(J)=DSQRT(XB(J)**2+YB(J)**2)
      END DO
      THETA(2*NBOUND+1)=THETA(1)+2.D0*PI
CCCCCCCCCCCC  SUMR??????Ĺ??????ɽ?????Ѥ???????  CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     REAL ELEMENT AREA
C
      SUMAREA=0.D0
      DO I=1,NEX
        DO L=1,3
          X(L)=XCOD(NODEX(L,I))
          Y(L)=YCOD(NODEX(L,I))
        END DO
        XAV(1,I)=(X(1)+X(2)+X(3))/3.D0
        YAV(1,I)=(Y(1)+Y(2)+Y(3))/3.D0
        RAV(1,I)=DSQRT(XAV(1,I)**2+YAV(1,I)**2)
        THAV(1,I)=DATAN2(YAV(1,I),XAV(1,I))
        R0MP=RAV(1,I)*RMP
        FRR=DSIN(R0MP)/R0MP
       area= X(1)*Y(2)-X(2)*Y(1)+X(2)*Y(3)-X(3)*Y(2)+X(3)*Y(1)-X(1)*Y(3)
        ARAS(I)=FRR*area/2.D0
        SUMAREA=SUMAREA+ARAS(I)
        DO L=1,3
          XAV(L+1,I)=(X(L)+XAV(1,I))/2.D0
          YAV(L+1,I)=(Y(L)+YAV(1,I))/2.D0
          RAV(L+1,I)=DSQRT(XAV(L+1,I)**2+YAV(L+1,I)**2)
          THAV(L+1,I)=DATAN2(YAV(L+1,I),XAV(L+1,I))
        END DO
      END DO
      WRITE(65,*) 'SUMAREA',SUMAREA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

