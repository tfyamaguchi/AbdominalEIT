CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     3D-FEM ELEMENTS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE B9C3D128(RADIUS,XXL,YYL,NEX,NEZ,NX2,NZ,NNODEX
     & ,NNODEZ,NB,IBND,R1,R2,R3,R4)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      DOUBLE PRECISION RADIUS,AREA,D1,D2,D3,D4,DDDZ,VOL,X,X1,X2,X3,XX
     & ,Y,Y1,Y2,Y3,YY,Z,ZERO,ZZ,XXL,YYL
      DOUBLE PRECISION R1,R2,R3,R4,RADAVE
      DOUBLE PRECISION PI,PIB3,PIB6,PIB8,PIB16,PIB32
      INTEGER MXN,I,ICCHI,II,J,JJ,JN,JR,K,KMAX,KMIN,KN,KWD,L,LIMIT
     & ,LL,LLL,m,MXB,MXE,MXZ,N,N1,N2,N3,NCN,NCN0,NCOM12,NCOM13,NCOM23
     & ,NCOR,NEX,NEXT,NEZ,NHF,NHWD,NINI,NLAST,NLIM,NNAIEND,NNODEX,NNODEZ
     & ,NR1,NR2,NREN,NRENREV,NREV,NSUM,NSUMX,NTYPE,NX,NZ,NX2,NB,IBND
      DIMENSION RADIUS(65)
      PARAMETER (MXN=419,MXE=708,MXB=128,MXZ=36)
C     NODE IS SPECIFIED IN N(XY) & N(Z)
C     ELEMENT IS SPECIFIED IN NODES 1-4, NX(XY) & NY(Z)
C               NODE IN X-Y  NODE IN Z   NODE NUMBERS FOR ELEMENT
      DIMENSION XX(MXN),YY(MXN),ZZ(MXZ),NX(4,MXE,MXZ),NZ(4,MXE,MXZ)
     & ,VOL(3),X(4),Y(4),Z(4),NCOR(MXN,MXN),NREN(MXN),NX2(4,MXE,MXZ)
     & ,NRENREV(MXN),NREV(4),NTYPE(MXE),XXL(MXN),YYL(MXN),IBND(MXB+2)
      COMMON NCN(10,MXN),NCN0(10,MXN)
C      PI=DACOS(-1.D0)
      PI = 3.141592653589793238D0
      ZERO=0.D0
      NB=MXB
      NNODEX=419
      NEX=708
      NNAIEND=291
C      OPEN(8,FILE='HR9NODE.PRN')
C      READ(8,*) R1,R2,R3,R4,R5,R6,R7,R8,R9
C       R1=1.75D-1
C       R2=2.95D-1
C       R3=4.0D-1
C       R5=7.6D-1
C       R5=8.0D-1
C       R6=8.5D-1
C       R6=8.6D-1
C       R7=9.2D-1
C       R8=9.7D-1
C       R9=1.0D0
C      DO I=1,64
C       RADIUS(I) = 1.0D0
C      END DO
C      READ(8,*) NNODEZ
      NNODEZ=13
      IF(NNODEZ.LT.2) THEN
C        WRITE(6,*) 'NNODEZ SHOULD BE >1'
C        PAUSE
        STOP
      ENDIF
C      CLOSE(8)
      NEZ=3*(NNODEZ-1)
CCCCCCCCCCCCCCCCCCCCC
C
C     NODE COODINATES ARE DEFINED HERE IN XY & IN Z INDEPENDENTLY
C
      RADAVE=0.0D0
      DO J=1,64
       RADAVE=RADAVE+RADIUS(J)
      END DO
      RADAVE=RADAVE/6.4D1
C      
      XX(1)=ZERO
      YY(1)=ZERO
      PIB3=PI/3.D0
C  1ST
      DO J=2,7
       XX(J)= R1*DSIN(PIB3*(J-2))
       YY(J)=-R1*DCOS(PIB3*(J-2))
      END DO
      PIB6=PI/6.D0
C  2ND
      DO J=8,19
       XX(J)= R2*DSIN(PIB6*(J-8))
       YY(J)=-R2*DCOS(PIB6*(J-8))
      END DO   
      PIB8=PI/8.D0
C  3RD
      DO J=20,35
C        XX(J)= R3*DSIN(PIB8*(J-20))
C        YY(J)=-R3*DCOS(PIB8*(J-20))
        XX(J)= R3/RADAVE*RADIUS(4*(J-19)-3)*DSIN(PIB8*(J-20))
        YY(J)=-R3/RADAVE*RADIUS(4*(J-19)-3)*DCOS(PIB8*(J-20))
      END DO
      PIB16=PI/16.D0
C  4TH & 5TH
      DO J=36,67
C        XX(J)= R4*DSIN(PIB16*(J-36))
C        YY(J)=-R4*DCOS(PIB16*(J-36))
        XX(J)= R4/RADAVE*RADIUS(2*(J-35)-1)*DSIN(PIB16*(J-36))
        YY(J)=-R4/RADAVE*RADIUS(2*(J-35)-1)*DCOS(PIB16*(J-36))
C        XX(J+32)= R5*DSIN(PIB16*(J-36))
C        YY(J+32)=-R5*DCOS(PIB16*(J-36))
C        XX(J+32)= R5/RADAVE*RADIUS(2*(J-35)-1)*DSIN(PIB16*(J-36))
C        YY(J+32)=-R5/RADAVE*RADIUS(2*(J-35)-1)*DCOS(PIB16*(J-36))
        XX(J+32)= (2.5D-1+RADIUS(2*(J-35)-1)/2.0D0)*DSIN(PIB16*(J-36))
        YY(J+32)=-(2.5D-1+RADIUS(2*(J-35)-1)/2.0D0)*DCOS(PIB16*(J-36))
      END DO
      PIB32=PI/32.D0
C  6TH & 7TH & 8TH & 9TH
      DO J=100,163
        XX(J)= RADIUS(J-99)*DSIN(PIB32*(J-100))
        YY(J)=-RADIUS(J-99)*DCOS(PIB32*(J-100))
        XX(J+64)= (1.0D0+RADIUS(J-99))/2.0D0*DSIN(PIB32*(J-100))
        YY(J+64)=-(1.0D0+RADIUS(J-99))/2.0D0*DCOS(PIB32*(J-100))
        XX(J+128)= (8.0D-1+2.0D-1*RADIUS(J-99))*DSIN(PIB32*(J-100))
        YY(J+128)=-(8.0D-1+2.0D-1*RADIUS(J-99))*DCOS(PIB32*(J-100))
        XX(2*J+92)= DSIN(PIB32*(J-100))
        YY(2*J+92)=-DCOS(PIB32*(J-100))
      END DO
      DO J=1,63
        XX(2*J+291)=(XX(2*J+290)+XX(2*J+292))/2.D0
        YY(2*J+291)=(YY(2*J+290)+YY(2*J+292))/2.D0
      END DO
        XX(419)=(XX(418)+XX(292))/2.D0
        YY(419)=(YY(418)+YY(292))/2.D0
CCCCCCCCCCCCCCCC
C      OPEN(8,FILE='TEST.TXT')
      DO J=1,NNODEX
C        WRITE(8,'(2F15.3)') XX(J),YY(J)
      END DO
CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC
C      CLOSE(8)
      DDDZ=1.D0/(NNODEZ-1)
      DO L=1,NNODEZ
        ZZ(L)=DDDZ*(L-1)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ELEMENTS IN WHICH THREE NODES ARE IN X=0,Y=0 ARE DEFINED HERE
C
      DO J=1,NEX
C         NZ(4,J,4)=2
C         NZ(4,J,6)=1
         NZ(4,J,4)=3
         NZ(4,J,6)=2
         DO K=1,3
C            NZ(K,J,4)=1
C            NZ(K,J,6)=2
            NZ(K,J,4)=2
            NZ(K,J,6)=3
         END DO
      END DO
CCCCCCCCCCCCCCCCCCCCCCC
C
C     INNER 1ST 1-6
C
      DO J=1,6
         NX(2,J,4)=1
         NX(3,J,4)=J+1
         NX(1,J,4)=J+2
      END DO
      NX(1,6,4)=2
CCCCCCCCCCCCCCCCCCCCCCC
C
C     INNER 2ND 7-24
C
      DO J=1,6
         NX(1,3*J+5,4)=J+1
         NX(2,3*J+5,4)=2*J+7
         NX(3,3*J+5,4)=J+2
         NX(2,3*J+4,4)=J+1
         NX(3,3*J+4,4)=2*J+6
         NX(1,3*J+4,4)=2*J+7
         NX(2,3*J+6,4)=J+2
         NX(3,3*J+6,4)=2*J+7
         NX(1,3*J+6,4)=2*J+8
      END DO
      NX(3,23,4)=2
      NX(2,24,4)=2
      NX(1,24,4)=8
C
CCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO L=1,4
CCCCCCCCCCCCCCCCCCCCCCC
C
C     3RD 25-52
C
         DO J=1,4
            NX(2,2*J+7*L+16,4)=J+3*L+4
            NX(3,2*J+7*L+16,4)=J+4*L+15
            NX(1,2*J+7*L+16,4)=J+4*L+16
         END DO
         DO J=1,3
            NX(1,2*J+7*L+17,4)=J+3*L+4
            NX(2,2*J+7*L+17,4)=J+4*L+16
            NX(3,2*J+7*L+17,4)=J+3*L+5
      END DO
      NX(2,52,4)=8
      NX(1,52,4)=20
      NX(3,51,4)=8
CCCCCCCCCCCCCCCCCCCCCCC
C
C     4TH 53-100
C
      DO J=1,4
         NX(2,3*J+12*L+38,4)=J+4*L+15
         NX(3,3*J+12*L+38,4)=2*J+8*L+26
         NX(1,3*J+12*L+38,4)=2*J+8*L+27
         NX(1,3*J+12*L+39,4)=J+4*L+15
         NX(2,3*J+12*L+39,4)=2*J+8*L+27
         NX(3,3*J+12*L+39,4)=J+4*L+16
         NX(2,3*J+12*L+40,4)=J+4*L+16
         NX(3,3*J+12*L+40,4)=2*J+8*L+27
         NX(1,3*J+12*L+40,4)=2*J+8*L+28
      END DO
      NX(3,99,4)=20
      NX(2,100,4)=20
      NX(1,100,4)=36
CCCCCCCCCCCCCCCCCCCCCCC
C
C     5TH 101-164
C
      DO J=1,4
         NX(2,4*J+16*L+81,4)=2*J+8*L+26
         NX(3,4*J+16*L+81,4)=2*J+8*L+58
         NX(1,4*J+16*L+81,4)=2*J+8*L+59
         NX(2,4*J+16*L+84,4)=2*J+8*L+28
         NX(3,4*J+16*L+84,4)=2*J+8*L+59
         NX(1,4*J+16*L+84,4)=2*J+8*L+60
         NX(1,4*J+16*L+82,4)=2*J+8*L+26
         NX(2,4*J+16*L+82,4)=2*J+8*L+59
         NX(3,4*J+16*L+82,4)=2*J+8*L+27
         NX(1,4*J+16*L+83,4)=2*J+8*L+27
         NX(2,4*J+16*L+83,4)=2*J+8*L+59
         NX(3,4*J+16*L+83,4)=2*J+8*L+28
      END DO
      NX(3,163,4)=36
      NX(2,164,4)=36
      NX(1,164,4)=68
C
CCCCCCCCCCCCCCCCCCCCCCC
C
C     6TH 165-260
C
      DO J=1,8
         NX(2,3*J+24*L+138,4)=J+8*L+59
         NX(3,3*J+24*L+138,4)=2*J+16*L+82
         NX(1,3*J+24*L+138,4)=2*J+16*L+83
         NX(1,3*J+24*L+139,4)=J+8*L+59
         NX(2,3*J+24*L+139,4)=2*J+16*L+83
         NX(3,3*J+24*L+139,4)=J+8*L+60
         NX(2,3*J+24*L+140,4)=J+8*L+60
         NX(3,3*J+24*L+140,4)=2*J+16*L+83
         NX(1,3*J+24*L+140,4)=2*J+16*L+84
      END DO
      NX(3,259,4)=68
      NX(2,260,4)=68
      NX(1,260,4)=100
CCCCCCCCCCCCCCCCCCCCCCC
C
C     7TH 261-388
C
      DO J=1,8
         NX(2,4*J+32*L+225,4)=2*J+16*L+82
         NX(3,4*J+32*L+225,4)=2*J+16*L+146
         NX(1,4*J+32*L+225,4)=2*J+16*L+147
         NX(2,4*J+32*L+228,4)=2*J+16*L+84
         NX(3,4*J+32*L+228,4)=2*J+16*L+147
         NX(1,4*J+32*L+228,4)=2*J+16*L+148
         NX(1,4*J+32*L+226,4)=2*J+16*L+82
         NX(2,4*J+32*L+226,4)=2*J+16*L+147
         NX(3,4*J+32*L+226,4)=2*J+16*L+83
         NX(1,4*J+32*L+227,4)=2*J+16*L+83
         NX(2,4*J+32*L+227,4)=2*J+16*L+147
         NX(3,4*J+32*L+227,4)=2*J+16*L+84
      END DO
      NX(3,387,4)=100
      NX(2,388,4)=100
      NX(1,388,4)=164
C
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     8TH 389-516
C
      DO J=1,64
         DO M=1,3
            NX(M,2*J+387,4)=NX(M,2*J+259,4)+64
            NX(M,2*J+388,4)=NX(M,2*J+260,4)+64
        END DO
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     9TH 517-708
C
      DO J=1,64
         NX(2,3*J+514,4)=J+227
         NX(3,3*J+514,4)=2*J+290
         NX(1,3*J+514,4)=2*J+291
         NX(1,3*J+515,4)=J+227
         NX(2,3*J+515,4)=2*J+291
         NX(3,3*J+515,4)=J+228
         NX(2,3*J+516,4)=J+228
         NX(3,3*J+516,4)=2*J+291
         NX(1,3*J+516,4)=2*J+292
      END DO
      NX(3,707,4)=228
      NX(2,708,4)=228
      NX(1,708,4)=292
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ELEMENT AREA CHECK
C
C      OPEN(8,FILE='TEST2.TXT')
      DO N=1,NEX
         N1=NX(1,N,4)
         N2=NX(2,N,4)
         N3=NX(3,N,4)
         X1=XX(N1)
         X2=XX(N2)
         X3=XX(N3)
         Y1=YY(N1)
         Y2=YY(N2)
         Y3=YY(N3)
         AREA=X1*Y2+X2*Y3+X3*Y1-X2*Y1-X3*Y2-X1*Y3
C         WRITE(8,'(4I5,7F15.5)') N,N1,N2,N3,X1,Y1,X2,Y2,X3,Y3,AREA
      END DO
C      CLOSE(8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE OF 1ST 1-6
C
      NTYPE(1)=4
      NTYPE(2)=1
      NTYPE(3)=1
      DO N=1,3
        NTYPE(7-N)=5-NTYPE(N)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE OF 2ND 7-24
C
      NTYPE(7)=4
      NTYPE(8)=2
      NTYPE(9)=1
      NTYPE(10)=1
      NTYPE(11)=3
      NTYPE(12)=4
      NTYPE(13)=4
      NTYPE(14)=3
      NTYPE(15)=1
      DO N=1,9
        NTYPE(25-N)=5-NTYPE(N+6)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE OF 3RD 25-52
C
      NTYPE(25)=4
      NTYPE(26)=2
      NTYPE(27)=1
      NTYPE(28)=3
      NTYPE(29)=4
      NTYPE(30)=3
      NTYPE(31)=1
      DO N=1,7
        NTYPE(39-N)=5-NTYPE(N+24)
        NTYPE(N+38)=NTYPE(N+24)
        NTYPE(53-N)=5-NTYPE(N+24)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE OF 4TH:53-100,6TH:165-260,9TH:517-708
C
      NTYPE(53)=4
      NTYPE(54)=2
      NTYPE(55)=1
      NTYPE(56)=4
      NTYPE(57)=3
      NTYPE(58)=1
      DO M=1,7
        DO N=1,6
          NTYPE(6*M+N+52)=NTYPE(N+52)
        END DO
      END DO
      DO M=1,16
        DO N=1,6
          NTYPE(6*M+N+158)=NTYPE(N+52)
        END DO
      END DO
      DO M=1,32
        DO N=1,6
          NTYPE(6*M+N+510)=NTYPE(N+52)
        END DO
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE OF 5TH:101-164, 7TH,8TH:261-516
C
      NTYPE(101)=4
      NTYPE(102)=2
      NTYPE(103)=3
      NTYPE(104)=1
      DO M=1,15
        DO N=1,4
          NTYPE(4*M+N+100)=NTYPE(N+100)
        END DO
      END DO
      DO M=1,64
        DO N=1,4
          NTYPE(4*M+N+256)=NTYPE(N+100)
        END DO
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TYPE CHECK
C
C      OPEN(8,FILE='TEST3.TXT')
      DO N=1,NEX
C        WRITE(8,'(6I5)') N,NTYPE(N)
      END DO
C      CLOSE(8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     AUTOMATIC NODE DEFINITION FOR 2ND 3RD ELEMENTS
C
      DO N=1,NEX
        NX(1,N,6)=NX(3,N,4)
        NX(3,N,6)=NX(1,N,4)
        NX(2,N,6)=NX(2,N,4)
        NX(1,N,5)=NX(1,N,4)
        NX(2,N,5)=NX(2,N,4)
        NX(3,N,5)=NX(3,N,4)
      IF(NTYPE(N).EQ.1) THEN
        NX(4,N,4)=NX(1,N,4)
        NX(4,N,6)=NX(2,N,6)
        NX(4,N,5)=NX(3,N,5)
        NZ(1,N,5)=3
        NZ(2,N,5)=2
        NZ(3,N,5)=2
        NZ(4,N,5)=3
      ELSE IF(NTYPE(N).EQ.4) THEN
        NX(4,N,4)=NX(3,N,4)
        NX(4,N,6)=NX(2,N,6)
        NX(4,N,5)=NX(1,N,5)
        NZ(1,N,5)=2
        NZ(2,N,5)=2
        NZ(3,N,5)=3
        NZ(4,N,5)=3
      ELSE IF(NTYPE(N).EQ.2) THEN
        NX(4,N,4)=NX(2,N,4)
        NX(4,N,6)=NX(1,N,6)
        NX(4,N,5)=NX(1,N,5)
        NZ(1,N,5)=2
        NZ(2,N,5)=3
        NZ(3,N,5)=2
        NZ(4,N,5)=3
      ELSE IF(NTYPE(N).EQ.3) THEN
        NX(4,N,4)=NX(2,N,4)
        NX(4,N,6)=NX(3,N,6)
        NX(4,N,5)=NX(3,N,5)
        NZ(1,N,5)=2
        NZ(2,N,5)=3
        NZ(3,N,5)=2
        NZ(4,N,5)=3
      ELSE
        WRITE(6,*) 'NTYPE MISSING',NTYPE(N)
      ENDIF
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    ELEMENT VOLUME CHECK
C
C      OPEN(8,FILE='TEST5.TXT')
      DO N=1,NEX
      DO L=4,6
        DO K=1,4
          X(K)=XX(NX(K,N,L))
          Y(K)=YY(NX(K,N,L))
          Z(K)=ZZ(NZ(K,N,L))
        END DO
          D1=X(2)*Y(3)*Z(4)+Y(2)*Z(3)*X(4)+Z(2)*X(3)*Y(4)
     &      -Z(2)*Y(3)*X(4)-X(2)*Z(3)*Y(4)-Y(2)*X(3)*Z(4)
          D2=X(3)*Y(4)*Z(1)+Y(3)*Z(4)*X(1)+Z(3)*X(4)*Y(1)
     &      -Z(3)*Y(4)*X(1)-X(3)*Z(4)*Y(1)-Y(3)*X(4)*Z(1)
          D3=X(4)*Y(1)*Z(2)+Y(4)*Z(1)*X(2)+Z(4)*X(1)*Y(2)
     &      -Z(4)*Y(1)*X(2)-X(4)*Z(1)*Y(2)-Y(4)*X(1)*Z(2)
          D4=X(1)*Y(2)*Z(3)+Y(1)*Z(2)*X(3)+Z(1)*X(2)*Y(3)
     &      -Z(1)*Y(2)*X(3)-X(1)*Z(2)*Y(3)-Y(1)*X(2)*Z(3)
          VOL(L-3)=D1-D2+D3-D4
      END DO
C        WRITE(8,'(I5,3F15.5)') N,(VOL(L),L=1,3)
      END DO
C      CLOSE(8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     COMMON NODES OF ADJACIANT ELEMETS
C
C      OPEN(8,FILE='TEST6.TXT')
      DO N=1,NEX
        NCOM12=0
        NCOM23=0
        NCOM13=0
        DO J=1,4
        DO K=1,4
          IF((NX(J,N,4).EQ.NX(K,N,5)).AND.(NZ(J,N,5).EQ.NZ(K,N,5))) THEN
            NCOM12=NCOM12+1
          ENDIF
          IF((NX(J,N,5).EQ.NX(K,N,6)).AND.(NZ(J,N,6).EQ.NZ(K,N,6))) THEN
            NCOM23=NCOM23+1
          ENDIF
          IF((NX(J,N,4).EQ.NX(K,N,6)).AND.(NZ(J,N,4).EQ.NZ(K,N,6))) THEN
            NCOM13=NCOM13+1
          ENDIF
        END DO
        END DO
C        WRITE(8,'(4I5)') N,NCOM12,NCOM23,NCOM13
      END DO
C      CLOSE(8)
C      OPEN(8,FILE='TEST4.TXT')
      DO N=1,NEX-1
      DO M=N+1,NEX
         DO I=4,6
         DO J=4,6
           ICCHI=0
           DO II=1,4
           DO JJ=1,4
           IF((NX(II,N,I).EQ.NX(JJ,M,J)).AND.(NZ(II,N,I).EQ.NZ(JJ,M,J)))
     &       ICCHI=ICCHI+1
           END DO
           END DO
C           IF(ICCHI.EQ.3) WRITE(8,'(4I5)') N,M,I,J
         END DO
         END DO
      END DO
      END DO
C      CLOSE(8)
C
c
C      open(8,file='nxc.txt')
      do n=1,nex
C        write(8,'(4i5)') (nx(k,n,4),k=1,4)
      end do
C      close(8)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     INVERSION TO J=1,3
C     modified
      DO N=1,NEX
        DO J=1,3
          NZ(1,N,J)=4-NZ(3,N,7-J)
          NZ(2,N,J)=4-NZ(2,N,7-J)
          NZ(3,N,J)=4-NZ(1,N,7-J)
          NZ(4,N,J)=4-NZ(4,N,7-J)
          NX(1,N,J)=NX(3,N,7-J)
          NX(2,N,J)=NX(2,N,7-J)
          NX(3,N,J)=NX(1,N,7-J)
          NX(4,N,J)=NX(4,N,7-J)
        END DO
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(NNODEZ.GT.3) THEN
C
C     EXTENDING TO THE UPPER ELEMENTS BY L
C     modified
      LL=(NNODEZ-1)/2-1
      DO L=1,LL
         DO J=1,6
            DO N=1,NEX
               DO K=1,4
                  NX(K,N,6*L+J)=NX(K,N,J)
                  NZ(K,N,6*L+J)=NZ(K,N,J)+2*L
               END DO
            END DO
         END DO
      END DO
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     OPTIMIZATION OF NODE NUMBERS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     CONNECTION MATRIX
C     modified
      DO N=1,NNODEX
         DO M=1,NNODEX
            NCOR(N,M)=0
         END DO
         NREN(N)=N
      END DO
      DO N=1,NEX
         DO J=1,3
            DO K=1,3
C               if((nz(j,n,4).ne.2).or.(nz(k,n,4).ne.2)) write(6,*) n,j,k
               NCOR(NX(J,N,4),NX(K,N,4))=1
            END DO
         END DO
      END DO
C      open(8,file='ncorc.txt')
      do n=1,nnodex
C         write(8,'(1000i2)') (ncor(m,n),m=1,nnodex)
      end do
C      close(8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      NINI=38
C      OPEN(10,FILE='TESTWD.TXT')
C
C     SORTING REPEAT FOR DIFFERENT NINI
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     modified
      DO NINI=38,38
C
         NHWD=0
         NSUMX=0
C     OPEN(88,FILE='NCN.TXT')
         DO N=1,NNODEX
            NSUM=0
            LL=1
            DO J=1,NNODEX
               IF((N.NE.J).AND.(NCOR(N,J).EQ.1)) THEN
                  NSUM=NSUM+1
                  LL=LL+1
                  NCN(LL,N)=J
                  NCN0(LL,N)=J
               ENDIF
            END DO
            NCN(1,N)=NSUM
            NCN0(1,N)=NSUM
C     WRITE(88,'(20I5)') N,(NCN(J,N),J=1,10)
            IF(NSUM.GT.NSUMX) NSUMX=NSUM
         END DO
C     CLOSE(88)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     NCN(1,N): ?m?[?hN?ɉ??̃m?[?h???q???��Ă??邩?H
C     NCN(K,N): K>1 ?m?[?hN?Ɍq???��Ă????m?[?h?̔ԍ?
C     NCN0 ?ŏ??̏??Ԃŕۑ??????邪?A
C     NCN ?͈??̃m?[?h???I?����????ƁA?????Ɍq???��Ă????m?[?hN?ł?
C         NCN(1,N) ?????????A?I?����ꂽ?m?[?h?ԍ????��????????B
C         ?��????Ƃ? SUBROUTINE DELNREN ?ōs???????B
C
CCCCCCCCCCCCCCCCCCCCCC
C
C     SORTING INITIALIZE
C
         NLAST=0
         DO J=1,NINI
            JR=NNAIEND+J
C
            NREN(J)=JR
            CALL DELNREN(JR)
            NLAST=JR
         END DO
C         OPEN(88,FILE='NCN2.TXT')
         DO N=1,NNODEX
C           WRITE(88,'(20I5)') N,(NCN(J,N),J=1,10)
         END DO
C         CLOSE(88)
         NR2=NCN(1,NLAST)
         LLL=NINI
C
CCCCCCCCCCCCCCCCCCCC
C
C     SORTING START
C
C     NREN(JJ)?FJJ?ԖڂɑI?����ꂽ?m?[?h?ԍ??ł????B???̏??ɔԍ????t???ւ????????B
C         ?ł���???I?����ꂽ?m?[?h?Ɍq???��Ă????m?[?h?����̃m?[?h?Ƃ??đI?����??B
C         ???̎??A???E?????̃m?[?h???????΂??????D?悷???B
C     NHWD?F?m?[?h???I?����????́A?K?????ԎႢ?m?[?h?Ɍq???��Ă????󂾂????A
C           ?I?����??????ԍ??Ƃ??̂Ƃ??̈??ԎႢ?m?[?h?̏????ԍ??Ƃ̍??????��Ă?????
C           ?q???��Ă????m?[?h?Ԃ̏????ԍ??ԋ????̍ő??l???????????B
C
C         OPEN(7,FILE='TEST8.TXT')
         DO JJ=1,NNODEX-1
C     WRITE(6,*) JJ
            JR=NREN(JJ)
C     WRITE(6,*) JR,NCN(1,JR)
            IF(NCN(1,JR).NE.0) THEN
               LIMIT=10
               NLIM=0
 1000          CONTINUE
               NLIM=NLIM+1
               NR1=NCN(1,JR)
               NEXT=NCN(2,JR)
               DO J=2,NR1+1
C
                  IF(NCN(J,JR).GT.(NNAIEND+1)) THEN
C
                     NEXT=NCN(J,JR)
                  ENDIF
               END DO
C               WRITE(7,'(12I5)') NEXT,(NCN0(K,NEXT),K=1,10)
               CALL DELNREN(NEXT)
C               WRITE(7,'(12I5)') JR,(NCN(K,JR),K=1,10)
               NLAST=NEXT
               LLL=LLL+1
               NHF=LLL-JJ+1
               IF(NHWD.LT.NHF) NHWD=NHF
               NREN(LLL)=NEXT
C     WRITE(6,*) JR,NCN(1,JR),NEXT,JJ,LLL
C               WRITE(7,'(10I5)') JR,NCN(1,JR),NEXT,JJ,LLL,NHF
               IF(NLIM.GE.LIMIT) THEN
                  WRITE(6,*) 'NLIM,LIMIT',NLIM,LIMIT
                  STOP
               ENDIF
               IF(NCN(1,JR).NE.0) GOTO 1000
            ENDIF
         END DO
C         CLOSE(7)
C         OPEN(8,FILE='TEST7.TXT')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     CORRELATION WIDTH
C
C     KWD?F?m?[?h???I?����ꂽ???ɕ��ׂĂ??????̊Ԃ̐ڑ????Ԃ?
C          ?}?g???b?N?X?ŕ\???????ɁA???̃m?[?h???q???��Ă???
C          ?ł���???m?[?h???ԍ??̏????????Ƒ傫?????Əo???邪
C          ???̍??i?g?????j?ł????B2*NHWD ?????͏??????B
C
C          OPEN(12,FILE='TEST9.TXT')
          KWD=0
          DO J=1,NNODEX
             JN=NREN(J)
             DO K=1,NNODEX
                KN=NREN(K)
                IF(NCOR(JN,KN).EQ.1) KMAX=K
                KN=NREN(NNODEX+1-K)
                IF(NCOR(JN,KN).EQ.1) KMIN=NNODEX+1-K
             END DO
C             WRITE(12,'(10I5)') J,JN,KMIN,KMAX,NREN(KMIN),NREN(KMAX)
             IF(KWD.LT.(KMAX-KMIN+1)) KWD=KMAX-KMIN+1
          END DO
C          WRITE(6,*) NINI, ' KWD= ',KWD,' NHWD=',NHWD
C          CLOSE(12)
C          WRITE(10,*) NINI,KWD,NHWD
          DO J=1,NNODEX
C             WRITE(8,'(2I5,2F15.5)') J,NREN(J),XX(NREN(J)),YY(NREN(J))
          END DO
C          CLOSE(8)
       END DO
C       CLOSE(10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       OPEN(10,FILE='TEST10.TXT')
       DO J=1,NNODEX
          DO K=1,NNODEX
         IF(NREN(K).EQ.J) NRENREV(J)=K
      END DO
C      WRITE(10,'(2I5)') J,NRENREV(J)
      END DO
C      CLOSE(10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     OUTPUT OF FEMDATA
C
C      OPEN(15,FILE='H9M5DATA.QQQ')
C      WRITE(15,'(2I6)') NEX,NEZ
      DO J=1,NEX
         DO K=1,NEZ
            DO L=1,4
              NREV(L)=   NRENREV(NX(L,J,K))
              NX2(L,J,K)=NRENREV(NX(L,J,K))
            END DO
C            WRITE(15,'(10I5)') J,K,(NREV(L),NZ(L,J,K),L=1,4)
         END DO
      END DO
C      WRITE(15,'(2I6)') NNODEX,NNODEZ
      DO J=1,NNODEX
C         IF (XX(NREN(J)).LT.0.1D-16) XX(NREN(J)) = 0D0
C         IF (YY(NREN(J)).LT.0.1D-16) YY(NREN(J)) = 0D0
C         WRITE(15,'(I5,2E24.9)') J,XX(NREN(J)),YY(NREN(J))
C         WRITE(15,'(I5,2E24.15)') J,XX(NREN(J)),YY(NREN(J))
        XXL(J)=XX(NREN(J))
        YYL(J)=YY(NREN(J))
      END DO
      DO J=1,NNODEZ
C         WRITE(15,'(I5,2E20.9)') J,ZZ(J)
      END DO
C      WRITE(15,'(I5)') MXB
      DO J=1,MXB
C         WRITE(15,'(I5)') NRENREV(J+NNODEX-MXB)
         IBND(J) = NRENREV(J+NNODEX-MXB)
      END DO
C      WRITE(15,'(12E15.5)') R1,R2,R3,R4
C      CLOSE(15)
C      STOP
      RETURN
      END SUBROUTINE B9C3D128
C
CCCCCCCCCCCCCCCCCCCCCCCC
C
C     ?I?����ꂽ???_?ł????Ɍq?????ߓ_?????̑I?????r??????
C
      SUBROUTINE DELNREN(NEXT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MXN=419)
      COMMON NCN(10,MXN),NCN0(10,MXN)
      NR0=NCN0(1,NEXT)
      DO J=2,NR0+1
         JNEXT=NCN0(J,NEXT)
         NR2=NCN(1,JNEXT)
         DO K=2,NR2+1
            IF(NCN(K,JNEXT).EQ.NEXT) THEN
               NCN(K,JNEXT)=NCN(NR2+1,JNEXT)
               NCN(1,JNEXT)=NCN(1,JNEXT)-1
            ENDIF
         END DO
      END DO
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

