C=======================================================================
C
C        STAMPANJE NAPONA TROUGAONE PLOCE 
C
C   SUBROUTINE PL5NAP
C              STANP5
C           *  STANG5
C              NAPL3M
C
C                NAPL3M (DEB,NAPV,INDNA,NDEB)
C                MATKEP (SKE,DEB,DEBNLM,NLM,MATP,NG,KOR)
C=======================================================================
      SUBROUTINE PL5NAP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C    GLAVNI PROGRAM ZA STAMPANJE NAPONA PLOCA
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      LAU=LMAX
C     UCITAVANJE PODATAKA O ELEMENTIMA I NAPONA SA DISKA
      CALL READP(A(LAU))
C     PODPROGRAM ZA STAMPANJE NAPONA
      IF(INDPR.EQ.0.OR.INDPR.EQ.2) CALL STANP5(A(LAU))
C     PODPROGRAM ZA STAMPANJE NAPONA ZA POSTPROCESIRANJE
C     IF(INDGR.EQ.0.OR.INDGR.EQ.2) CALL STANG5(A(LAU))
C     IF(IATYP.EQ.0) RETURN
C     IF(NMODM.LE.4) RETURN
C     NPROS=NE*NGAUSX*NGAUSY*16*IDVA
C     LMA8=LSTAZA(3)-1
C     CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMA8,LDUZI)
      RETURN
      END
C======================================================================
      SUBROUTINE STANP5(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      include 'paka.inc'
      
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON/PL3CEL/ NPP,NGAUS,INULAZ,KDT,
     1IPP,IPMK,IPP2,MAXSK,INDPF,MATP,INDNA,NAPONG(3,5),LM(18),NDEB
      COMMON /PL3REP/ LSKEB,LSKEM,LEMB,LSIGB,LEN,LFUE
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TEMPCV/ LTECV,ITEMP
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
C
      DIMENSION AU(*)
      REAL AU
C
      GO TO (1,1,3,3,5,5,5,5),NMODM
    1 CALL NAPL3M(AU(LTHID),A(LSIGB),A(LEN),AU(LISNA),INDNA,NDEB,
     1            AU(LNSLOJ),AU(LDSLOJ),AU(LNMAT))
      RETURN
C   3 CALL STAN31(A(LSIGMA),AU(LISNA),AU(LNEL),A(LTECV))
    3 RETURN
C   5 CALL STAN51(A(LPLAST),AU(LISNA),AU(LIPRC))
    5 RETURN
      END
C======================================================================
      SUBROUTINE NAPL3M (DEB,SIGB,EN,NAPV,INDNA,NDEB,NSLOJ,DSLOJ,IMAT)
C=======================================================================
C     ........ RAC. I STAMPA NAPONE U TEZ. ELEM.  .............
C=======================================================================
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION EN(NE,3,*),SIGB(NE,3,*)
      DIMENSION SGL(3),SGME(3), DEB(*), NAPV(*)
      DIMENSION NSLOJ(*),DSLOJ(MSLOJ,*),IMAT(*)
      PI=1.D0
      PI=4.D0*DATAN(PI)
      STEP=180.D0/PI
C
C....  PETLJA PO ELEMENTIMA
      INDSN=0
      DO 44 NLM=1,NE
      DEB0=DEB(NLM)
      DEBNLM=DEB0
      MST=IMAT(NLM)
      NNSL=1
      IF(MSET.GT.0)NNSL=NSLOJ(MST)
      DO 44 MSL=1,NNSL
      IF(NAPV(NLM).EQ.1) GO TO 44
      IF(INDSN.EQ.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      INDSN=1
      ENDIF
        IF(MSET.GT.0) THEN
           DDD =DSLOJ(MSL,MST)
           DEBNLM=DEB0*DDD
        ENDIF
      DO 45 I=1,3
      SGME(I)=EN(NLM,I,MSL)/DEBNLM
      IF(I.NE.3) GO TO 51
      SGL(3)=SGME(3)+SIGB(NLM,3,MSL)
      GO TO 45
  51  IF(SGME(I)) 46,46,47
  46  IF(SIGB(NLM,I,MSL).LE.0.) GO TO 48
  49  SGL(I)=SGME(I)-SIGB(NLM,I,MSL)
      GO TO 45
  47  IF(SIGB(NLM,I,MSL)) 49,49,48
  48  SGL(I)=SGME(I)+SIGB(NLM,I,MSL)
  45  CONTINUE
      IF(SGL(1)-SGL(2))150,151,152
  150 IF(SGL(3))153,154,154
  153 ALFP=.5*DATAN(2.*SGL(3)/(SGL(1)-SGL(2)))
      ALFA=90.-ALFP*STEP
      GO TO 160
  154 ALFP=.5*DATAN(2.*SGL(3)/(SGL(1)-SGL(2)))
      ALFA=-90.+ALFP*STEP
      GO TO 160
  152 IF(SGL(3))155,156,156
  155 ALFA=-.5*DATAN(2.*SGL(3)/(SGL(2)-SGL(1)))*STEP
      GO TO 160
  156 ALFA=.5*DATAN(2.*SGL(3)/(SGL(1)-SGL(2)))*STEP
  160 TMAX=DSQRT(.25*(SGL(1)-SGL(2))**2+SGL(3)*SGL(3))
      SG1=.5*(SGL(1)+SGL(2))+TMAX
      SG2=SG1-2.*TMAX
      GO TO 63
C     HIDROSTATICKO STANJE NAPONA
  151 SG1=SGL(1)
      SG2=SGL(2)
      ALFA=0.D0
      TMAX=0.D0
  63  IF(MSL.EQ.1)WRITE(IZLAZ,1000) NLM
      WRITE(IZLAZ,1010)(SIGB(NLM,I,MSL),I=1,3),(SGME(I),I=1,3),SG1,SG2
     1,TMAX,ALFA
  44  CONTINUE
      RETURN
C
 1000 FORMAT(I5)
 1010 FORMAT(2X,2(2X,3(1X,1PD11.4))/
     1       20X,2X,2(1X,1PD11.4),4X,2(2X,1PD11.4))
C-----------------------------------------------------------------------
 2000 FORMAT(///21X,'N A P O N I  U TEZISTU ELEMENTA PLOCA'//
     1' ','ELEMENT',9X,'  S A V O J N I  ',18X,'M E M B R A N S K I',/
     3 11X,'SX',10X,'SY',9X,'SXY',12X,'SX',10X,'SY',9X,'SXY'/
     28X,'   GL. NAPONI',8X,'S1',10X,'S2',10X,'TAUMAX',7X,'ALFA')
C-----------------------------------------------------------------------
 6000 FORMAT(///21X,'S T R E S S   I N   P L A T E   E L E M E N T S'//
     1' ','ELEMENT',9X,'  B E N D I N G  ',18X,'  M E M B R A N E  ',/
     3 11X,'SX',10X,'SY',9X,'SXY',12X,'SX',10X,'SY',9X,'SXY'/
     28X,'PRINC. STRESS',8X,'S1',10X,'S2',10X,'TAUMAX',7X,'ALFA')
C-----------------------------------------------------------------------
      END