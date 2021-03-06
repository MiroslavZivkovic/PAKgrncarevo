C=======================================================================
C
C        UCITAVANJE PODATAKA O OSNOSIMETRICNOJ LJUSCI
C
C   SUBROUTINE UL7EGL
C              UL7EK
C              UL7EK7
C              LM7MHT
C              UCGRP7
C              PSPOL7
C              TGRAF7
C
C=======================================================================
      SUBROUTINE UL7EGL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
C
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA ZA
C     OSNOSIMETRICNU LJUSKU
C
      include 'paka.inc'
      
      COMMON /ELEMA7/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1                LBETA
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /AXISYM/ INDAX
C
C     OSNOVNI PODACI O ELEMENTIMA
C
      INDAX=1
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) MSLOJ,NGAUSX,NGAUSY
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) MSLOJ,NGAUSX,NGAUSY
      IF(NGAUSX.EQ.0) NGAUSX=2
      IF(NGAUSY.EQ.0) NGAUSY=2
      IF(MSLOJ.EQ.0) MSLOJ=1
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1CALL WBROJK(KARTIC,0)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,2000) MSLOJ,NGAUSX,NGAUSY
C
C     REPERI U VEKTORU A ZA ULAZNE PODATKE
C
      NCVE=4
      LTHID=1
      LBETA=LTHID+NE*MSLOJ*NCVE*IDVA
      LNMAT=LBETA+NE*MSLOJ*IDVA
      LIPGC=LNMAT+NE*MSLOJ
      LISNA=LIPGC+NE
      LIPRC=LISNA+NE
      LNEL=LIPRC+NE
      LAU=LMAX
C
C     POZIVANJE PROGRAMA ZA ULAZNE PODATKE U VEKTOR AU
C
      CALL UL7EK(A(LAU))
      NCVE3=NCVE*3
      MXAE = NCVE3+(MSLOJ+2*NCVE+9*NCVE3+NCVE3*(NCVE3+1)/2)*IDVA
      LMAX = LMAX + MXAE
C
      RETURN
C
 1000 FORMAT(14I5)
C-----------------------------------------------------------------------
 2000 FORMAT(
     111X,'NUMBER OF LAYERS IN AXSISYMETRIC SHELLS ........ MSLOJ =',I5/
     116X,'EQ.0; MSLOJ = 1'///
     611X,'NUMBER OF GAUSS POINTS IN R-DIRECTION ......... NGAUSX =',I5/
     616X,'EQ.0; NGAUSX = 2'///
     611X,'NUMBER OF GAUSS POINTS IN S-DIRECTION ......... NGAUSY =',I5/
     616X,'EQ.0; NGAUSY = 2')
CS 2000 FORMAT(
CS     111X,'BROJ SLOJEVA OSNOSIMETRICNE LJUSKE ............. MSLOJ =',I5/
CS     116X,'EQ.0; MSLOJ = 1'///
CS     611X,'BROJ GAUSOVIH TACAKA U PRAVCU OSE KSI ......... NGAUSX =',I5/
CS     616X,'EQ.0; NGAUSX = 2'///
CS     611X,'BROJ GAUSOVIH TACAKA U PRAVCU OSE ETA ......... NGAUSY =',I5/
CS     616X,'EQ.0; NGAUSY = 2')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UL7EK(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
C
      include 'paka.inc'
      
      COMMON /ELEMA7/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1                LBETA
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
C
      DIMENSION AU(*)
      REAL AU
      NCVE3=NCVE*3
C
C     POZIVANJE PROGRAMA ZA ULAZNE PODATKE O OSNOSIM. LJUSCI
C
      CALL ICLEAR(AU(LNEL),NE*NCVE)
C
      CALL UL7EK7(AU(LNEL),AU(LNMAT),AU(LIPGC),AU(LISNA),AU(LIPRC),
     1AU(LTHID),AU(LBETA))
C ....    GRAFIKA
      IF(NBLGR.GE.0) CALL TGRAF7(AU(LNEL))
C
      LLMEL=LNEL+NE*NCVE
      LMXAU=LLMEL+NE*NCVE3
      MXAU = LMXAU - 1
      CALL DELJIV(MXAU,2,INDL)
      IF(INDL.EQ.1) MXAU=MXAU+1
      LMAX=LAU+MXAU
      IF (LMAX.LE.MTOT) GO TO 5
      WRITE(IZLAZ,2005) LMAX,MTOT
      STOP
C
C     FORMIRANJE VEKTORA LM I VISINA STUBOVA
C
    5 NCVE3=NCVE*3
      CALL ICLEAR(AU(LLMEL),NE*NCVE3)
      CALL LM7MHT(A(LID),AU(LNEL),AU(LLMEL),NCVE3)
C
      LMAX8=LMAX8+1
      WRITE(IELEM,REC=LMAX8)
     1 MSLOJ,NGAUSX,NGAUSY,NCVE,ITERME,
     1 MXAU,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,LLMEL,LBETA
      CALL WRITDD(AU(LTHID),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      IF(NMODM.LE.4) GO TO 10
      STOP '***  PALSTICNOST NIJE UGRADJEN ZA TIP ELEMENTA  7  ***'
C     LPLAST=LMAX
C     LPLAS1=LMAX
C     LSIGMA=LMAX
C     NPROS=NE*NGAUSX*NGAUSY*16*IDVA
C     LPLAS1=LPLAST+NPROS
C     LMAX=LPLAS1+NPROS
C     IF(LMAX.GT.MTOT) CALL ERROR(1)
C     CALL CLEAR(A(LPLAST),NPROS*2/IDVA)
C     CALL WRITDD(A(LPLAST),NPROS/IDVA,IELEM,LMAX8,LDUZI)
C     CALL WRITDD(A(LPLAS1),NPROS/IDVA,IELEM,LMAX8,LDUZI)
C     RETURN
   10 LSIGMA=LMAX
C ...   6 KOMPONENTI NAPONA
      NPROS=NE*MSLOJ*NGAUSX*NGAUSY*6*IDVA
      LMAX=LSIGMA+NPROS
      IF(LMAX.GT.MTOT) CALL ERROR(1)
      CALL CLEAR(A(LSIGMA),NPROS/IDVA)
      CALL WRITDD(A(LSIGMA),NPROS/IDVA,IELEM,LMAX8,LDUZI)
      RETURN
C
C-----------------------------------------------------------------------
 2005 FORMAT(///' NOT ENOUGH SPACE IN WORKING VECTOR  A'
     1/' REQUESTED DIMENSION , LMAX=',I10
     2/' AVAILABLE DIMENSION , MTOT=',I10)
CS 2005 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA PODATKE O OSNOSIME
CS     1TRICNOJ LJUSCI'/' POTREBNA DIMENZIJA , LMAX=',I10/' RASPOLOZIVA DI
CS     2MENZIJA NTOT =',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UL7EK7(NEL,MATV,IPGCV,ISNAP,IPRCV,THICK,BETA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
C
C     PODPROGRAM ZA UCITAVANJE PODATAKA O IZOPARM. KRIVOLIN. ELEM. K09
C
      include 'paka.inc'
      
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      DIMENSION NEL(NE,*),MATV(NE,*),IPGCV(*),ISNAP(*),IPRCV(*),
     1          THICK(NE,MSLOJ,*),BETA(NE,*),THICS(4)
C
C     PODACI O ELEMENTIMA
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,2000)
      NMATS=  1
      IPRCOS=-1
      IPGSS= -1
      ISNAA= -1
      DO 4 I=1,4
    4 THICS(I)=  1.
      I = 0
      NAUT=0
    5 I = I + 1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NN,IPRCO,ISNA,IPGS,KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NN,IPRCO,ISNA,IPGS,KORC
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (NEL(NN,J),J=1,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) (NEL(NN,J),J=1,4)
C
      DO 10 MS=1,MSLOJ
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (THICK(NN,MS,J),J=1,4),NMAT,BETA(NN,MS)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) (THICK(NN,MS,J),J=1,4),NMAT,BETA(NN,MS)
        DO 7 J=1,4
        IF (NEL(NN,J).EQ.0) GO TO 7
        IF (DABS(THICK(NN,MS,J)).LT.1.D-20) THICK(NN,MS,J)=THICS(J)
        THICS(J)=THICK(NN,MS,J)
   7  CONTINUE
      IF(NMAT.EQ.0) NMAT=NMATS
      NMATS=NMAT
      MATV(NN,MS)= NMAT
10    CONTINUE
C
      IF (IPRCO.EQ.0) IPRCO=IPRCOS
      IPRCOS=IPRCO
      IF(IPRCO.LT.0) IPRCO = 0
      IPRCV(NN)=IPRCO
      IF(ISNA.EQ.0) ISNA=ISNAA
      ISNAA=ISNA
      IF(ISNA.LT.0) ISNA = 0
      ISNAP(NN)=ISNA
      IF(IPGS.EQ.0) IPGS=IPGSS
      IPGSS=IPGS
      IF(IPGS.LT.0) IPGS = 0
      IPGCV(NN)=IPGS
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3) THEN
      WRITE(IZLAZ,2001) NN,(NEL(NN,J),J=1,4),IPRCO,ISNA,IPGS,KORC
      DO 15 MS=1,MSLOJ
15    WRITE(IZLAZ,2002)MS,(THICK(NN,MS,J),J=1,4),MATV(NN,MS),BETA(NN,MS)
                                   ENDIF
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(I.EQ.NE) GO TO 50
      GO TO 5
C
   20 NAUT=1
      NN1=NN
      KORA=KORC
      GO TO 5
C
C  AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
C
   30 NN2=NN
      N11=NEL(NN1,1)
      N22=NEL(NN2,1)
      N12=N22-N11
      CALL DELJIV(N12,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      DO 40 NC=2,2
      N11=NEL(NN1,NC)
      N22=NEL(NN2,NC)
      N21=N22-N11
      IF(N12.NE.N21) GO TO 100
   40 CONTINUE
      DO 45 NC=3,4
      N11=NEL(NN1,NC)
      N22=NEL(NN2,NC)
      IF(N11.EQ.0.AND.N22.EQ.0)GO TO 45
      N21=N22-N11
      IF(N12.NE.N21)GO TO 100
   45 CONTINUE
      NNN=NN2-NN1-1
      NGG=N12/KORA-1
      IF(NNN.NE.NGG) GO TO 150
      IAUT=N12/KORA-1
      DO 34 J=1,IAUT
      JJ=NN1+J
      IPRCV(JJ)=IPRCV(NN1)
      ISNAP(JJ)=ISNAP(NN1)
      IPGCV(JJ)=IPGCV(NN1)
        DO 33 MS=1,MSLOJ
          DO 32 MC=1,4
   32     THICK(JJ,MS,MC)=THICK(NN1,MS,MC)
          BETA(JJ,MS) = BETA(NN1,MS)
   33     MATV(JJ,MS) = MATV(NN1,MS)
      DO 35 NN=1,4
      NODP=NEL(JJ-1,NN)
      NEL(JJ,NN)=NODP+KORA
      IF(NODP.EQ.0) NEL(JJ,NN) = 0
   35 CONTINUE
   34 CONTINUE
      I=I+IAUT
      IF(I.EQ.NE) GO TO 50
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      NN1=NN2
      GO TO 5
C
C     ODREDJIVANJE MAKSIMALNOG BROJA
C
   50 NCVE=2
      DO 60 I=1,NE
      DO 60 J=3,4
      IJ=NEL(I,J)
      IF(IJ.GT.0.AND.J.GT.NCVE) NCVE=J
   60 CONTINUE
C
C  STAMPANJE UCITANIH I GENERISANIH PODATAKA O CVOROVIMA ELEMENATA
C
      IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      WRITE(IZLAZ,2140)
      DO 70 I=1,NE
      WRITE(IZLAZ,2001) I,(NEL(I,J),J=1,4),
     1IPRCV(I),ISNAP(I),IPGCV(I)
      DO 65 MS=1,MSLOJ
   65 WRITE(IZLAZ,2002) MS,(THICK(I,MS,J),J=1,4),NMAT,BETA(I,MS)
   70 CONTINUE
C.... PRETVORI BETA U RADIJANE
      DO 90 NN=1,NE
      DO 90 MS=1,MSLOJ
   90 BETA(NN,MS)=BETA(NN,MS)*3.1415927D0/180.
      RETURN
C
  150 WRITE(IZLAZ,2110) NN1,NN2,NGG
      STOP
  100 WRITE(IZLAZ,2100) N22,N11,KORA
      STOP
C
 1000 FORMAT(14I5)
 1010 FORMAT(4F10.0,I5,F10.0)
 2001 FORMAT(11X,I5,8X,4I6,11X,I6,3X,I6,4X,I6,5X,I6)
 2002 FORMAT(13X,'(',I2,')',7X,4F6.3,2X,I6,33X,F6.3)
C-----------------------------------------------------------------------
 2000 FORMAT(///6X,'R E A D I N G    D A T A    A B O U T    A X I S Y M
     1 E T R I C    S H E L L'/6X,75('-')///11X,
     3'ELEMENT',7X,'NODE1 NODE2 NODE3 NODE4',5X,'    ',4X,'COORD.',4X,
     4'STRESS',4X,'STRESS',3X,'STEP'  /11X,'/LAYER',9X,
     4'/NODE THICKNESS       ',3X,'/MATER.',3X,' PRINT',4X,' PRINT',4X,
     4'AT NODE',2X,'/BETA')
CS 2000 FORMAT(///6X,'U C I T A N I    P O D A C I    O    O S N O S I M E
CS     1 T R I C N O J    LJ U S C I'/6X,80('-')///
CS     111X,'(MOGUCE GENERISANJE PODATAKA O ELEMENTIMA,'/
CS     212X,'UCITATI POTREBAN BROJ KARTICA ZA DEF. SVIH ELEMENATA)'///11X,
CS     3'ELEMENT',7X,'CVOR1 CVOR2 CVOR3 CVOR4',5X,'    ',4X,'KOORD.',4X,
CS     4'NAPON',5X,'NAPON',4X,'KORAK'  /11X,'/SLOJ',10X,
CS     4'/DEBLJINE U CVOROVIMA ',3X,'/MATER.',3X,'STAMPA',4X,'STAMPA',4X,
CS     4'CVORA',4X,'/BETA')
 2140 FORMAT(6X,'G E N E R A T E D    D A T A    A B O U T    A X S I S
     1Y M E T R I C    S H E L L'/6X,81('-')///11X,
     3'ELEMENT',7X,'NODE1 NODE2 NODE3 NODE4',5X,'    ',4X,'COORD.',4X,
     4'STRESS',4X,'STRESS',           /11X,'/LAYER',9X,
     4'/NODE THICKNESS       ',3X,'/MATER.',3X,' PRINT',4X,' PRINT',4X,
     4'AT NODE',2X,'/BETA')
CS 2140 FORMAT(6X,'G E N E R I S A N I    P O D A C I    O    O S N O S I
CS     1M E T R I C N O J    LJ U S C I'/6X,86('-')///11X,
CS     3'ELEMENT',7X,'CVOR1 CVOR2 CVOR3 CVOR4',5X,'    ',4X,'KOORD.',4X,
CS     4'NAPON',5X,'NAPON'                         /11X,'/SLOJ',10X,
CS     4'/DEBLJINE U CVOROVIMA ',3X,'/MATER.',3X,'STAMPA',4X,'STAMPA',4X,
CS     4'CVORA',4X,'/BETA')
 2110 FORMAT(///' AMONG ELEMENTS N1=',I5,' AND N2=',I5,' THERE IS NO  ',
     1'NG=',I5,' ELEMENTS TO GENERATE')
CS 2110 FORMAT(///' IZMEDJU ELEMENATA N1=',I5,' I N2=',I5,' NEMA NG=',I5,
CS     1' ELEMENATA KOJE TRBA GENERISATI')
 2100 FORMAT(///' NODE NUMBER N2=',I5,' CAN NOT BE DEFINED BY SUPERPOSIT 
     1ION OF NODE N1=',I5,' WITH A NUMBER OF STEPS  KORA=',I5)
CS 2100 FORMAT(///' BROJ CVORA N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJ
CS     1A CVORA N1=',I5,' I KONACNOG BROJA KORAKA KORA=',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE LM7MHT(ID,NEL,LMEL,NCVE3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORMIRANJE VEKTORA LM I VISINA STUBOVA
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
C
      DIMENSION ID(NP,*),NEL(NE,*),LMEL(NCVE3,*)
C     PETLJA PO ELEMENTIMA
      DO 100 NLM=1,NE
      DO 10 NC=1,NCVE
      IF(NEL(NLM,NC).EQ.0) GO TO 10
      NNC=3*(NC-1)
      JJ=NEL(NLM,NC)
      DO 20 I=1,3
      II=I
      IF(I.EQ.3) II=6
   20 LMEL(NNC+I,NLM)=ID(JJ,II)
   10 CONTINUE
C
C     FORMIRANJE VISINA STUBOVA
C
      CALL VISINE(A(LMHT),NCVE3,LMEL(1,NLM))
      IF (IABS(ICCGG).EQ.1) THEN
         ND=NCVE3
         WRITE(IDRAKCE) ND,(LMEL(I,NLM),I=1,ND)
         NELUK=NELUK+1
      ENDIF
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE UCGRP7
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     REPERI  ZA PRITISKE NA GRANICI OSNOSIMETRICNE LJUSKE
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /OPTERE/ NCF,NPP2,NPP3,NPGR,NPGRI,NPLJ,NTEMP
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
C
      LFAKP = LRAD
      LNFUN = LFAKP + NPLJ*2*IDVA
      LIPRAV = LNFUN + NPLJ
      LNODPR = LIPRAV + NPLJ
      LMAX = LNODPR + NPLJ*4
      IF(ISOPS.GT.0) LMAX=LMAX+JEDN*IDVA
      IF(NDIN.EQ.0) GO TO 10
      IF(MDVI.EQ.1.AND.ISOPS.EQ.0) LMAX=LMAX+JEDN*IDVA
   10 IF(LMAX.GT.MTOT) CALL ERROR(1)
C
C     UCITAVANJE POVRSINSKIH PRITISAKA NA GRANICI OSNOSIM. LJUSKE
C
      CALL PSPOL7(A(LNFUN),A(LIPRAV),A(LFAKP),A(LNODPR),NPLJ)
      NPROS=LMAX-LRAD
      IF(NRAD.LT.NPROS) NRAD=NPROS
      CALL UCELP(A(LIPODS),10)
      CALL WRITDD(A(LFAKP),NPLJ*2,IPODS,LMAX13,LDUZI)
      CALL IWRITD(A(LNFUN),NPLJ*6,IPODS,LMAX13,LDUZI)
      IF(MAX13.LT.LMAX13) MAX13=LMAX13
      RETURN
      END
C======================================================================
      SUBROUTINE PSPOL7(NFUN,IPRAV,FAKP,NODPR,NPP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
C
C     PODPROGRAM ZA UCITAVANJE ULAZNIH PODATAKA O POVRSINSKIM
C     PRITISCIMA NA GRANICI 2/D ELEMENTA
C
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
C
      DIMENSION NFUN(*),IPRAV(*),FAKP(NPP,*),NODPR(NPP,*)
C
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,2000)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,2010)
      I = 0
      NAUT=0
    5 I = I + 1
      CALL ISPITA(ACOZ)
      IF(I.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NFUN(I),IPRAV(I),(FAKP(I,J),J=1,2),KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NFUN(I),IPRAV(I),(FAKP(I,J),J=1,2),KORC
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (NODPR(I,J),J=1,4)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1030) (NODPR(I,J),J=1,4)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,2020) (NODPR(I,J),J=1,4),
     1NFUN(I),IPRAV(I),(FAKP(I,J),J=1,2),KORC
      IF(NAUT.GT.0) GO TO 30
      IF(KORC.NE.0) GO TO 20
      IF(I.EQ.NPP) GO TO 50
      GO TO 5
C
   20 NAUT=1
      N1=NODPR(I,1)
      KORA=KORC
      GO TO 5
C
C  AUTOMATSKO GENERISANJE PODATAKA IZMEDJU CVOROVA N1 I N2
C
   30 N2=NODPR(I,1)
      N11=NODPR(I-1,1)
      N22=NODPR(I,1)
      N12=N22-N11
      CALL DELJIV(N12,KORA,INDD)
      IF(INDD.EQ.1) GO TO 100
      N11=NODPR(I-1,2)
      N22=NODPR(I,2)
      N21=N22-N11
      IF(N12.NE.N21) GO TO 100
      N11=NODPR(I-1,3)
      N22=NODPR(I,3)
      IF(N11.EQ.0.AND.N22.EQ.0) GO TO 31
      N21=N22-N11
      IF(N12.NE.N21) GO TO 100
   31 N11=NODPR(I-1,4)
      N22=NODPR(I,4)
      IF(N11.EQ.0.AND.N22.EQ.0) GO TO 32
      N21=N22-N11
      IF(N12.NE.N21) GO TO 100
   32 RKORA = KORA
      RN1N2 = N2-N1
      DD = RKORA/RN1N2
      DF1F1 = DD*(FAKP(I,1) - FAKP(I-1,1))
      DF2F2 = DD*(FAKP(I,2) - FAKP(I-1,2))
      IAUT=(N2-N1)/KORA-1
      II=I-1
      I=I+IAUT
      NFUN(I) = NFUN(II+1)
      IPRAV(I) = IPRAV(II+1)
      FAKP(I,1) = FAKP(II+1,1)
      FAKP(I,2) = FAKP(II+1,2)
      DO 33 NN=1,4
   33 NODPR(I,NN)=NODPR(II+1,NN)
      DO 34 J=1,IAUT
      JJ=II+J
      NFUN(JJ) = NFUN(II)
      IPRAV(JJ) = IPRAV(II)
      FAKP(JJ,1) = FAKP(JJ-1,1)+DF1F1
      FAKP(JJ,2) = FAKP(JJ-1,2)+DF2F2
      DO 35 NN=1,4
      NODP=NODPR(JJ-1,NN)
      NODPR(JJ,NN)=NODP+KORA
      IF(NODP.EQ.0) NODPR(JJ,NN) = 0
   35 CONTINUE
   34 CONTINUE
      IF(I.EQ.NPP) GO TO 50
      NAUT=0
      IF(KORC.EQ.0) GO TO 5
      KORA=KORC
      NAUT=1
      N1=N2
      GO TO 5
C
C  STAMPANJE UCITANIH I GENERISANIH PODATAKA O CVOROVIMA NA GRANICI
C
   50 IF(NULAZ.NE.1.AND.NULAZ.NE.3) RETURN
      CALL WBROJK(KARTI,0)
      WRITE(IZLAZ,2040)
      DO 70 I=1,NPP
      WRITE(IZLAZ,2020) (NODPR(I,J),J=1,4),
     1NFUN(I),IPRAV(I),(FAKP(I,J),J=1,2)
   70 CONTINUE
      RETURN
C
  100 WRITE(IZLAZ,2100) N22,N11,KORA
      STOP
C
 1000 FORMAT(2I5,2F10.2,I5)
 1030 FORMAT(4I5)
 2020 FORMAT(8X,4I5,3X,2I6,5X,E12.5,E13.5,2X,I7)
C-----------------------------------------------------------------------
 2000 FORMAT(///'1'/6X,'D A T A    A B O U T    S U R F A C E    L O A D   
     1I N G   (PRESSURE)    A X I S Y M E T R I C    S H E L L S'/6X,
     1110('-'))
CS 2000 FORMAT('1'/6X,'P O D A C I    O    P O V R S I N S K I M    P R I
CS     1T I S C I M A    O S N O S I M E T R I C N E   LJ U S K E '/6X,
CS     1110('-'))
 2010 FORMAT(//////6X,'R E A D I N G    D A T A    '///13X,
     1'   N O D E S   ',6X,'NFUN',2X,'IPRAV',3X,
     1'PROPORTIONALITY FACTORS  ',5X,' STEP   ')
CS 2010 FORMAT(//////6X,'U C I T A N I    P O D A C I'///13X,
CS     1' C V O R O V I ',6X,'NFUN',2X,'IPRAV',3X,
CS     1'FAKTORI PROPORCIONALNOSTI',5X,'KOR.GEN.')
 2040 FORMAT(6X,'G E N E R A T E D    D A T A'///13X,
     1'   N O D E S   ',6X,'NFUN',2X,'IPRAV',3X,
     1'PROPORTIONALITY FACTORS  ')
CS 2040 FORMAT(6X,'G E N E R I S A N I    P O D A C I'///13X,
CS     1' C V O R O V I ',6X,'NFUN',
CS     12X,'IPRAV',3X,'FAKTORI PROPORCIONALNOSTI')
 2100 FORMAT(///' NODE NUMBER N2=',I5,' CAN NOT BE DEFINED BY SUPERPOSIT
     1ION OF NODE N1=',I5,' WITH A NUMBER OF STEPS  KORA=',I5)
CS 2100 FORMAT(///' BROJ N2=',I5,' NE MOZE SE DOBITI SABIRANJEM BROJA N1='
CS     1,I5,' I KONACNOG BROJA KORAKA KORA=',I5)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE TGRAF7(NEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C.    PROGRAM
C.       ZA STAMPANJE OSNOSIMETRICNE LJUSKE U UNIVERZALNI FILE
C.
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SUMELE/ ISUMEL,ISUMGR
      COMMON /NIDEAS/ IDEAS
C
      DIMENSION NEL(NE,*)
C
C     GRAFICKI OPIS STAPA: SA 2 CVORA = 1
      IF(ideas.eq.-1) return
      IFGD=1
C     
      IFDI=14
C     TABELA FIZICKIH OSOBINA
      IPTN=NGE
C     TABELA MATERIJALA
      MPTN=NMODM
C     BOJA  
      ICOL=8
C     BROJ CVOROVA NA ELEMENTU
      NNODS=2
      IND=-1
      ITYP=71
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 I=1,NE
C     REDNI BROJ ELEMENTA
      IEL=I+ISUMEL
      WRITE(IGRAF,1000) IEL,IFGD,IFDI,IPTN,MPTN,ICOL,NNODS
      WRITE(IGRAF,1000)(NEL(I,J),J=1,2)
   10 CONTINUE
      WRITE(IGRAF,1100) IND
      ISUMEL=ISUMEL+NE
      RETURN
 1000 FORMAT(8I10)
 1100 FORMAT(I6)
      END
