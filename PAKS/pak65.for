C=======================================================================
C
C        INTEGRALJENJE MATRICA TANKOZIDNE GREDE
C
C   SUBROUTINE TPGMAS 
C              READE6
C              SISTT6
C              ELTE6
C=======================================================================
      SUBROUTINE TPGMAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C    GLAVNI PROGRAM ZA POZIVANJE PROGRAMA ZA RACUNANJE MATRICA ELEMENATA
C
      include 'paka.inc'
      
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMA6/ MXAU,LAU,LLMEL,LNEL,LNMAT,LCTR,LMGLV,LLUV,LINDPR,
     1                LELKG,LNAPGV,LPRR,LNTIPG
      COMMON /DUPLAP/ IDVA
C
      LAU=LMAX
      CALL READ6M(A(LAU))
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LAE=LMAX
      NCVE6=NCVE*6
      NN6 = NCVE6*NCVE6
      MXAE = NCVE6+(2*NN6+NCVE6*(NCVE6+1)/2)*IDVA
      LMAX = LAE + MXAE
      IF(LMAX.LT.MTOT) GO TO 70
      WRITE(IZLAZ,2009) LMAX,MTOT
 2009 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA MATRICE ELEMENATA'
     1/' POTREBNA DIMENZIJA, LMAX=',I10/
     2' RASPOLOZIVA DIMENZIJA, MTOT=',I10)
      STOP
C
C     FORMIRANJE MATRICE KRUTOSTI ELEMENATA I PAKOVANJE U SISTEM
C
   70 CALL SIST6M(A(LAE),A(LAU))
C
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE READ6M(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
C
      include 'paka.inc'
      
      COMMON /ELEMA6/ MXAU,LAU,LLMEL,LNEL,LNMAT,LCTR,LMGLV,LLUV,LINDPR,
     1                LELKG,LNAPGV,LPRR,LNTIPG
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON/TRAN/X1,Y1,Z1,X2,Y2,Z2,DUZ,TRR(3,3),AKOR(4)
      COMMON /COMPRF/ YC,ZC,YI,ZI,ALFB,TK,CAPAY,CAPAZ,AREA,POVR,YM,ZM,
     1                DDELT,YP,SWO(1),OM(1),DELTAP(3),TAUTK(3),OY,OZ,
     2                NT,NELM,NTIP,NKAR,IND,INDRPR(3),KR1,NEXTR,NPRPR
C
      DIMENSION AU(*)
      REAL AU
C     POZIVANJE PROGRAMA ZA ULAZNE PODATKE .
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     1 NCVE,MXAU,LNEL,LNMAT,LNTIPG,LCTR,LMGLV,LLUV,LINDPR,LELKG,
     2 LNAPGV,LPRR,LCEL,LELC,NMA,NMI,
     3 YC,ZC,YI,ZI,ALFB,TK,CAPAY,CAPAZ,YM,ZM,AREA,POVR,DDELT,YP
     1 ,(DELTAP(I),I=1,3),(TAUTK(I),I=1,3),(AKOR(I),I=1,4),
     2 NT,NELM,NTIP,NKAR,KR1,NEXTR,NPRPR
C
      LSTAZA(2)=LMAX8+1
      CALL READDD(AU(LNEL),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAU+MXAU
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
C     PLASTICNOST
      RETURN
      END
C======================================================================
      SUBROUTINE SIST6M(AE,AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA(KCAL=1)
C     RACUNANJE NAPONA (KCAL=2)
C
      include 'paka.inc'
      
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEMA6/ MXAU,LAU,LLMEL,LNEL,LNMAT,LCTR,LMGLV,LLUV,LINDPR,
     1                LELKG,LNAPGV,LPRR,LNTIPG
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
C
      DIMENSION AE(*),AU(*)
      REAL AE,AU
C
C     REPERI U VEKTORU ELEMENATA AE
C
      NCVE6=NCVE*6
      NWE=NCVE6*(NCVE6+1)/2
      NN6 = NCVE6*NCVE6
      LSQ=1
      LSKG=LSQ+NN6*IDVA
      LSKE=LSKG+NN6*IDVA
      LLM=LSKE+NWE*IDVA
      MXAE1=LLM+NCVE6
      IF(MXAE1.EQ.MXAE+1)GO TO 30
      WRITE(IZLAZ,2000)MXAE1,MXAE
 2000 FORMAT(/' ','NE POKLAPAJU SE MXAE1 I MXAE'/2I10)
      STOP 'SISTT6'
C
C     OSNOVNA PETLJA PO ELEMENTIMA
C
   30 DO 100 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
CE          LTBTH: POINTER FOR TIME OF ELEMENT BIRTH
CE          LTDTH: POINTER FOR TIME OF ELEMENT DEATH
CE          VREME: CURRENT TIME
CE          NLM:   CURRENT ELEMENT
CE          IBD:   INDICATOR FOR EXISTENCE OF ELEMENT (=0-YES; =1-NO)
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 100
      IZBR=MXAE/IDVA
      CALL CLEAR(AE,IZBR)
C
C     RACUNANJE MATRICE ELEMENATA I/ILI NAPONA
C
      CALL ELT6M (A(LSK),A(LID),A(LCORD),AU(LNEL),AU(LNTIPG),
     1 AU(LCTR),AU(LELKG),AU(LNAPGV),AU(LMGLV),AU(LLUV),AE(LSQ),
     2 AE(LSKG),AE(LSKE),AE(LLM),AU(LLMEL),AU(LPRR),A(LGUSM),NLM,NCVE6)
C
C     PAKOVANJE MATRICA ELEMENATA U MATRICU SISTEMA
C
      IF(IMASS.EQ.1) CALL SPAKUJ(A(LSK),A(LMAXA),AE(LSKE),AE(LLM),NCVE6)
  100 CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE ELT6M(SK,ID,CORD,NEL,NTIPG,CTR,ELKG,NAPGV,MGLV,LUV,
     1 SQ,SKG,SKE,LM,LMEL,PR,GUSM,NLM,NCVE6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     OVAJ PODPROGRAM SLUZI ZA FORMIRANJE MATRICE KRUTOSTI 12X12
C     LINIJSKOG ELEMENTA TANKOZIDNE GREDE
C
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /REPERM/ MREPER(4)
      COMMON /TEMPCV/ LTECV,ITEMP
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
CCCC  COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /PAKDFL/ IND16,IR16,IND40,IND49,IND69,IND99
      COMMON /PGREDG/ ELCON(3),EE,PNE,G
      COMMON /PGREDC/ NPUN,NPRBR,MAXPRR,NPRES,ICV1,ICV2,MATG,
     1                NKARD,MG,NVR,NWRITE,NAPG,INDOF
      COMMON /COMPRF/ YC,ZC,YI,ZI,ALFB,TK,CAPAY,CAPAZ,AREA,POVR,YM,ZM,
     1                DDELT,YP,SWO(1),OM(1),DELTAP(3),TAUTK(3),OY,OZ,
     2                NT,NELM,NTIP,NKAR,IND,INDRPR(3),KR1,NEXTR,NPRPR
      COMMON/TRAN/X1,Y1,Z1,X2,Y2,Z2,DUZ,TRR(3,3),AKOR(4)
      DIMENSION SK(*),ID(NP,*),CORD(NP,*),NEL(NE,*),NTIPG(NE,*),
     1CTR(NE,*),ELKG(MATG,*),NAPGV(*),MGLV(*),LUV(NE,*)
C
      DIMENSION SQ(NCVE6,*),SKG(NCVE6,*),LM(*),
     1SKE(*),PR(*),GUSM(50,*),U(12),SKD(12)
C
C      SQRT(X)=DSQRT(X)
C      WRITE(3,*)'POVR,IY,IZ,TK,ALFB',POVR,YI,ZI,TK,ALFB
      NAP = NAPG
C
      NSZM = 0
      NVERT = 0
      NPUN = 0
C
C     OSNOVNA PETLJA PO ELEMENTIMA
C     ODREDITI BROJEVE $VOROVA I KARAKTERISTIKE ELEMENTA
      ICV1=NEL(NLM,1)
      ICV2=NEL(NLM,2)
C
      IF(NLM.EQ.1)TEZ=0.D0
C     RO=ELCON(3)
      X1=CORD(ICV1,1)
      Y1=CORD(ICV1,2)
      Z1=CORD(ICV1,3)
      X2=CORD(ICV2,1)
      Y2=CORD(ICV2,2)
      Z2=CORD(ICV2,3)
      X12 = X1 - X2
      Y12 = Y1 - Y2
      Z12 = Z1 - Z2
      DUZ=DSQRT(X12*X12 + Y12*Y12 + Z12*Z12)
      DUZ2=DUZ*DUZ
C
      MAT = 1
      NMODM=1
      GO TO (1,991,991,991,991,991,991,991),NMODM
    1 LFUN=MREPER(1)
  991 CONTINUE
C***
      AAA=POVR
      PX=TK
      PZ=YI
      PY=ZI
      RO=GUSM(NMODM,MAT)
C
      TEZ1=AAA*DUZ*RO
C     UNULITI MATRICU SQ(12*12)
      CALL CLEAR(U,12)
      CALL CLEAR(SQ,144)
      TRSP=13.D0/35.
      SPT =6.D0/5.
      A11D=11.D0/210.
      DVSE=9.D0/70.
      TRCD=13.D0/420.
      DVPT=2.D0/15. 
      DUZ3 = DUZ2*DUZ
      SQ(1,1)=0.5D0
      SQ(1,7)=0.D0
C      SQ(1,1)=0.33333333D0
C      SQ(1,1)=5.D0/12.
C      SQ(1,7)=SQ(1,1)/2.
C      SQ(1,7)=1.D0/12.
      SQ(7,7)=SQ(1,1)
      SQ(2,2)=TRSP+SPT*PZ/(AAA*DUZ2)
      SQ(2,6)=A11D*DUZ+PZ/(10.*AAA*DUZ)
      SQ(2,8)=DVSE-SPT*PZ/(AAA*DUZ2)
      SQ(2,12)=PZ/(10.*AAA*DUZ)-TRCD*DUZ
      SQ(3,3)=TRSP+SPT*PY/(AAA*DUZ2)
      SQ(3,5)=-A11D*DUZ-PY/(10.*AAA*DUZ)
      SQ(3,9)=DVSE-SPT*PY/(AAA*DUZ2)
      SQ(3,11)=-PY/(10.*AAA*DUZ)+TRCD*DUZ
      SQ(4,4) =PX/AAA/3.
      SQ(4,10)=SQ(4,4)/2.
      SQ(10,10)=SQ(4,4)
      DSP = DUZ2/105.
      SQ(5,5)=DSP+DVPT*PY/AAA
      SQ(5,9)=-SQ(3,11)
      SQ(5,11)=-DUZ2/140.-PY/AAA/30.
      SQ(6,6)=DSP+DVPT*PZ/AAA
      SQ(6,8)=-SQ(2,12)
      SQ(6,12)=-DUZ2/140.-PZ/AAA/30.
      SQ(8,8)=SQ(2,2)
      SQ(8,12)=-SQ(2,6)
      SQ(9,9)=SQ(3,3)
      SQ(9,11)=-SQ(3,5)
      SQ(11,11)=SQ(5,5)
      SQ(12,12)=SQ(6,6)
C     IZJEDNACENJE SIMETRICNIH CLANOVA SQ(I,J)
      DO 71 I=1,12
      DO 71 J=1,12
      IF(I-J)72,71,71
   72 SQ(J,I)=SQ(I,J)
   71 CONTINUE
C     MNOZENJE MATRICE SQ(I,J) KONSTANTOM EE
      DO 73 I=1,12
      DO 73 J=1,12
   73 SQ(I,J)=TEZ1*SQ(I,J)
C***
C  KOREKCIJA MATRICE SQ AKO IMA SLOBODNIH POMERANJA
C  ODREDIVANJE MATRICE TRANSFORMACIJE TR(12*12)
C     UNULITI MATRICU TR(3,3)
C
      CALL TRANG (CTR,NE,U,NLM,1)
C
      DO 100 I=1,6
      LM(I)=ID(ICV1,I)
  100 LM(I+6)=ID(ICV2,I)
C
      KK=0
      DO 10 II=1,4
      DO 11 LL=1,3
      DO 11 JJ=1,12
      SKG(LL+KK,JJ)=0.D0
      DO 11 K1=1,3
   11 SKG(LL+KK,JJ)=SKG(LL+KK,JJ)+TRR(K1,LL)*SQ(K1+KK,JJ)
      KK=KK+3
   10 CONTINUE
      KK=0
      DO 83 II=1,12
      DO 83 JJ=1,12
      SQ(II,JJ)=SKG(II,JJ)
   83 SKG(II,JJ)=0.D0
      DO 12 II=1,4
      DO 13 LL=1,3
      DO 13 JJ=1,12
      DO 13 K1=1,3
   13 SKG(JJ,LL+KK)=SKG(JJ,LL+KK)+SQ(JJ,K1+KK)*TRR(K1,LL)
      KK=KK+3
   12 CONTINUE
C
C  FORMIRANJE MATRICE KRUTOSTI U OBLIKU VEKTIRA-SYM. DEO
C
      II=0
      DO 101 I=1,12
      DO 101 J=I,12
      II=II+1
      IF(I.EQ.J) SKD(I)=SKG(I,I)
  101 SKE(II)=SKG(I,J)
C
      IF(IMASS.EQ.2) THEN
         CALL CLEAR(SKD,12)
         SKD(1)=.5*TEZ1
         SKD(2)=SKD(1)
         SKD(3)=SKD(1)
         SKD(7)=SKD(1)
         SKD(8)=SKD(1)
         SKD(9)=SKD(1)
         DO 102 I=1,12
            II=LM(I)
            IF(II.GT.0) SK(II)=SK(II)+SKD(I)
  102    CONTINUE
      ENDIF
C
      RETURN
      END