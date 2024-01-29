MODULE math
!-----------------------------------------------------------------------
USE parameters
!-----------------------------------------------------------------------
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(0:301) :: FACT
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

!=======================================================================
SUBROUTINE factrl
!-----------------------------------------------------------------------
  INTEGER :: i
!-----------------------------------------------------------------------
! Calculates the first 301 factorials
  fact(0)=1.d0
  do i=1,301
    fact(i)=fact(i-1)*dfloat(i)
  enddo
!-----------------------------------------------------------------------
END SUBROUTINE factrl
!=======================================================================


!=======================================================================
SUBROUTINE comprod(a,b,c,d,e,f)
!-----------------------------------------------------------------------
! PERFORMS THE PRODUCT OF TWO COMPLEX NUMBERS: (A+iB)*(C+iD)=(E+iF)
! ON INPUT:   A,B,C,D
! ON OUTPUT:  E,F
!-----------------------------------------------------------------------
  DOUBLE PRECISION, INTENT(IN) :: a, b, c, d
  DOUBLE PRECISION, INTENT(OUT) :: e, f
!-----------------------------------------------------------------------
  e=a*c-b*d
  f=a*d+b*c
  RETURN
!-----------------------------------------------------------------------
END SUBROUTINE comprod
!=======================================================================


!=======================================================================
FUNCTION W3JS(J1,J2,J3,M1,M2,M3)
! Calculates 3J symbols. Input arguments must be multiplied by 2
!-----------------------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER, INTENT(IN) :: J1, J2, J3, M1, M2, M3
  INTEGER :: IA,IB,IC,ID,IE,IF1,IG,IH, JSUM
  INTEGER :: Z,ZMIN,ZMAX
  call factrl 
  W3JS=0.0
  IF(M1+M2+M3.NE.0) GOTO 1000
  IA=J1+J2
  IF(J3.GT.IA) GOTO 1000
  IB=J1-J2
  IF(J3.LT.IABS(IB)) GOTO 1000
  JSUM=J3+IA
  IC=J1-M1
  ID=J2-M2
  IF(MOD(JSUM,2).NE.0) GOTO 1000
  IF(MOD(IC,2).NE.0) GOTO 1000
  IF(MOD(ID,2).NE.0) GOTO 1000
  IF(IABS(M1).GT.J1) GOTO 1000
  IF(IABS(M2).GT.J2) GOTO 1000
  IF(IABS(M3).GT.J3) GOTO 1000
  IE=J3-J2+M1
  IF1=J3-J1-M2
  ZMIN=MAX0(0,-IE,-IF1)
  IG=IA-J3
  IH=J2+M2
  ZMAX=MIN0(IG,IH,IC)
  CC=0.0
  DO 200 Z=ZMIN,ZMAX,2
   DENOM=FACT(Z/2)*FACT((IG-Z)/2)*FACT((IC-Z)/2) &
        *FACT((IH-Z)/2)*FACT((IE+Z)/2)*FACT((IF1+Z)/2)
   IF(MOD(Z,4).NE.0) DENOM=-DENOM
   CC=CC+1.0/DENOM
  200 CONTINUE
  CC1=FACT(IG/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2) &
   *1/FACT((JSUM+2)/2)
    CC2=FACT((J1+M1)/2)*FACT(IC/2)*FACT(IH/2) &
   *FACT(ID/2)*FACT((J3-M3)/2)*FACT((J3+M3)/2)
    CC=CC*DSQRT(CC1*CC2)
    IF(MOD(IB-M3,4).NE.0) CC=-CC
    W3JS=CC
    IF(ABS(W3JS).LT.1.0E-8) W3JS=0.0
 1000 RETURN
END
!=======================================================================


!=======================================================================
FUNCTION W6JS(J1,J2,J3,L1,L2,L3)
!-----------------------------------------------------------------------
! Calculates 6J symbols. Input arguments must be multiplied by 2
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: J1, J2, J3, L1, L2, L3
  INTEGER :: IA, IB, IC, ID, IE, IIF, IG, IH, II, IJ, IK
  INTEGER :: SUM1,SUM2,SUM3,SUM4
  INTEGER :: W,WMIN,WMAX
  DOUBLE PRECISION :: OMEGA, DENOM
  DOUBLE PRECISION :: THETA, THETA1, THETA2, THETA3, THETA4
  DOUBLE PRECISION :: W6JS
!-----------------------------------------------------------------------
call factrl 
 W6JS=0.0
  IA=J1+J2
  IF(IA.LT.J3) RETURN
  IB=J1-J2
  IF(IABS(IB).GT.J3) RETURN
  IC=J1+L2
  IF(IC.LT.L3) RETURN
  ID=J1-L2
  IF(IABS(ID).GT.L3) RETURN
  IE=L1+J2
  IF(IE.LT.L3) RETURN
  IIF=L1-J2
  IF(IABS(IIF).GT.L3) RETURN
  IG=L1+L2
  IF(IG.LT.J3) RETURN
  IH=L1-L2
  IF(IABS(IH).GT.J3) RETURN
  SUM1=IA+J3
  SUM2=IC+L3
  SUM3=IE+L3
  SUM4=IG+J3
  WMIN=MAX0(SUM1,SUM2,SUM3,SUM4)
  II=IA+IG
  IJ=J2+J3+L2+L3
  IK=J3+J1+L3+L1
  WMAX=MIN0(II,IJ,IK)
  OMEGA=0.D0
  DO W=WMIN,WMAX,2
    DENOM=FACT((W-SUM1)/2)*FACT((W-SUM2)/2)*FACT((W-SUM3)/2)&
         *FACT((W-SUM4)/2)*FACT((II-W)/2)*FACT((IJ-W)/2)*FACT((IK-W)/2)
    IF(MOD(W,4).NE.0) DENOM=-DENOM
    OMEGA=OMEGA+FACT(W/2+1)/DENOM
  ENDDO
  THETA1=FACT((IA-J3)/2)*FACT((J3+IB)/2)*FACT((J3-IB)/2)/FACT(SUM1/2+1)
  THETA2=FACT((IC-L3)/2)*FACT((L3+ID)/2)*FACT((L3-ID)/2)/FACT(SUM2/2+1)
  THETA3=FACT((IE-L3)/2)*FACT((L3+IIF)/2)*FACT((L3-IIF)/2)/FACT(SUM3/2+1)
  THETA4=FACT((IG-J3)/2)*FACT((J3+IH)/2)*FACT((J3-IH)/2)/FACT(SUM4/2+1)
  THETA=THETA1*THETA2*THETA3*THETA4
  W6JS=OMEGA*DSQRT(THETA)
  IF(ABS(W6JS).LT.1.0D-8) W6JS=0.0
!-----------------------------------------------------------------------
  RETURN
!-----------------------------------------------------------------------
END FUNCTION W6JS
!=======================================================================
 
!=======================================================================
FUNCTION wul(ju2,jl2)
!-----------------------------------------------------------------------
 DOUBLE PRECISION :: wul
 INTEGER, INTENT(IN) :: jl2, ju2
 INTEGER :: add
!-----------------------------------------------------------------------
  wul= dsqrt(3*dble(ju2+1))*w6js(2,2,4,ju2,ju2,jl2)
  add = ju2+jl2
  if(mod(add,2).ne.0) then !Difference Ju-Jl must be 0,\pm 1
   wul = 0.d0
   print*, 'The polarizability factor is undefined for this pair of J_u and J_l'
  elseif(mod(add,4).eq.0) then
   wul = -wul
  endif
!-----------------------------------------------------------------------
return
!-----------------------------------------------------------------------
END FUNCTION wul
!=======================================================================

!=======================================================================
SUBROUTINE ROTMATSU(KMAX,ALFA,BETA,GAMMA,DR,DI)
!-----------------------------------------------------------------------
!CALCULATES THE ROTATION MATRICES OF ORDER 0, 1, 2,....,KMAX
!THE RESULTS ARE GIVEN IN THE MATRICES DR(K,IQ,IQP) AND 
!DI(K,IQ,IQP) REPRESENTING THE REAL AND IMAGINARY PART OF THE 
!OTATION MATRIX, RESPECTIVELY.
!THE ANGLES ALFA, BETA, AND GAMMA ARE GIVEN IN RADIANS.
!THE MAXIMUM ALLOWED VALUE FOR KMAX IS 6
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: KMAX
 DOUBLE PRECISION, INTENT(IN) :: ALFA,BETA,GAMMA  
 DOUBLE PRECISION, DIMENSION(0:6,-6:6,-6:6), INTENT(OUT) :: DR, DI
 INTEGER :: K, IQ, IQP,IT1,IT2,ITMAX,ITMIN,IT3,IT,IE1,IE2
 DOUBLE PRECISION :: A1,A2,CA,SA,CG,SG,BO2,CB2,SB2
 DOUBLE PRECISION :: REDM,S,DEN
!-----------------------------------------------------------------------
CALL FACTRL
DO 1 K=0,6
 DO 1 IQ=-6,6
  DO 1 IQP=-6,6
   DR(K,IQ,IQP)=0.
   DI(K,IQ,IQP)=0.
1    CONTINUE
 DO 2 K=0,KMAX
  DO 3 IQ=-K,K
   DO 4 IQP=-K,K
    A1=ALFA*DFLOAT(IQ)
    A2=GAMMA*DFLOAT(IQP)
    CA=DCOS(A1)
    SA=DSIN(A1)
    CG=DCOS(A2)
    SG=DSIN(A2)
    BO2=BETA/2.
    CB2=DCOS(BO2)
    SB2=DSIN(BO2)
    IT1=K+IQ
    IT2=K-IQP
    ITMAX=IT1
    IF(IT2.LT.ITMAX) ITMAX=IT2
    ITMIN=0
    IT3=IQ-IQP
    IF(IT3.GT.0) ITMIN=IT3
    REDM=0.D0
    IF(BETA.EQ.0.D0) GO TO 6
    DO 5 IT=ITMIN,ITMAX
     S=1.D0
     IF(MOD(IT,2).NE.0) S=-1.D0
     IE1=2*(K-IT)+IQ-IQP
     IE2=2*IT+IQP-IQ
     DEN=FACT(IT1-IT)*FACT(IT2-IT)*FACT(IT)*FACT(IT-IT3)
     REDM=REDM+S*(CB2**IE1)*(SB2**IE2)/DEN
5    CONTINUE
    REDM=REDM*DSQRT(FACT(K+IQ)*FACT(K-IQ)*FACT(K+IQP)*FACT(K-IQP))
    GO TO 7
 6    CONTINUE
    IF(IQ.EQ.IQP) REDM=1.
 7    CONTINUE
    DR(K,IQ,IQP)=REDM*(CA*CG-SA*SG)
    DI(K,IQ,IQP)=-REDM*(CA*SG+SA*CG)
 4    CONTINUE
 3    CONTINUE
 2    CONTINUE
!-----------------------------------------------------------------------
RETURN
!-----------------------------------------------------------------------
END SUBROUTINE ROTMATSU
!=======================================================================

!=======================================================================
SUBROUTINE TKQ2(ALFA,BETA,GAMMA,ALFAP,BETAP,GAMMAP,TR,TI)
!-----------------------------------------------------------------------
!CALCULATES THE TENSOR TKQ FOR A DOUBLE ROTATION
!ON INPUT: ALFA,BETA,GAMMA=EULER ANGLES IN RADIANS OF FIRST ROT.
!        : ALFAP,BETAP,GAMMAP=SAME FOR SECOND ROTATION
!ON OUTPUT: TR, TI= REAL AND IMAG PARTS OF TENSOR TKQ
!-----------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: alfa, beta, gamma,alfap,betap,gammap
DOUBLE PRECISION, INTENT(OUT) :: TR(0:2,-2:2,0:3),TI(0:2,-2:2,0:3)
DOUBLE PRECISION, DIMENSION(0:6,-6:6,-6:6) :: DR1, DI1, DR2, DI2
DOUBLE PRECISION, DIMENSION(0:3) :: TPR,TPI
INTEGER :: i,k,iq,ip,ir
DOUBLE PRECISION :: e,f,g,h, r,s
!-----------------------------------------------------------------------
CALL ROTMATSU(2,ALFA,BETA,GAMMA,DR1,DI1)
CALL ROTMATSU(2,ALFAP,BETAP,GAMMAP,DR2,DI2)
DO i=0,3
 DO k=0,2
  DO iq=-k,k
   TR(k,iq,i)=0.d0
   TI(k,iq,i)=0.d0
   r=0.d0
   s=0.d0
   DO ip=-k,k
    CALL TKP(k,ip,TPR,TPI)
    DO ir=-k,k
     CALL COMPROD(TPR(i),TPI(i),DR1(k,ip,ir),DI1(k,ip,ir),e,f)
     CALL COMPROD(e,f,DR2(k,ir,iq),DI2(k,ir,iq),g,h)
     r=r+g 
     s=s+h
    ENDDO !IR
   ENDDO! IP
   TR(k,iq,i) = r
   TI(k,iq,i) = s
  ENDDO! IQ
 ENDDO! K
ENDDO! I
!------------------------------------------------------------------------
 RETURN
!------------------------------------------------------------------------
END SUBROUTINE TKQ2
!=======================================================================


!=======================================================================
 SUBROUTINE TKP(k,ip,TR,TI)
!------------------------------------------------------------------------
!-----CALCULATES THE SYMBOL TKP------------------------
 INTEGER, INTENT(IN) ::k, ip
 INTEGER :: i
 DOUBLE PRECISION :: A
 DOUBLE PRECISION, INTENT(OUT) :: TR(0:3), TI(0:3)
!------------------------------------------------------------------------
DO i = 0,3
 TR(I)=0.d0
 TI(I)=0.d0
ENDDO !i
IF(k.EQ.0.AND.ip.EQ.0) TR(0) = 1.d0
IF(k.EQ.1.AND.ip.EQ.0) TR(3) = DSQRT(1.5d0)
IF(k.EQ.2.AND.ip.EQ.0) TR(0) = DSQRT(0.5d0)
IF(k.EQ.2.AND.ip.EQ.2) THEN
A=DSQRT(3.D0)/2.D0
TR(1)=-A
TI(2)=-A    
ENDIF
IF(k.EQ.2.AND.ip.EQ.-2) THEN
 A=DSQRT(3.D0)/2.D0
 TR(1)=-A
 TI(2)=A    
ENDIF
!------------------------------------------------------------------------
 RETURN
!------------------------------------------------------------------------
END SUBROUTINE TKP
!=======================================================================


!=======================================================================
SUBROUTINE wav_vac_air(wvac,wair)
!-----------------------------------------------------------------------
!  implicit double precision (a-h,o-z)
  INTEGER :: ir, i1, i2
  DOUBLE PRECISION, INTENT(IN) :: wvac
  DOUBLE PRECISION, INTENT(OUT) :: wair
  DOUBLE PRECISION :: r, p, d, w1
  DOUBLE PRECISION, dimension(0:90) :: delta
!-----------------------------------------------------------------------
  data delta/0.648d0,0.667d0,0.687d0,0.708d0,0.731d0,&
             0.754d0,0.777d0,0.801d0,0.825d0,0.850d0,&
             0.875d0,0.900d0,0.925d0,0.950d0,0.976d0,&
             1.001d0,1.027d0,1.053d0,1.079d0,1.105d0,&
             1.131d0,1.157d0,1.183d0,1.210d0,1.236d0,&
             1.262d0,1.289d0,1.315d0,1.342d0,1.368d0,&
             1.395d0,1.421d0,1.448d0,1.475d0,1.501d0,&
             1.528d0,1.555d0,1.581d0,1.608d0,1.635d0,&
             1.662d0,1.689d0,1.715d0,1.742d0,1.769d0,&
             1.796d0,1.823d0,1.850d0,1.877d0,1.904d0,&
             1.931d0,1.957d0,1.984d0,2.011d0,2.038d0,&
             2.065d0,2.092d0,2.119d0,2.146d0,2.173d0,&
             2.200d0,2.227d0,2.254d0,2.281d0,2.308d0,&
             2.335d0,2.362d0,2.389d0,2.417d0,2.444d0,&
             2.471d0,2.498d0,2.525d0,2.552d0,2.579d0,&
             2.606d0,2.633d0,2.660d0,2.687d0,2.714d0,&
             2.741d0,2.769d0,2.796d0,2.823d0,2.850d0,&
             2.877d0,2.904d0,2.931d0,2.958d0,2.985d0,&
             3.012d0/
!-----------------------------------------------------------------------
  r=(wvac-2000.d0)/100.d0
  ir=int(r)
  i1=ir
  if(wvac.lt.2000.d0) i1=0
  if(wvac.gt.11000.d0) i1=89
  i2=i1+1
  w1=2000.d0+100.d0*dfloat(i1)
  p=(wvac-w1)/100.d0
  d=(1.d0-p)*delta(i1)+p*delta(i2)
  wair=wvac-d
  return
!-----------------------------------------------------------------------
END SUBROUTINE wav_vac_air
!=======================================================================

!=======================================================================
SUBROUTINE wav_air_vac(wair,wvac)
!-----------------------------------------------------------------------
  INTEGER :: ir, i1, i2
  DOUBLE PRECISION, INTENT(IN) :: wair
  DOUBLE PRECISION, INTENT(OUT) :: wvac
  DOUBLE PRECISION :: r, p, d, w1
  DOUBLE PRECISION, dimension(0:90) :: delta
!-----------------------------------------------------------------------
  data delta/0.648d0,0.667d0,0.687d0,0.708d0,0.731d0,&
             0.754d0,0.777d0,0.801d0,0.825d0,0.850d0,&
             0.875d0,0.900d0,0.925d0,0.950d0,0.976d0,&
             1.001d0,1.027d0,1.053d0,1.079d0,1.105d0,&
             1.131d0,1.157d0,1.183d0,1.210d0,1.236d0,&
             1.262d0,1.289d0,1.315d0,1.342d0,1.368d0,&
             1.395d0,1.421d0,1.448d0,1.475d0,1.501d0,&
             1.528d0,1.555d0,1.581d0,1.608d0,1.635d0,&
             1.662d0,1.689d0,1.715d0,1.742d0,1.769d0,&
             1.796d0,1.823d0,1.850d0,1.877d0,1.904d0,&
             1.931d0,1.957d0,1.984d0,2.011d0,2.038d0,&
             2.065d0,2.092d0,2.119d0,2.146d0,2.173d0,&
             2.200d0,2.227d0,2.254d0,2.281d0,2.308d0,&
             2.335d0,2.362d0,2.389d0,2.417d0,2.444d0,&
             2.471d0,2.498d0,2.525d0,2.552d0,2.579d0,&
             2.606d0,2.633d0,2.660d0,2.687d0,2.714d0,&
             2.741d0,2.769d0,2.796d0,2.823d0,2.850d0,&
             2.877d0,2.904d0,2.931d0,2.958d0,2.985d0,&
             3.012d0/
!-----------------------------------------------------------------------
  r=(wair-2000.d0)/100.d0
  ir=int(r)
  i1=ir
  if(wair.lt.2000.d0) i1=0
  if(wair.gt.11000.d0) i1=89
  i2=i1+1
  w1=2000.d0+100.d0*dfloat(i1)
  p=(wair-w1)/100.d0
  d=(1.d0-p)*delta(i1)+p*delta(i2)
  wvac=wair+d
  return
!-----------------------------------------------------------------------
END SUBROUTINE wav_air_vac
!=======================================================================


!=======================================================================
SUBROUTINE LGDR(N,A,B,X,W)
!*************************************************************
!     Subroutine for Gauss-Legendre Quadrature               *
!     (in OUTPUT: Weights, Abscissae)                        *
!     Veronique Bommier, Meudon, February 1993               *
!*************************************************************
!----- Integral from A to B of F(X)dX = 
!----- Sum from I=1 to N of W(I)*F(X(I)); A < B
!----- INPUT:  N: Number of points, A, B
!----- OUTPUT: X(I=1,...,N): Abscissae
!-----         W(I=1,...,N): Weights
!----- Basis change in this Subroutine: Z=(X-(a+b)/2)*2/(a-b)
!----- without basis change: a=-1. and b=+1.
!-----     (in this case, a basis change remains however: Z=-X)
!----- This Subroutine: "Numerical Recipes", P.125
!-----   | Recurrence and Derivative Formulae:
!-----   | Abramowitz and Stegun P.782-783
!-----   | Weights: "Numerical Analysis", Z.Kopal, P.369
!-----   | Approximation for Zeros of Legendre Polynomials:
!-----   | Abramowitz and Stegun, P.787
!----- Possible TEST: The Weights sum is equal to b-a
!-----------------------------------------------------------------------
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      DIMENSION X(N),W(N)
!      DOUBLE PRECISION :: Z, PI
!      INTRINSIC COS,ACOS
!      PARAMETER (EPS=3.D-16)
!      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: N
  INTEGER :: I, J, M
  DOUBLE PRECISION, INTENT(IN) :: A, B
  DOUBLE PRECISION, INTENT(OUT) :: X(N), W(N)
  DOUBLE PRECISION :: Z, PI, XM, XL, P1, P2, P3, PP, Z1
  DOUBLE PRECISION, PARAMETER :: EPS=3.D-16
  DOUBLE PRECISION, PARAMETER :: ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0
!-----------------------------------------------------------------------
      PI=ACOS(-ONE)
      M=(N+1)/2
      XM=(A+B)/TWO
      XL=(B-A)/TWO
      DO 12 I=1,M
      Z=COS(PI*(FOUR*I-ONE)/(FOUR*N+TWO))
    1 CONTINUE
      P1=ONE
      P2=ZERO
      DO 11 J=1,N
      P3=P2
      P2=P1
      P1=((TWO*J-ONE)*Z*P2-(J-ONE)*P3)/J
   11 CONTINUE
      PP=N*(Z*P1-P2)/(Z*Z-ONE)
      Z1=Z
      Z=Z1-P1/PP
      IF(ABS(P1/PP).GT.EPS) GO TO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=TWO*XL/(ONE-Z*Z)/PP/PP
      W(N+1-I)=W(I)
   12 CONTINUE
      RETURN
!-----------------------------------------------------------------------
END SUBROUTINE LGDR
!=======================================================================

!=======================================================================
FUNCTION planck(nu0,T)
!-------------------------------------------------------------------------
! Returns the Planck function for a given temperature and frequency
!-------------------------------------------------------------------------
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: nu0, T
  DOUBLE PRECISION :: planck
  DOUBLE PRECISION :: a, b, esp
!-------------------------------------------------------------------------
  a=(2.d0*ph*nu0**3)/(pc**2)
  b=(ph*nu0)/(pk*T)
  esp=dexp(b)
  planck=a/(esp-1.d0)
!-------------------------------------------------------------------------
  return
!------------------------------------------------------------------------
END FUNCTION planck
!=======================================================================

!=======================================================================
FUNCTION wien(nu0,T)
!-------------------------------------------------------------------------
! Returns the Planck function in the Wien limit for a given temperature
!  and frequency 
!-------------------------------------------------------------------------
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: nu0, T
  DOUBLE PRECISION :: wien
  DOUBLE PRECISION :: a, b, esp
!-------------------------------------------------------------------------
  a=(2.d0*ph*nu0**3)/(pc**2)
  b=-1.d0*(ph*nu0)/(pk*T)
  esp=dexp(b)
  wien=a*esp
!-------------------------------------------------------------------------
  return
END FUNCTION wien
!=======================================================================

!=======================================================================
FUNCTION erfc(x)
!-------------------------------------------------------------------------
! Computing the Error function 
!-------------------------------------------------------------------------
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x
  DOUBLE PRECISION :: erfc
  DOUBLE PRECISION, PARAMETER :: a1=0.278393, a2=0.230389, a3=0.000972, a4=0.078108
!-------------------------------------------------------------------------
  erfc=1.d0/(1.d0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)**4
!-------------------------------------------------------------------------
  return
!-------------------------------------------------------------------------
END FUNCTION erfc
!=======================================================================


!=======================================================================
SUBROUTINE lubksb(a,n,indx,b)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: n, indx(n)
  DOUBLE PRECISION, INTENT(IN) :: a(n,n)
  DOUBLE PRECISION, INTENT(INOUT) :: b(n)
  INTEGER :: i, ii, j, ll
  DOUBLE PRECISION :: sum
!-----------------------------------------------------------------------
  ii=0
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if(ii.ne.0) then
      do j=ii,i-1
        sum=sum-a(i,j)*b(j)
      enddo
      else if(sum.ne.0.) then
      ii=i
    endif
    b(i)=sum
  enddo
  do i=n,1,-1
    sum=b(i)
    if(i.lt.n) then
      do j=i+1,n
        sum=sum-a(i,j)*b(j)
      enddo
    endif
    b(i)=sum/a(i,i)
  enddo
  return
!-----------------------------------------------------------------------
END SUBROUTINE LUBKSB
!=======================================================================


!=======================================================================
SUBROUTINE ludcmp(a,n,indx,d)
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: indx(n)
  DOUBLE PRECISION, INTENT(INOUT) :: a(n,n)
  DOUBLE PRECISION, INTENT(OUT) :: d
  INTEGER, PARAMETER :: nmax=3200
  DOUBLE PRECISION, PARAMETER :: tiny=1.d-20
  INTEGER :: i, imax, j, k
  DOUBLE PRECISION :: aamax, dum, sum, vv(nmax)
!-----------------------------------------------------------------------
  d=1.d0
  do i=1,n
    aamax=0.d0
    do j=1,n
      if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
    enddo
    if (aamax.eq.0.) STOP 'Singular matrix.'
    vv(i)=1.d0/aamax
  enddo  
  do j=1,n
    if (j.gt.1) then
      do i=1,j-1
        sum=a(i,j)
        if (i.gt.1)then
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        endif
      enddo
    endif
    aamax=0.d0
    do i=j,n
      sum=a(i,j)
      if (j.gt.1) then
        do k=1,j-1
          sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
      endif
      dum=vv(i)*dabs(sum)
      if (dum.ge.aamax) then
        imax=i
        aamax=dum
      endif
    enddo
    if (j.ne.imax) then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(j.ne.n) then
      if(a(j,j).eq.0.) a(j,j)=tiny
      dum=1.d0/a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      enddo
    endif
  enddo
  if(a(n,n).eq.0.) a(n,n)=tiny
  RETURN
!-----------------------------------------------------------------------
END SUBROUTINE LUDCMP
!=======================================================================


!=======================================================================
SUBROUTINE MatInv(a)
!------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(INOUT) :: a(0:3,0:3)
DOUBLE PRECISION :: absmax, det
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: b
INTEGER :: i,j
DOUBLE PRECISION, PARAMETER :: snap = 1.d-35
!------------------------------------------------------------------------
absmax = 0.d0
do i = 0,3
  do j = 0,3
  if(absmax .le. DABS(a(i,j))) absmax = a(i,j)
  enddo !j
enddo!i
!
if(DABS(absmax).le.1.d-100) stop
absmax = 1.d0/ absmax
!
do i= 0,3
 do j= 0,3
 a(i,j) = a(i,j)*absmax
 enddo!j
enddo!i
!
!Cofactors (reduced determinants)
b(0,0) = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
        -a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(2,2)*a(3,1) 
b(1,0) = a(1,2)*a(2,0)*a(3,3)+a(1,3)*a(2,2)*a(3,0)+a(1,0)*a(2,3)*a(3,2) &
        -a(1,2)*a(2,3)*a(3,0)-a(1,3)*a(2,0)*a(3,2)-a(1,0)*a(2,2)*a(3,3)
b(2,0) = a(1,3)*a(2,0)*a(3,1)+a(1,0)*a(2,1)*a(3,3)+a(1,1)*a(2,3)*a(3,0) &
        -a(1,3)*a(2,1)*a(3,0)-a(1,0)*a(2,3)*a(3,1)-a(1,1)*a(2,0)*a(3,3)
b(3,0) = a(1,0)*a(2,2)*a(3,1)+a(1,1)*a(2,0)*a(3,2)+a(1,2)*a(2,1)*a(3,0) &
        -a(1,0)*a(2,1)*a(3,2)-a(1,1)*a(2,2)*a(3,0)-a(1,2)*a(2,0)*a(3,1)
b(0,1) = a(2,1)*a(3,3)*a(0,2)+a(2,2)*a(3,1)*a(0,3)+a(2,3)*a(3,2)*a(0,1) &
        -a(2,1)*a(3,2)*a(0,3)-a(2,2)*a(3,3)*a(0,1)-a(2,3)*a(3,1)*a(0,2)
b(1,1) = a(2,2)*a(3,3)*a(0,0)+a(2,3)*a(3,0)*a(0,2)+a(2,0)*a(3,2)*a(0,3) &
        -a(2,2)*a(3,0)*a(0,3)-a(2,3)*a(3,2)*a(0,0)-a(2,0)*a(3,3)*a(0,2)
b(2,1) = a(2,3)*a(3,1)*a(0,0)+a(2,0)*a(3,3)*a(0,1)+a(2,1)*a(3,0)*a(0,3) &
        -a(2,3)*a(3,0)*a(0,1)-a(2,0)*a(3,1)*a(0,3)-a(2,1)*a(3,3)*a(0,0)
b(3,1) = a(2,0)*a(3,1)*a(0,2)+a(2,1)*a(3,2)*a(0,0)+a(2,2)*a(3,0)*a(0,1) &
        -a(2,0)*a(3,2)*a(0,1)-a(2,1)*a(3,0)*a(0,2)-a(2,2)*a(3,1)*a(0,0)
b(0,2) = a(3,1)*a(0,2)*a(1,3)+a(3,2)*a(0,3)*a(1,1)+a(3,3)*a(0,1)*a(1,2) &
        -a(3,1)*a(0,3)*a(1,2)-a(3,2)*a(0,1)*a(1,3)-a(3,3)*a(0,2)*a(1,1)
b(1,2) = a(3,2)*a(0,0)*a(1,3)+a(3,3)*a(0,2)*a(1,0)+a(3,0)*a(0,3)*a(1,2) &
        -a(3,2)*a(0,3)*a(1,0)-a(3,3)*a(0,0)*a(1,2)-a(3,0)*a(0,2)*a(1,3)
b(2,2) = a(3,3)*a(0,0)*a(1,1)+a(3,0)*a(0,1)*a(1,3)+a(3,1)*a(0,3)*a(1,0) &
        -a(3,3)*a(0,1)*a(1,0)-a(3,0)*a(0,3)*a(1,1)-a(3,1)*a(0,0)*a(1,3)
b(3,2) = a(3,0)*a(0,2)*a(1,1)+a(3,1)*a(0,0)*a(1,2)+a(3,2)*a(0,1)*a(1,0) &
        -a(3,0)*a(0,1)*a(1,2)-a(3,1)*a(0,2)*a(1,0)-a(3,2)*a(0,0)*a(1,1)
b(0,3) = a(0,1)*a(1,3)*a(2,2)+a(0,2)*a(1,1)*a(2,3)+a(0,3)*a(1,2)*a(2,1) &
        -a(0,1)*a(1,2)*a(2,3)-a(0,2)*a(1,3)*a(2,1)-a(0,3)*a(1,1)*a(2,2)
b(1,3) = a(0,2)*a(1,3)*a(2,0)+a(0,3)*a(1,0)*a(2,2)+a(0,0)*a(1,2)*a(2,3) &
        -a(0,2)*a(1,0)*a(2,3)-a(0,3)*a(1,2)*a(2,0)-a(0,0)*a(1,3)*a(2,2)
b(2,3) = a(0,3)*a(1,1)*a(2,0)+a(0,0)*a(1,3)*a(2,1)+a(0,1)*a(1,0)*a(2,3) &
        -a(0,3)*a(1,0)*a(2,1)-a(0,0)*a(1,1)*a(2,3)-a(0,1)*a(1,3)*a(2,0)
b(3,3) = a(0,0)*a(1,1)*a(2,2)+a(0,1)*a(1,2)*a(2,0)+a(0,2)*a(1,0)*a(2,1) &
        -a(0,0)*a(1,2)*a(2,1)-a(0,1)*a(1,0)*a(2,2)-a(0,2)*a(1,1)*a(2,0) 
!Dividing by determinant
det = a(0,0)*b(0,0)+a(0,1)*b(1,0)+a(0,2)*b(2,0)+a(0,3)*b(3,0)
IF(det.eq.0.d0) THEN
 print*, 'WARNING:NON-INVERTIBLE MATRIX'
 STOP
ENDIF 
absmax = absmax/det
a = b*absmax
do i = 0,3
 do j=0,3
 if(DABS(a(i,j)).le.snap) a(i,j) = 0.d0
 enddo!j
enddo !i
!------------------------------------------------------------------------
return
!------------------------------------------------------------------------
END SUBROUTINE MatInv
!=======================================================================
END MODULE math
