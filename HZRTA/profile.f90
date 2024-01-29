MODULE profile
!-----------------------------------------------------------------------
USE parameters
!-----------------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------


!===========================================================================
SUBROUTINE prof(gam,arg,phi,psi)
!-----------------------------------------------------------------------
! Subroutine for calculating the Lorentzian profile (phi) and the dispersion 
! profile, using as inputs the Gamma value (gam) and the argument arg=nu0-nu
!-----------------------------------------------------------------------
  DOUBLE PRECISION, INTENT(IN) :: gam, arg
  DOUBLE PRECISION, INTENT(OUT) :: phi, psi
  DOUBLE PRECISION :: gam2, arg2
!-----------------------------------------------------------------------
  gam2=gam*gam
  arg2=arg*arg
  phi=gam/(pi*(gam2+arg2))
  psi=arg/(pi*(gam2+arg2))
!-----------------------------------------------------------------------
END SUBROUTINE prof
!========================================================================

!========================================================================
SUBROUTINE PROFILARCH(A,V,VO,GA)
!------------------------------------------------------------------------
! FROM REICHEL,JQSRT 8,1601 (FOR A.LT.1E-8, ASSUME A=1E-8);
! FOR V.GT.10 OR A.GT.3, TAKE  ABRAMOWITZ&STEGUN,P.328 
!------------------------------------------------------------------------
 INTEGER :: N
 DOUBLE PRECISION, INTENT(OUT) :: VO, GA
 DOUBLE PRECISION  :: A, V
 DOUBLE PRECISION :: A1, A2, EX, U, U0, U1, S, S0, S1, T, Q, V2
!------------------------------------------------------------------------
 IF(A.LT.1.E-8) A=1.E-8
 IF(A.GT.3..OR.ABS(V).GT.10.) GOTO 210
 V2=V*V
 A2=A*A*2.
 U0=A*1.7724538509
 U1=-A2
 S0=U0+U1
 N=2
 206 U0=U0*A2/N
 N=N+1
 U1=U1*A2/N
 N=N+1
 U=U0+U1
 S0=S0+U
 IF(ABS(S0).LT.1E-20) GO TO 206
 IF(ABS(U/S0).GT.1E-8) GO TO 206
 N=1
 T=1.
 S=S0
 S1=0.
  207 S0=(1-S0)*A2/(2*N-1)
  Q=T*S0
  S1=Q+S1
  T=T*V2/N
  U=T*S0
  S=S+U
  N=N+1
  IF(ABS(S1).LT.1E-20) GO TO 207
  IF(ABS(Q/S1).GT.1E-8) GO TO 207
  IF(ABS(S).LT.1E-20) GO TO 207
  IF(ABS(U/S).GT.1E-8) GO TO 207
  EX=EXP(-V2)
  S=S*EX
  S1=S1*EX*V/A
  A1=A*1.7724538509
  VO=S/A1
  GA=S1/A1
  RETURN
  210 CALL PROFAS(A,V,VO,GA)
  RETURN
!------------------------------------------------------------------------
END SUBROUTINE PROFILARCH
!========================================================================

!========================================================================
SUBROUTINE PROFAS(A,V,VO,GA)
!------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: A, V 
DOUBLE PRECISION, INTENT(OUT) :: VO, GA
DOUBLE PRECISION :: ZERO, UNO, A1, A2, A3, B1, B2, B3
DOUBLE COMPLEX :: Z,Z2,A0C,A1C,A2C,A3C,B1C,B2C,B3C,RES
DATA A1/0.4613135/
DATA A2/0.09999216/
DATA A3/0.002883894/
DATA B1/0.1901635/
DATA B2/1.7844927/
DATA B3/5.5253437/
!------------------------------------------------------------------------
ZERO=0.
UNO=1.
Z=DCMPLX(V,A)
Z2=Z*Z
A0C=DCMPLX(ZERO,UNO)
A1C=DCMPLX(A1,ZERO)
A2C=DCMPLX(A2,ZERO)
A3C=DCMPLX(A3,ZERO)
B1C=DCMPLX(B1,ZERO)
B2C=DCMPLX(B2,ZERO)
B3C=DCMPLX(B3,ZERO)
RES=A0C*Z*(A1C/(Z2-B1C)+A2C/(Z2-B2C)+A3C/(Z2-B3C))
VO=DREAL(RES)
GA=DIMAG(RES)
!------------------------------------------------------------------------
RETURN
!------------------------------------------------------------------------
END SUBROUTINE PROFAS
!========================================================================

END MODULE profile
