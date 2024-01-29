MODULE parameters
!------------------------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------------
DOUBLE PRECISION, PARAMETER :: pc=2.99792458d10
DOUBLE PRECISION, PARAMETER :: ph=6.626068d-27
DOUBLE PRECISION, PARAMETER :: pk=1.3806503d-16
DOUBLE PRECISION, PARAMETER :: pi=3.1415926535d0
DOUBLE PRECISION, PARAMETER :: me=9.1093826d-28
DOUBLE PRECISION, PARAMETER :: mh=1.6605d-24
DOUBLE PRECISION, PARAMETER :: sthom=6.653d-25
!------------------------------------------------------------------------------
!------------------------------  FLAGS  ---------------------------------------
! fl_red: Redistribution. 0: CRD, 1: Coherent scattering, 2: PRD
INTEGER, PARAMETER :: fl_red = 2
! fl_anom: Elements of the redistribution matrix. 0: Accounting for all elements
! 1: Setting anomalous dispersion terms to zero.
INTEGER, PARAMETER :: fl_anom = 0 
!
INTEGER, PARAMETER :: fl_calc = 0
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: niter= 3300  !max number of iterations
!--LINE
DOUBLE PRECISION, PARAMETER :: vstepL = 0.02d0   !fine step int. nu' [Dopp. width units;DWU]
DOUBLE PRECISION, PARAMETER :: vstepCrs = 0.5d0!0.5d0  !Coarse step int. nu' [DWU]
DOUBLE PRECISION, PARAMETER :: deltaL=6.d0  !semi-interv. int. nu' [DWU]
DOUBLE PRECISION, PARAMETER :: deltaV=4.d0  !critical separat. from line [DWU]
DOUBLE PRECISION, PARAMETER :: deltaVR=3.5d0!critical separat. from line [DWU, with Raman]
DOUBLE PRECISION, PARAMETER :: deltaT=5.d0
DOUBLE PRECISION, PARAMETER :: thdG2=1.d-20 !threshold value of weights G2 
                                  !(weights set to 0 if below)
!--CONTINUUM
DOUBLE PRECISION, PARAMETER :: vstepC=0.2   !fine step int. nu' [Dopp. width units]
DOUBLE PRECISION, PARAMETER :: deltaC=4.d0  !semi-interv. int. nu' [Dopp. width units]
DOUBLE PRECISION, PARAMETER :: thdC=1.d-20  !threshold value of weights Cij
                                  !(weights set to zero if below)
!--INTEGRAL OVER DIRECTIONS
INTEGER, PARAMETER :: ndir = 18!18 !20     !Quadr. order Gauss-Leg. int. over mu

INTEGER, PARAMETER :: ndirx = 8 !32 Trapezoidal integration for azimuth
!--INTEGRAL FOR ANGLE-AVERAGED R_II 
!(these parameters are not used in case of Guttebroze approx.)
INTEGER, PARAMETER :: nqgl = 15!60    !Quadr. order Gauss-Leg. int. over theta 
                                  !(non resonance frequency)
INTEGER, PARAMETER :: nqglc = 100!150   !Quadr. order Gauss-Leg. int. over theta 
                                  !(resonance frequency)
INTEGER, PARAMETER :: nqglf= 45!60
!------------------------------------------------------------------------------
DOUBLE PRECISION :: mu(ndir), wtdir(ndir)   !abscis. and weig. quadr. Gauss-Leg. mu
DOUBLE PRECISION :: az(ndirx), waz(ndirx)      !Angle and weig. for azimuth using trap. int.
DOUBLE PRECISION :: ygl(nqgl),wgl(nqgl)     !abscis. and weig. quadr. Gauss-Leg. theta
DOUBLE PRECISION :: yglc(nqglc),wglc(nqglc) !abscis. and weig. quadr. Gauss-Leg. theta
DOUBLE PRECISION :: yglf(nqglf),wglf(nqglf)
!------------------------------------------------------------------------------
END MODULE parameters
