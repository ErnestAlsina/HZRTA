PROGRAM hz_rta
!----------------------------------------------------------------------
! The program solves the non-LTE problem of the second kind for a 
! two-level atom with unpolarized and infinitely sharp lower level, in 
! the presence of deterministic weak magnetic fields, but without using
! the weak-field approximation.
! It considers redistribution effects due to:
! 1) elastic collisions. 
!    The redistribution matrix of Bommier (1997) is considered.
!     As far as R_III is concerned, the 
!    assumption of CRD in the observer's frame is made.
! The non-LTE problem is solved through jacobian iterative method.
!-----------------------------------------------------------------------
USE global_vars
USE parameters
USE math
USE profile
USE abs_prof
USE redis
USE matrmod
USE rt_formal
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER :: iz, i, ip, j, k, iter, ivac, in, ku,ik, iq,in1,ll2,lu2,sl2,su2
INTEGER :: im, ax,qmax,mu3,mup3,ml3,mlp3,p1,p2,p3,p4
INTEGER :: fl1, fl_lred, fl_cred,fl_lam,miter,kp,kpp,tg,tgm,si
DOUBLE PRECISION :: wl0, wl0vac, wlvac, v0, norm
DOUBLE PRECISION :: etaI_l, wte, wth!,mi,chi
DOUBLE PRECISION :: x0, z0
DOUBLE PRECISION :: polb,r2v,tmf
DOUBLE PRECISION :: fac, t1,t2,time
DOUBLE PRECISION :: emin, emax,val1,val2,val3,val4
DOUBLE PRECISION, DIMENSION(:), allocatable :: mi, chi
DOUBLE PRECISION, DIMENSION(:), allocatable :: wl, Bpb, Wm
DOUBLE PRECISION, DIMENSION(:), allocatable :: se, sray, dnde, dndh
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: etaI_c, chi_s,&
                                         eps_c, k_c, s, sh, ss
DOUBLE PRECISION, DIMENSION(:,:,:,:), allocatable :: etaI, r
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: setaI, sr
DOUBLE PRECISION, DIMENSION(:), allocatable :: T, popl, dld, k_l, ne, nh, vmic
DOUBLE PRECISION, DIMENSION(:), allocatable :: Cul, Qel
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: D2
DOUBLE PRECISION, DIMENSION(:), allocatable :: epsln!, epsprm
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: Bp
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: emer,emer_lin
!DOUBLE PRECISION, DIMENSION(:,:), allocatable :: O20, O21r, O21i, O22r, O22i
INTEGER, DIMENSION(:), allocatable :: indx,indxb
INTEGER, DIMENSION(:,:), allocatable :: indxsto!,indxtob
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: m
DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), allocatable :: lstar
DOUBLE PRECISION, DIMENSION(:,:,:), allocatable :: msto, CT
DOUBLE PRECISION, DIMENSION(:), allocatable :: rho
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: rhol,add
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: jint, sline, scont, dline
DOUBLE PRECISION :: jt,phi,psi,arg!Frequency-integrated radiation field tensor
!
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: keepline,keepcont
DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), allocatable :: mmat, invk!, kmat
DOUBLE PRECISION, DIMENSION(:,:,:,:), allocatable :: smmat, sinvk,em_coef
DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:,:), allocatable :: r2,r3
DOUBLE PRECISION :: d, DS0, RDS0, MRDS0, ht_MRDS0, wl_MRDS0
DOUBLE PRECISION :: DS2, RDS2, MRDS2, ht_MRDS2, wl_MRDS2, MRDS0prev
DOUBLE PRECISION :: DS1, RDS1,MRDS1
DOUBLE PRECISION :: DSR21, RDSR21, DSI21, RDSI21, DSR22, RDSR22, DSI22, RDSI22
DOUBLE PRECISION :: MRDSR11,MRDSI11,RDSR11,RDSI11,DSR11,DSI11
DOUBLE PRECISION :: MRDSR21, MRDSI21
DOUBLE PRECISION :: MRDSR22,MRDSI22
DOUBLE PRECISION :: Hcf,eneru,enerl,eneriv,eneruv,ae
DOUBLE PRECISION :: alpha,beta, gamma
DOUBLE PRECISION :: timi,timf,sour_re,sour_im,em_a,ctrl,tre,tim
DOUBLE PRECISION :: wthd
DOUBLE PRECISION, DIMENSION(0:2,-2:2,0:3) :: TR, TI
DOUBLE PRECISION, DIMENSION(:,:,:), allocatable :: eta_ref,rho_ref
CHARACTER(len=1024) :: fileis,fileib,fileic,fillin
CHARACTER(len=1024) :: format_string,format_strb,format_strc
CHARACTER(len=1024) :: filtst,format_tst,filskh,format_skh
!-------
DOUBLE PRECISION :: mu1(50), mu2(50), wtdir1(50), wtdir2(50), theta
INTEGER :: ndirh
!-------
! ndir, mu(ndir), wtdir(ndir) :: parameters. Number of directions, mu
! and weights
! nqgl, ygl(nqgl), wgl(nqgl) :: parameters. Quadrature for R_II-AA. Not
! used for CRD or for Gouttebroze approx
! nqglc, yglc(nqgl), wglc(nqgl) :: parameters. Idem
! niter :: parameters. Number of iterations
! nwl, nu(nwl), Aul :: global_vars. Total number of freq, freq, and line
! width
! nz, ht(nz), dnd(nz), ad(nz) :: global_vars. 
! prof_vo(in,nz,nwl), G3(2,2,5,nz,nwl,0:2) :: global_vars. To be allocated
!=======================================================================
!============================= INITIALIZE ==============================
!=======================================================================
 call cpu_time(timi)
 call factrl
!-----------------------------------------------------------------------
! Calculate weights and abscissas for gaussian quadrature over 
! directions (needed for calculating J^K_Q).
! Once the total number of directions (ndir, parameter) is fixed, 
! calculate separately weights and abscissas in the two quadrants 
! (mu=[-1,0], mu=[0,1]) for a number of directions equal to half the 
! total directions. Consider first incoming and then outgoing directions
!-----------------------------------------------------------------------
ndirh=ndir/2
call LGDR(ndirh,-1.d0,0.d0,mu1,wtdir1)
call LGDR(ndirh,0.d0,1.d0,mu2,wtdir2)
do i=1,ndir
  if(i.le.ndirh) then
    mu(i)=mu1(i)
    wtdir(i)=wtdir1(i)
  else
    mu(i)=mu2(i-ndirh)
    wtdir(i)=wtdir2(i-ndirh)
  endif
enddo
!--------------------------------------------------------------------
! For the azimuthal angle we use the trapezoidal method of integration.
!--------------------------------------------------------------------
 do i = 1, ndirx
  az(i) = 2.d0*pi*DBLE(i-1)/DBLE(ndirx)
  waz(i) = 2.d0*pi/DBLE(ndirx) 
 enddo !i (azimuthal angle)
!---------------------------------------------------------------------
! Angular quadratures for the RII weights
!---------------------------------------------------------------------
 call LGDR(nqgl,0.d0,1.d0,ygl,wgl)
 call LGDR(nqglc,0.d0,1.d0,yglc,wglc)
 call LGDR(nqglf,0.d0,1.d0,yglf,wglf)
!---------------------------------------------------------------------
! Write the values of the T^K_Q tensors for all possible angles (to be
! used as a global variable)
!---------------------------------------------------------------------
!print*, 'Find the T^K_Q tensors for the whole angular grid'
!call cpu_time(t1)
allocate(TTR(0:2,-2:2,0:3,1:ndir,1:ndirx),TTI(0:2,-2:2,0:3,1:ndir,1:ndirx))
gamma = pi/2
do i = 1, ndir
 beta = DACOS(mu(i))
 do j = 1, ndirx
  alpha = az(j)
  call TKQ2(-gamma,-beta,-alpha,0.d0,0.d0,0.d0,TR,TI)
  TTR(:,:,:,i,j) = TR(:,:,:)
  TTI(:,:,:,i,j) = TI(:,:,:)
 enddo!j
enddo!i
!call cpu_time(t2)
!print*, 'Tensors found. Time (in seconds):', t2-t1
!---------------------------------------------------------------------
! Define the global matrix iden (the identity matrix)
!---------------------------------------------------------------------
do i = 0, 3
 do j = 0, 3
  iden(i,j) = 0.d0
  if(i.eq.j) iden(i,j) = 1.d0
 enddo
enddo
!=======================================================================
!========================== READ INPUT FILES ===========================
!=======================================================================
! Read input file 'lambda.dat'
! -Read parameter 'ivac'
!     ivac=0 -> air wavelengths in input
!     ivac=1 -> vacuum wavelengths in input
! -Read line center wavelength 'wl0' [A]
! -Read number of wavelengths 'nwl'
! -Read values of wavelengths 'wl' [A]
! Wavelength grid calculated by Han's code
!-----------------------------------------------------------------------
open(unit=1,file='lambda.dat')
read(1,*) ivac
read(1,*) wl0
read(1,*) nwl
allocate(wl(1:nwl))
do i=1,nwl
 read(1,*) wl(i)
enddo
close(unit=1)
!-----------------------------------------------------------------------
! Read input file 'radfield.dat'
! -Read values of 'J00' and 'J20' at each height and frequency.
!  [erg cm^-2 s^-1 Hz^-1]
! These values correspond to the converged solution of Han's code, and 
! are used as initial guess.
!----------------------------------------------------------------------
open(unit=2,file='radfield.dat')
read(2,*) nwl,nz
! For radout, read wavl, comment for readfield
allocate(jrad(1:2,0:2,0:2,1:nz,1:nwl))
jrad = 0.d0! Initialize as zero
do iz = 1,nz
 read(2,*) (jrad(1,0,0,iz,i),i=1,nwl)
enddo
do iz=1,nz
 read(2,*) (jrad(1,2,0,iz,i),i=1,nwl)
enddo
close(unit=2)
!-----------------------------------------------------------------------
! Read input file 'continuum.dat'
! -Read values of continuum opacity and emissivity at each height and 
!  frequency:
!  total opacity (absorpt. + scatt.) 'etaI_c' [cm^-1]
!  scattering 'chi_s' [cm^-1]
!  emissivity (without scatt.) 'eps_c' [erg cm^-3 s^-1 Hz^-1 sr^-1]
! These quantities are calculated by Han's code.
!-----------------------------------------------------------------------
allocate(etaI_c(1:nz,1:nwl),chi_s(1:nz,1:nwl),eps_c(1:nz,1:nwl))
open(unit=3,file='continuum.dat')
read(3,*) nwl,nz
do iz=1,nz
  read(3,*) (etaI_c(iz,i),i=1,nwl)
  read(3,*) (chi_s(iz,i),i=1,nwl)
  read(3,*) (eps_c(iz,i),i=1,nwl)
enddo
close(unit=3)
!----------------------------------------------------------------------
open(unit=13,file='isot_2L.dat')
read(13,*) ae
read(13,*) jl2,ll2,sl2
read(13,*) enerl
read(13,*) ju2,lu2,su2
read(13,*) eneru
read(13,*) atwe
read(13,*) eneriv,eneruv
close(unit=13)
v0=eneru-enerl
!-----------------------------------------------------------------------
! Read input file 'atmosphere.dat'
! -Read value of Einstein coeff. for spontaneous emission 'Aul' [s^-1]
! -At each height read:
!     height 'ht' [cm]
!     temperature 'T' [K]
!     lower level population 'popl' [cm^-3]
!     doppler broadening 'dld' [mA]
!     damping constant 'ad'
!     rate superelastic collisions 'Cul' [s^-1]
!     rate elastic collisions 'Qel' [s^-1]
!     rate depolarizing collisions 'D2' [s^-1]
!     electron number density 'ne' [cm^-3]
!     neutral hydrogen number density 'nh' [cm^-3]
!     microturbulent velocity 'vmic' [cm s^-1]
!-----------------------------------------------------------------------
allocate(ht(1:nz),T(1:nz),popl(1:nz),dld(1:nz),ad(1:nz))
allocate(Cul(1:nz),Qel(1:nz),D2(0:ju2,1:nz),ne(1:nz),nh(1:nz),vmic(1:nz))
open(unit=4,file='atmosphere.dat')
read(4,*) Aul
read(4,*) nz
do iz=1,nz
 read(4,*) ht(iz), T(iz), popl(iz), dld(iz), ad(iz), Cul(iz),&
 Qel(iz), D2(2,iz), ne(iz), nh(iz), vmic(iz)
enddo
close(unit=4)
open(unit=11,file='height.res') 
 write(11,*) nz
 write(11,*) (ht(iz),iz = 1,nz)
close(unit=11)
!-----------------------------------------------------------------------
!Modifying the branching ratios for CRD and coherent scattering limits
!-----------------------------------------------------------------------
r2v = ((5.d0/2.d0)*(13.6d0/(eneriv-eneruv))**2)**0.4d0
do iz =1,nz
 tmf = (T(iz)*(1.d0 +(1.d0/atwe)))**0.3d0
 D2(0,iz) = 0.d0!Should always be zero
do ik = 1,ju2
  D2(ik,iz) = 0.d0!0.5d0*Qel(iz)
 enddo
enddo
!-----------------------------------------------------------------------
select case(fl_red)
 case(0)! CRD
  do iz = 1,nz 
   Qel(iz) = 1.d16!
  enddo
 case(1)! Coherent scattering
  do iz = 1,nz
   Qel(iz) = 0.d0
  enddo
 case default
end select
!-----------------------------------------------------------------------
! Read input file 'param.dat'
! Read the following parameters of the problem:
! -mi: cosine heliocentric angle for which the Stokes parameters of the 
!      emergent radiation are calculated
! -fl_lred: flag on the calculation type (RII, RIII, RII+RIII)
! -fl_cred: flag on continuum calculation (CS unpol, CS pol, non-CS pol)
!-----------------------------------------------------------------------
open(unit=10,file='param.dat')
 read(10,*) tgm
allocate(mi(1:tgm),chi(1:tgm))
do i=1,tgm
 read(10,*) mi(i)
 read(10,*) chi(i)
enddo
read(10,*) fl_lred
read(10,*) fl_cred
close(unit=10)
!-----------------------------------------------------------------------
! Read input file magnetic.dat
! Reads the following parameters
! -Bint: Magnetic field strength, in G
! -them: Polar angle of the direction of the magnetic field with respect 
!        to the vertical, in radians. Used for the considered emergent
!        intensity, Stokes Q, and Stokes U
! -chim: Azimuthal angle of the direction of the magnetic field, in radians.
!        Used for the considered emergent intensity, Stokes Q and Stokes U
! -gu: Lande factor of the upper level
! -gl: Lande factor for the lower level
!----------------------------------------------------------------------
open(unit=9,file ='magnetic.dat')
read(9,*) Bint
read(9,*) them
read(9,*) chim
read(9,*) fl_lam
close(unit=9)
gu = 1.d0 + &
 (ju2*(ju2+2.d0)+su2*(su2+2.d0)-lu2*(lu2+2.d0))/(2.d0*ju2*(ju2+2.d0))
if(jl2.ne.0) then
 gl = 1.d0 + &
 (jl2*(jl2+2.d0)+sl2*(sl2+2.d0)-ll2*(ll2+2.d0))/(2.d0*jl2*(jl2+2.d0))
else
 gl = 1.d0 
endif
Hu = 8.79d6*gu*Bint/Aul
print*, 'Hanle critical field', Aul/(8.79d6*gu)
nular = 1.39962334235d6*Bint
!==============================================================================
! =========   CALCULATE QUANTITIES THAT WILL BE USED LATER ====================
!==============================================================================
! Transforms wavelengths into frequencies 'vp' [Hz] (decreasing order)
!------------------------------------------------------------------------------
allocate(nu(1:nwl))
if(ivac.eq.1) then
 wl0vac=wl0
else
 call wav_air_vac(wl0,wl0vac)
endif
do i=1,nwl
 if(ivac.eq.1) then
  wlvac=wl(i)
 else
  call wav_air_vac(wl(i),wlvac)
 endif
  nu(i)=pc/(wlvac*1.d-8)
enddo
!-------------------------------------------------------------------------------
! Calculate factor fr3=(nu/nu0)^3 (global variable)
! Do not use! May create problems! Origin of problems should be invest.
!--------------------------------------------------------------------------------
allocate(fr3(1:nwl))
do i=1,nwl
!    fr3(i)=(nu(i)/v0)**3
 fr3(i)=1.d0
enddo
!------------------------------------------------------------------------
! At each height calculate:
!   -Doppler broadening in frequency 'dnd' [Hz]
!   -Wien function at frequency v0 'Wm' [erg cm^-2 s^-1 Hz^-1 sr^-1]
!   -Planck function 'Bp' [erg cm^-2 s^-1 Hz^-1 sr^-1]
!   -epsilon 'epsln'
!   -epsilon prime 'epsprm'
!   -Planck function at deepest height 'Bpb'
!-----------------------------------------------------------------------
allocate(dnd(1:nz),Bp(1:nz,1:nwl),Bpb(1:nwl),Wm(1:nz),epsln(1:nz),epsprm(1:nz))
allocate(Hup(0:ju2,1:nz))
allocate(Hup2(1:nz))!Only needed for K=2
do iz=1,nz
 dnd(iz)=(dld(iz)*1.d-3)*(v0**2)/(pc*1.d8)
 Wm(iz)=wien(v0,T(iz))
 do i=1,nwl
  Bp(iz,i)=planck(nu(i),T(iz))
 enddo
 epsln(iz)=Cul(iz)/(Aul+Cul(iz))
 epsprm(iz)=Cul(iz)/Aul
enddo
!
do i=1,nwl
 Bpb(i)=Bp(nz,i) !Planck function for bottom boundary
enddo
!-----------------------------------------------------------------------
! At each height calculate:
!   -line absorption profile at each frequency 'prof_vo' [s]
!   -frequency integrated line absorption profile 'k_l' [cm^-1 s^-1]
!-----------------------------------------------------------------------
allocate(k_l(1:nz))
print*,' -- Calculating absorption profile...'
call prof_2L(z0)
! Subroutine first calculates z0 and then the Voigt profile is 
! returned and normalized
do iz=1,nz
 k_l(iz)=z0*popl(iz) !Gives k_L
enddo
print*, 'k_L found. Calculating generalized profile'
allocate(gprof_re(0:ju2,0:2,-2:2,1:nz,1:nwl),gprof_im(0:ju2,0:2,-2:2,1:nz,1:nwl))
call gener(gprof_re,gprof_im)! Calculates generalized profile (real & imag)
print*, 'Generalized profile found'
!--------------------------------------------------------------------------
! Also find the C-factors
!--------------------------------------------------------------------------
print*,'Cfactors',ju2,jl2
allocate(cfac(0:2,0:2,-2:2,-ju2:ju2,-ju2:ju2,-jl2:jl2,-jl2:jl2))
do kp = 0,2
 do kpp = 0,2
  qmax = min(kp,kpp)
  do iq=-qmax,qmax
   do mu3 = -ju2,ju2,2
    do mup3 = -ju2,ju2,2
     do ml3 = -jl2,jl2,2
      do mlp3 = -jl2,jl2,2
       cfac(kp,kpp,iq,mu3,mup3,ml3,mlp3) = 0.d0
       do p1 =-2,2,2
        if(p1.ne.mu3-ml3) CYCLE
        do p2 = -2,2,2
         if(p2.ne.mup3-ml3) CYCLE
         do p3 =-2,2,2
          if(p3.ne.mu3-mlp3) CYCLE
          do p4 = -2,2,2
           if(p4.ne.mup3-mlp3) CYCLE
           in1 = 2*ju2-ml3-mlp3
           if(MOD(in1,2).ne.0) then
            print*,'Error in indices'
            stop
           endif
           fac = 1.d0
           if(MOD(in1,4).ne.0) fac = -1.d0
           val1 = 3.d0*dble(ju2+1)*sqrt(dble(2*kp+1)*dble(2*kpp+1))
           val2 = W3JS(ju2,jl2,2,mu3,-ml3,-p1)*W3JS(ju2,jl2,2,mup3,-ml3,-p2)
           val3 = W3JS(ju2,jl2,2,mu3,-mlp3,-p3)*W3JS(ju2,jl2,2,mup3,-mlp3,-p4)
           val4 = W3JS(2,2,2*kp,-p1,p2,2*iq)*W3JS(2,2,2*kpp,-p3,p4,2*iq)
           cfac(kp,kpp,iq,mu3,mup3,ml3,mlp3)=cfac(kp,kpp,iq,mu3,mup3,ml3,mlp3)+&
           fac*val1*val2*val3*val4
          enddo!p4
         enddo!p3
        enddo!p2
       enddo!p1
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
!-----------------------------------------------------------------------
! Calculate the RT coefficient for the propagation matrix
! Requires calling the rotation matrices
!-----------------------------------------------------------------------
call ROTMATSU(2,0.d0,-them,-chim,DR,DI)
print*, 'Rotation matrix found. Calculating radiative transfer coefficients'
allocate(eta_re(0:3,1:nz,1:nwl,1:ndir,1:ndirx), eta_im(0:3,1:nz,1:nwl,1:ndir,1:ndirx))
allocate(rho_re(0:3,1:nz,1:nwl,1:ndir,1:ndirx), rho_im(0:3,1:nz,1:nwl,1:ndir,1:ndirx))
call etarho(k_l,eta_re,eta_im,rho_re,rho_im) ! eta and rho RT coefficients found
print*, 'Radiative transfer coefficients found (depend only on atmospheric parameters and magnetic field)'
if(fl_anom.ne.0) then
 do iz = 1,nz
  do i = 1,nwl
   do im = 1,ndir
    do ax = 1,ndirx
     do si = 1,3
      rho_re(si,iz,i,im,ax) = 0.d0
      rho_im(si,iz,i,im,ax) = 0.d0
! Add condition for neglecting also dichroism terms
     enddo 
    enddo
   enddo
  enddo
 enddo
endif
! --------------------------------------------------------------------------------------
! Calculate for each each height, frequency, and angle:
!  - eta_I total: adding the previously calculated line eta and the continuum eta
!  - r factor
! --------------------------------------------------------------------------------------
allocate(etaI(1:nz,1:nwl,1:ndir,1:ndirx), r(1:nz,1:nwl,1:ndir,1:ndirx))
do iz=1,nz
 do i=1,nwl
  do im=1,ndir
   do ax = 1, ndirx
    etaI(iz,i,im,ax) = eta_re(0,iz,i,im,ax) + etaI_c(iz,i)!continuum part doesnt depend on angle
    r(iz,i,im,ax) = eta_re(0,iz,i,im,ax)/etaI(iz,i,im,ax)
   enddo !ax
  enddo !im
 enddo !i
enddo !iz
print*,'Eta'
! -----------------------------------------------------------------------------------
! Calculate the angle-dependent part of the redistribution matrix (only R_III)
! since we are considering CRD.
! K and K' go from 0 to 2. |Q| can not exceed either K or K'
! -----------------------------------------------------------------------------------
allocate(r3(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz))
allocate(r2(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz))
!if(fl_calc.eq.2) then
!  call stored_redis(r2,r3)
!else
do iz = 1, nz
 do ik = 0,ju2 
  Hup(ik,iz) = Hu/(1.d0+D2(ik,iz)/Aul+epsprm(iz))
 enddo
 Hup2(iz) = Hu/(1.d0+Qel(iz)/Aul+epsprm(iz))
enddo !iz
print*,'Before R3'
call crd_red(epsprm,D2,Qel,Hup,Hup2,gprof_re,gprof_im,r3)
print*, 'R3 found'
do kp = 0,2
 do kpp = 0,2
  qmax = min(kp,kpp)
  do iz = 1,nz
   do i = 1,nwl
    do j = 1,nwl
     do iq=-qmax,qmax
      r2(1,kp,kpp,iq,i,j,iz) = 0.d0
      r2(2,kp,kpp,iq,i,j,iz) = 0.d0
     enddo
    enddo
   enddo 
  enddo
 enddo
enddo
!
call coh_red(epsprm,Qel,Hup2,gprof_re,gprof_im,r2)
!-------------------------------------------------------------------------------------- 
! Continuum part. Since we are only applying the frequency-by-frequency method
! to the S^0_0 of the line part, we will not be needing Cij for the M_ij matrix.
! This might change if we are not using coherent scattering
!-------------------------------------------------------------------------------------- 
allocate(s(1:nz,1:nwl))
do iz = 1, nz
 do i = 1, nwl
  s(iz,i) = (etaI_c(iz,i)-chi_s(iz,i))/etaI_c(iz,i)
 enddo !i
enddo !iz
!----------------------------------------------------------------------------
! Set dimensions of variables that will be used in the following steps
!----------------------------------------------------------------------------
allocate(jint(1:2,0:2,0:2,1:nz,1:nwl))
allocate(sline(1:2,0:2,0:2,1:nz,1:nwl),scont(1:2,0:2,0:2,1:nz,1:nwl))
allocate(keepline(1:2,0:2,0:2,1:nz,1:nwl),keepcont(1:2,0:2,0:2,1:nz,1:nwl))
allocate(dline(1:2,0:2,0:2,1:nz,1:nwl))
allocate(indx(1:nwl),indxsto(1:nz,1:nwl))
allocate(m(1:nwl,1:nwl),msto(1:nz,1:nwl,1:nwl))
allocate(lstar(0:3,0:3,1:nz,1:nwl,1:ndir,1:ndirx))
allocate(rho(1:nwl),rhol(1:nz,1:nwl))
!----------------------------------------------------------------------------
! Start the iterative scheme
!----------------------------------------------------------------------------
open(unit=11,file='convergence.res')
open(unit=41,file='MRC_points.res')
fl1=0; df=0 
keepline = 0.d0
keepcont = 0.d0
!----------------------------------------------------------------------------
! 0) Initialize source functions other than S^0_0. Used in order to find 
!     variation in the second iteration
!----------------------------------------------------------------------------
do iz = 1, nz 
 do i = 1, nwl
  do ik =0, 2
   do iq = 0,ik
    do in=1,2
     sline(in,ik,iq,iz,i) = 0.d0
     scont(in,ik,iq,iz,i) = 0.d0
    enddo!in
   enddo !iq
  enddo !ik 
 enddo !i
enddo !iz
MRDS0prev = 0.d0!In order to declare it. Not actually important until iteration 100
allocate(add(1:nz,1:nwl))
open(unit=67,file="emis_comp_BESSER.res")
write(67,*) "First iteration"
do iter = 1, niter
 write(*,'(A,I6)'), 'Iteration:', iter
!-----------------------------------------------------------------------------
! 1) Calculate the frequency-integrated radiation tensor. We need the 
!     radiation tensors, and the frequency-dependent part of the redistribution
!     matrix, as well as the rotation matrices (present in global_vars)
!-----------------------------------------------------------------------------
 call red_rad(k_l,jrad,r2,r3,jint)
!-----------------------------------------------------------------------------
! 2) Find the multipolar components of the source functions: line and 
!    continuum parts. From the second iteration onwards we shall use the 
!    frequency-by-frequency method in order to find the S^0_0 component of the
!    line part
!-----------------------------------------------------------------------------
 IF(iter.eq.1) THEN
  iz=56
  write(67,*) "Line scattering E00"
  write(67,*) (dnd(iz)*jint(1,0,0,iz,i),i=1,nwl)
  fac = k_l(iz)*Wm(iz)*epsprm(iz)/(epsprm(iz) + 1.d0)
  write(67,*) "Line thermal E00"
  write(67,*) (dnd(iz)*fac*gprof_re(0,0,0,iz,i),i=1,nwl)
  call lambd_line(k_l,epsprm,Wm,gprof_re,gprof_im,jint,sline)
  call lambd_cont(s,eps_c,etaI_c,jrad,scont)!
  iz = 56
  write(67,*) "Continuum E00"
  write(67,*) (dnd(iz)*scont(1,0,0,iz,i),i=1,nwl)
!------------------------------------------------------------------------------
  do iz = 1, nz
   fac = k_l(iz)*Wm(iz)*epsprm(iz)/(epsprm(iz) + 1.d0)
   do i = 1, nwl
    sline(1,0,0,iz,i) = fac*gprof_re(0,0,0,iz,i) + jint(1,0,0,iz,i)!In following iterations
   enddo!i
  enddo!iz
!
  allocate(invk(0:3,0:3,1:nz,1:nwl,1:ndir,1:ndirx))!
  allocate(mmat(0:3,0:3,1:nz,1:nwl,1:ndir,1:ndirx))! 1*e^-Dt - K' (K/etaI)
  call lam_star(etaI,eta_re,rho_re,invk,mmat,lstar) !
 ELSE 
!----------------------------------------------------------------------------
! Following iterations: using Jacobi method
!----------------------------------------------------------------------------
  IF(iter.eq.2) THEN
!----------------------------------------------------------------------------
! Only for iteration 2: find the m-matrix for the FBF method
!----------------------------------------------------------------------------
   print*, '----- Obtaining the M-matrix and LU decomposition -----'
   do iz=1,nz
    if(fl_lam.ne.2) then
     call m_mat(iz,k_l,r,r2,r3,lstar,eta_re,m)!
    else
     do i=1,nwl
      do j=1,nwl
       m(i,j) = 0.d0
       if(i.eq.j) m(i,j) = 1.d0
      enddo
     enddo
    endif
!----------------------------------------------------------------------------
! Apply LUD method and save new (decomposed) matrices
! --------------------------------------------------------------------------
    call ludcmp(m,nwl,indx,d)
    do i = 1,nwl
     indxsto(iz,i) = indx(i)
     do ip = 1,nwl
      msto(iz,i,ip) = m(i,ip)
     enddo!ip
    enddo!i
   enddo !iz
  ENDIF!ITERATION 2
  MRDS0 = 0.d0; MRDS2 =0.d0!
  ht_MRDS0 = 0.d0; ht_MRDS2 = 0.d0
  wl_MRDS0 = 0.d0; wl_MRDS2 = 0.d0
  MRDS1 = 0.d0; MRDSR11 = 0.d0; MRDSI11 = 0.d0
  MRDSR21 = 0.d0; MRDSI21 = 0.d0
  MRDSR22 = 0.d0; MRDSI22 = 0.d0
!---------------------------------------------------------------------------
! Use lambda-iteration for all continuum source function multipolar components
! and for all line components except S^0_0. For S_0^0 we use the frequency-
! -by-frequency method. The order of one or the other doesn't matter. We'll
! do lambda iterations first because they carry a loop over height in their subroutines
!---------------------------------------------------------------------------
  keepline = sline
  keepcont = scont
  call lambd_line(k_l,epsprm,Wm,gprof_re,gprof_im,jint,sline)
  call lambd_cont(s,eps_c,etaI_c,jrad,scont)
  IF(fl_lam.eq.1) THEN
   dline = sline - keepline
   do iz = 1,nz
    fac = k_l(iz)*Wm(iz)*epsprm(iz)/(epsprm(iz) + 1.d0)
    do i = 1,nwl
     add(iz,i) = jint(1,0,0,iz,i) - sline(1,0,0,iz,i) + fac*gprof_re(0,0,0,iz,i)
    enddo
   enddo
   call find_rho(k_l,etaI,lstar,dline,add,r2,r3,rhol)
  ENDIF
!---------------------------------------------------------------------------
  do iz=1,nz
   do i=1,nwl
    indx(i)=indxsto(iz,i)!reobtaining the stored indices 
    do ip=1,nwl
     m(i,ip)=msto(iz,i,ip)!and matrices
    enddo
    IF(fl_lam.ne.1) THEN
     rho(i) = jint(1,0,0,iz,i) + k_l(iz)*epsprm(iz)/(1.d0+epsprm(iz))*&
     Wm(iz)*gprof_re(0,0,0,iz,i) - sline(1,0,0,iz,i)
    ELSE
     rho(i) = rhol(iz,i)
    ENDIF
   enddo!Here S00 correction is found through the jacobi method
   call lubksb(m,nwl,indx,rho)! rho becomes Delta S^0_0
!-----------------------------------------------------------------------
! Start loop over frequency
!-----------------------------------------------------------------------
  do i = 1,nwl
   DS0 = dabs(rho(i))!variation of line S^0_0
   sline(1,0,0,iz,i) = sline(1,0,0,iz,i) + rho(i)
   ! space left for variations of other source functions, and the continuum part.
   DS2 = dabs(sline(1,2,0,iz,i) - keepline(1,2,0,iz,i))
   DS1 = dabs(sline(1,1,0,iz,i) - keepline(1,1,0,iz,i))
   DSR11 = dabs(sline(1,1,1,iz,i) - keepline(1,1,1,iz,i))
   DSI11 = dabs(sline(2,1,1,iz,i) - keepline(2,1,1,iz,i))
   DSR21 = dabs(sline(1,2,1,iz,i) - keepline(1,2,1,iz,i))
   DSI21 = dabs(sline(2,2,1,iz,i) - keepline(2,2,1,iz,i))
   DSR22 = dabs(sline(1,2,2,iz,i) - keepline(1,2,2,iz,i))
   DSI22 = dabs(sline(2,2,2,iz,i) - keepline(2,2,2,iz,i))
!-----------------------------------------------------------------------
! Calculate relative change of S00 (and leave space for others)
!-----------------------------------------------------------------------
   RDS0 = DS0/dabs(sline(1,0,0,iz,i))!; RDS2=DS2/dabs(S20(iz,i))
   RDS2 = DS2/dabs(sline(1,2,0,iz,i)); RDS1 = DS1/dabs(sline(1,1,0,iz,i))
   RDSR11 = DSR11/dabs(sline(1,1,1,iz,i)); RDSI11=DSI11/dabs(sline(2,1,1,iz,i))
   RDSR21 = DSR21/dabs(sline(1,2,1,iz,i)); RDSI21=DSI21/dabs(sline(2,2,1,iz,i))
   RDSR22 = DSR22/dabs(sline(1,2,2,iz,i)); RDSI22=DSI22/dabs(sline(2,2,2,iz,i))
!-----------------------------------------------------------------------
! Calculate maximum relative change of S00 (and, later on, others)
! Store frequency and height at which maximum relative change was found
!-----------------------------------------------------------------------
   if(RDS0.gt.MRDS0) then
    MRDS0=RDS0
    ht_MRDS0=ht(iz)
    wl_MRDS0=wl(i)
   endif
   if(RDS2.gt.MRDS2.and.ht(iz).gt.4.5d7) then
     MRDS2=RDS2
     ht_MRDS2=ht(iz)
     wl_MRDS2=wl(i)
   endif
   if(RDS1.gt.MRDS1.and.ht(iz).gt.4.5d7) MRDS1=RDS1 
   if(RDSR11.gt.MRDSR11.and.ht(iz).gt.-1.d8) MRDSR11=RDSR11
   if(RDSI11.gt.MRDSI11.and.ht(iz).gt.-1.d8) MRDSI11=RDSI11
   if(RDSR21.gt.MRDSR21.and.ht(iz).gt.-1.d8) MRDSR21=RDSR21
   if(RDSI21.gt.MRDSI21.and.ht(iz).gt.-1.d8) MRDSI21=RDSI21
   if(RDSR22.gt.MRDSR22.and.ht(iz).gt.-1.d8) MRDSR22=RDSR22
   if(RDSI22.gt.MRDSI22.and.ht(iz).gt.-1.d8) MRDSI22=RDSI22
!-----------------------------------------------------------------------
  enddo !i
  if(iz.eq.57) then
   fac = k_l(iz)*Wm(iz)*epsprm(iz)/(epsprm(iz) + 1.d0) 
  endif
 enddo !height loop
!-----------------------------------------------------------------------
 write(11,*) iter,MRDS0,MRDS1,MRDS2,MRDSR11,MRDSI11,MRDSR21,MRDSI21,&
  MRDSR22,MRDSI22
 write(41,*) iter,ht_MRDS0,wl_MRDS0,ht_MRDS2,wl_MRDS2
 print*, 'Maximum relative change:',MRDS0, MRDS2
!-----------------------------------------------------------------------
! If the maximum relative change of S00 is larger than that calculated 
! at the previous iteration, and more than 100 iterations have been 
! performed (so that initial 'oscillations' are over), the exit the 
! loop over iterations
!-----------------------------------------------------------------------
 if(iter.gt.18.and.MRDS0.ge.MRDS0prev) then
  fl1=1
  miter=iter
 endif
 MRDS0prev=MRDS0
 if(iter.gt.18.and.MRDS0.lt.5d-12) then
  print*,'!'
  fl1=1
  miter = iter
 endif
ENDIF
!=======================================================================
! 3) Perform formal solution transfer equation
!=======================================================================
 if(iter.eq.niter.or.fl1.eq.1) EXIT
!-----------------------------------------------------------------------
! Perform formal solution of the transfer equation and calculate new 
! values of the radiation tensors
!-----------------------------------------------------------------------
 print*,' -- Calculating formal solution RT equations...'
  call rt_solver(iter,Bpb,etaI,r,invk,mmat,sline,scont)
 print*, 'Done.'! Iteration:', iter
 if (df.eq.1) then !Disaster test
!  close(unit=66)
  print*, '--------------- DISASTER! ------------------'
  stop
 else
 endif!Disaster test
enddo!Iteration loop
close(unit=11)
close(unit=41)
!=======================================================================
! 4) Output radiation field and source functions.
!=======================================================================
 open(unit=14,file='radout_full.dat')
 write(14,*) nwl, nz
 do i=1,nwl
  write(14,*) wl(i)
 enddo
 do iz = 1,nz
  write(14,*) (jrad(1,0,0,iz,i),i=1,nwl)
 enddo
 do iz=1,nz
  write(14,*) (jrad(1,2,0,iz,i),i=1,nwl)
 enddo
 do iz=1,nz
  write(14,*) (jrad(1,1,0,iz,i),i=1,nwl)
 enddo
 do ik=1,2
  do iq=1,ik
   do iz=1,nz
    write(14,*) (jrad(1,ik,iq,iz,i),i=1,nwl)
   enddo
  do iz=1,nz
   write(14,*) (jrad(2,ik,iq,iz,i),i=1,nwl)
  enddo
 enddo!ik
enddo!iq
close(unit=14)
!===================================================================================
! 5) Perform last iteration calculations to find angle-dependent quantities such as
! emergent radiation and RT coefficients
!===================================================================================
do tg = 1,tgm
 if(tg .lt. 10) then
  format_string = "(A7,I1)"
  format_strb = "(A4,I1)"
  format_strc = "(A4,I1)"
  format_tst = "(A4,I1)"
  format_skh = "(A4,I1)"
 else
  format_string = "(A7,I2)"
  format_strb = "(A4,I2)"
  format_strc = "(A4,I2)"
  format_tst = "(A4,I2)"
  format_skh = "(A4,I2)"
 endif
 write(fileis,format_string) "stokes_",tg
 fileis = trim(fileis)//".res"
 write(fillin,format_string) "stolin_",tg
 fillin = trim(fillin)//".res"
 write(fileib,format_strb) "eta_",tg
 fileib = trim(fileib)//".res"
 write(fileic,format_strc) "rho_",tg
 fileic = trim(fileic)//".res"
  write(filskh,format_skh) "skh_",tg! Stokes parameters at each height
  filskh = trim(filskh)//".res"
 print*,' -- Calculating emergent radiation for mu=',mi(tg),'chi=',chi(tg)
 beta = DACOS(mi(tg))
 alpha = chi(tg)
 gamma = pi/2
 call TKQ2(-gamma,-beta,-alpha,0.d0,0.d0,0.d0,TR,TI)
! print*, 'TKQ2 works'
 allocate(eta_ref(0:3,1:nz,1:nwl),rho_ref(0:3,1:nz,1:nwl))
 call etasimple(k_l,TR,TI,eta_ref,rho_ref)
 if(fl_anom.ne.0) then
  do iz = 1,nz
   do i = 1,nwl
    do si = 1,3
     rho_ref(si,iz,i) = 0.d0
    enddo
   enddo
  enddo
 endif
! ------------------  ------------------  ------------------
 open(unit=12,file=trim(fileib))!file='eta_tg.res'
 write(12,*) mi(tg)
 do i = 1,nwl
  write(12,*) wl(i)
 enddo
 do si = 0,3
  do i=1,nwl
   if(si.eq.0) then
    write(12,*) (eta_ref(0,iz,i)+etaI_c(iz,i),iz = 1,nz)
   else
    write(12,*) (eta_ref(si,iz,i),iz = 1,nz)
   endif
  enddo
 enddo
 close(unit=12)
! ------------------  ------------------  ------------------
 open(unit=12,file=trim(fileic))!file='rho_tg.res'
 write(12,*) mi(tg)
 do i = 1,nwl
  write(12,*) wl(i)
 enddo
 do si = 0,3
  do i=1,nwl
    write(12,*) (rho_ref(si,iz,i),iz=1,nz)
  enddo!i 
 enddo!si
 close(unit=12)
!
 allocate(emer(0:3,1:nwl))
 allocate(emer_lin(0:3,1:nwl))
 call rt_emer(mi(tg),chi(tg),TR,TI,Bpb,etaI_c,sline,scont,eta_ref,rho_ref,emer)
!==============================================================================
! 6) Store emergent radiation (and Stokes parameters at each point)
!==============================================================================
 open(unit=12,file=trim(fileis))!file='stokes_1.res'
 write(12,*) nwl,miter
 do i=1,nwl
  wlvac = pc/(nu(i)*1.d-8)
  write(12,100) wl(i),emer(0,i),emer(1,i),emer(2,i),emer(3,i)
 enddo
 close(unit=12)
!
 deallocate(eta_ref,rho_ref)
 deallocate(emer)
 deallocate(emer_lin)
enddo!tg loop (in order to get 4 directions for emergent radiation)
!
call cpu_time(timf)
print*,'Time:',(timf-timi)/3600,'hours'
!print*, 'BESSER'
!------------------------------------------------------------------------
 100 format(e18.11,2x,4(2x,e13.6))
 102 format(2(e10.3,2x))
 105 format(e18.11,2x)
 333 format(e11.4)
!------------------------------------------------------------------------
END PROGRAM hz_rta
