MODULE rt_formal
!-----------------------------------------------------------------
! Ernest Alsina, 2015
! Contains subroutines for the calculation of the lambda star
! operators, for the formal solution routine, and for calculating
! the emergent Stokes profiles.
!-----------------------------------------------------------------
USE parameters
USE global_vars
USE math
!-----------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------

!=================================================================
SUBROUTINE lam_star(etaI,eta_re,rho_re,invk,mmat,lstar)
DOUBLE PRECISION, INTENT(IN) :: etaI(nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: eta_re(0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: rho_re(0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(OUT) :: invk(0:3,0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(OUT) :: mmat(0:3,0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(OUT) :: lstar(0:3,0:3,nz,nwl,ndir,ndirx)
INTEGER :: i, im, ax, k, kO, k1, kdel, kU, kD,si,sj,sk1,sk2
DOUBLE PRECISION :: CM, CO, CP, ClM, ClO
DOUBLE PRECISION :: dlU, dlD,dtU,dtD, exu
DOUBLE PRECISION :: mi, mi2, alph,den,b
DOUBLE PRECISION :: tst
DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: av
DOUBLE PRECISION, DIMENSION(1:nz,0:3,0:3) :: Iprev
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: kcur,kprev,opm,opo,omat
DOUBLE PRECISION, PARAMETER :: trick = 1.d-3,rtiny = 1.d-30
LOGICAL :: flov
!------------------------------------------------------------------
! Initializing Iprev (Stokes vector at the previous point due to a
! source function unit pulse Sj at point O). Needed for the
! ray coming from that point due to the source function at point O
!-----------------------------------------------------------------
do k = 1,nz
 do si = 0,3
  do sj = 0,3
   Iprev(k,si,sj) = 0.d0
  enddo
 enddo!si
enddo !k
!---------------------------------------------------------------
! Begin frequency loop
!---------------------------------------------------------------
do i = 1, nwl
!---------------------------------------------------------------
! Start loops over polar angle (using Gauss-Legendre quadrature)
! and azimuthal angle (trapezoidal quadrature). The weights are
! not considered here because we aren't performing a numerical
! integration here. For the polar angle, the mu values are found
! so they go from -1 to 1
!---------------------------------------------------------------
 do im = 1,ndir !polar angle goes from 0 to pi
  mi = mu(im)
  mi2 = mi*mi
  do ax = 1,ndirx ! goes from 0 to 2*pi
   alph = az(ax)
!---------------------------------------------------------------
! Start loop over heights. Range will depend on the polar angle;
! i.e. whether we are going inward or outward
!---------------------------------------------------------------
   if(mi.le.0) then !INWARD
    kO = 1
    k1 = nz
    kdel = 1
   else              !OUTWARD
    kO = nz
    k1 = 1
    kdel = -1
   endif
!--------------------------------------------------------------------
! Initializing Iprev (Stokes vector at the previous point due to a
! source function unit pulse Sj at point O). Needed for the
! ray coming from that point due to the source function at point O
!--------------------------------------------------------------------
   do k = 1, nz
    do si = 0,3
     do sj = 0,3
      Iprev(k,si,sj) = 0.d0
     enddo
    enddo!si
   enddo !k
!--------------------------------------------------------------------
! Now start the loop over height points
!--------------------------------------------------------------------
   do k = kO+kdel, k1, kdel
    ! UPWIND
    kU = k-kdel
    dlU = DABS((ht(kU) - ht(k))/mi)
    dtU = 0.5d0*(etaI(k,i,im,ax) + etaI(kU,i,im,ax))*dlU !Change in optical depth
    ! DOWNWIND
    if(k.ne.k1) then
     kD = k+kdel
     dlD = DABS((ht(k) - ht(kD))/mi)
     dtD = 0.5d0*(etaI(k,i,im,ax) + etaI(kD,i,im,ax))*dlD !Change in optical depth
    else !If we are on the surface kD is undefined (out of the array)
     dlD = dlU
     dtD = dtU
    endif
    ! Calculate e^-{delta tau_M}
    if(dtU.ge.trick) then
     exu = DEXP(-dtU)
    else
     exu =  1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
      (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
      (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
      (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
      (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
      (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0
    endif! dtU
    if(k.eq.kO + kdel) then
     kprev = 0.d0
     den = etaI(kO,i,im,ax)
     kprev(0,1) = eta_re(1,kO,i,im,ax)/den! Line contribution is the only one there is
     kprev(0,2) = eta_re(2,kO,i,im,ax)/den! since continuum part is zero 
     kprev(0,3) = eta_re(3,kO,i,im,ax)/den
     kprev(1,0) = eta_re(1,kO,i,im,ax)/den
     kprev(1,2) = rho_re(3,kO,i,im,ax)/den
     kprev(1,3) = -rho_re(2,kO,i,im,ax)/den
     kprev(2,0) = eta_re(2,kO,i,im,ax)/den
     kprev(2,1) = -rho_re(3,kO,i,im,ax)/den
     kprev(2,3) = rho_re(1,kO,i,im,ax)/den
     kprev(3,0) = eta_re(3,kO,i,im,ax)/den
     kprev(3,1) = rho_re(2,kO,i,im,ax)/den
     kprev(3,2) = -rho_re(1,kO,i,im,ax)/den
    else
     kprev = kcur!K_M
    endif! k = k+kdel
    kcur = 0.d0
    den = etaI(k,i,im,ax)
    kcur(0,1) = eta_re(1,k,i,im,ax)/den! Line contribution is the only one there is
    kcur(0,2) = eta_re(2,k,i,im,ax)/den! since continuum part is zero 
    kcur(0,3) = eta_re(3,k,i,im,ax)/den
    kcur(1,0) = eta_re(1,k,i,im,ax)/den
    kcur(1,2) = rho_re(3,k,i,im,ax)/den
    kcur(1,3) = -rho_re(2,k,i,im,ax)/den
    kcur(2,0) = eta_re(2,k,i,im,ax)/den
    kcur(2,1) = -rho_re(3,k,i,im,ax)/den
    kcur(2,3) = rho_re(1,k,i,im,ax)/den
    kcur(3,0) = eta_re(3,k,i,im,ax)/den
    kcur(3,1) = rho_re(2,k,i,im,ax)/den
    kcur(3,2) = -rho_re(1,k,i,im,ax)/den
!-----------------------------------------------------------------    
    call LINEAR(exu,dtU,ClO,ClM)
    opm = exu*iden - ClM*kprev![1*e^-dtU - Psi_M K'_M] 
    opo = iden + ClO*kcur
    call MatInv(opo) ![1 + Psi_O' K'_O]
    do si = 0,3 !!!TEST!!!
     do sj = 0,3
      if(dabs(opo(si,sj)).le.rtiny) opo(si,sj) = 0.d0
     enddo!in2
    enddo!in1
    invk(:,:,k,i,im,ax) = opo(:,:)
!------------------------------------------------------------------
    call MatMat(opo,opm,omat)
    do si = 0,3 !!!TEST!!!
     do sj = 0,3
      if(dabs(omat(si,sj)).le.rtiny) omat(si,sj) = 0.d0
     enddo!sj
    enddo!si
    mmat(:,:,k,i,im,ax) = omat(:,:)![1 + Phi_O*K']^-1 [1*e^Dtau - Phi_M*K_M]
!
    if(k.ne.k1) then!parabolic
    call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
!Contribution from "previous" point O, given unit pulse at P
     do si = 0,3
      do sj = 0,3
       call Beziercm(dtU,dtD,0.d0,0.d0,1.d0,cpm,flov)
       Iprev(k,si,sj) = opo(si,sj)*(omega_c*cpm)! I_i obtained from a unit pulse S_j  
      enddo
     enddo
     do sj = 0,3!S_j unit pulse
     call Beziercm(dtU,dtD,0.d0,1.d0,0.d0,cpm,flov)
      do si = 0,3 !Gives I_0,i
       av(si,sj) = opo(si,sj)*(omega_o + omega_c*cpm)
       b = 0.d0!Contribution to I_0,i due to the I_M produced by a unit pulse
               ! S_j at the current point
        do sk2 = 0,3
          b = b + omat(si,sk2)*Iprev(k-kdel,sk2,sj)
        enddo
       lstar(si,sj,k,i,im,ax) = av(si,sj) + b! ??
      enddo
     enddo 
    else !point at the surface: we can't use CP, which must stay zero
     do sj = 0,3!S_j unit pulse
      do si = 0,3 !Gives I_0,i
       av(si,sj) = opo(si,sj)*ClO
       b = 0.d0!Contribution to I_0,i due to the I_M produced by a unit pulse S_j at the current point
       do sk2 = 0,3
        b = b + omat(si,sk2)*Iprev(k-kdel,sk2,sj)
       enddo
       lstar(si,sj,k,i,im,ax) = av(si,sj) + b! ??
      enddo!si
     enddo!sj    
    endif!Any point
!    if(im.eq.9.and.ax.eq.1.and.i.eq.82) then
!     write(61,*) "Height:",k
!     write(62,*) "Height:",k
!     write(63,*) "Height:",k
!     do si = 0,3
!      write(61,*) (lstar(si,sj,k,i,im,ax),sj=0,3)
!      write(62,*) (invk(si,sj,k,i,im,ax),sj=0,3)
!      write(63,*) (mmat(si,sj,k,i,im,ax),sj=0,3)
!     enddo
!    endif
!-------------------------------------------------------------------------------
! Using linear interpolation, there is no contribution to the upwind due to S, 
! and so no I coming from it. Therefore, no PSP*exu
!-------------------------------------------------------------------------------
!    lstarp(k,i,im,ax) = PLINEA(dtU,1.d0,0.d0)
   enddo!k 
  enddo !ax
 enddo !m
enddo !i
!close(21)
!close(unit=61)
!close(unit=62)
!close(unit=63)
!-------------------------------------------------------------------------------
 return
!-------------------------------------------------------------------------------
END SUBROUTINE lam_star
!=================================================================

!=================================================================
SUBROUTINE rt_solver(iter,Bpb,etaI,r,invk,mmat,sline,scont)
! Formal solution routine using the DELOPAR method.
INTEGER :: iter
DOUBLE PRECISION, INTENT(IN) :: Bpb(nwl)
DOUBLE PRECISION, INTENT(IN) ::  etaI(nz,nwl,ndir,ndirx),r(nz,nwl,ndir,ndirx)!Keep the r for a test. 
DOUBLE PRECISION, INTENT(IN) :: invk(0:3,0:3,nz,nwl,ndir,ndirx) ![1 + Phi_0*K']^-1
DOUBLE PRECISION, INTENT(IN) :: mmat(0:3,0:3,nz,nwl,ndir,ndirx) ![1 + Phi_0*K']^-1 [1e^-dt - Phu_M*K']
DOUBLE PRECISION, INTENT(IN) :: sline(2,0:2,0:2,nz,nwl), scont(2,0:2,0:2,nz,nwl)
INTEGER :: k, k1, kO, kdel, i, iz, im, ax, kU, kD
INTEGER :: ik, iq, in
DOUBLE PRECISION :: dlU, dtU, dlD, dtD, exu
DOUBLE PRECISION :: ClM, ClO, CM, CO, CP
DOUBLE PRECISION :: a,b,c, eM, eO, eP!, elM, elO, elP
DOUBLE PRECISION :: mi, mi2
DOUBLE PRECISION :: cre,cim,lre,lim
DOUBLE PRECISION :: ctre, wgt, ar, ai
DOUBLE PRECISION, DIMENSION(0:3) :: SM, SO, SP
DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
DOUBLE PRECISION, DIMENSION(0:3) :: stokes, sto_prev, svector
DOUBLE PRECISION, DIMENSION(0:3) :: avec, bvec !not an allocatable
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: amat, bmat
DOUBLE PRECISION, PARAMETER :: trick = 1.d-3
CHARACTER(len=1024) :: filop,format_op
LOGICAL :: flov
!Initializing the Stokes parameters as zero:
do in=0,3!TEST
 stokes(in) = 0.d0
 sto_prev(in) = 0.d0
 svector(in) = 0.d0
enddo!in
! Set the radiation tensor initially to zero
jrad = 0.d0
! Initialize constant
ctre = 0.d0
df = 0!disaster flag
if(iter.eq.1) open(unit=29, file='rt_test.dat')
!TEST
!---------------------------------------------------------------
! Begin frequency loop
!---------------------------------------------------------------
 do i = 1, nwl
!---------------------------------------------------------------
! Start loops over polar angle (using Gauss-Legendre quadrature)
! and azimuthal angle (trapezoidal quadrature). The weights are
! not considered here because we aren't performing a numerical
! integration here. For the polar angle, the mu values are found
! so they go from -1 to 1
!---------------------------------------------------------------
  do im = 1, ndir !polar angle goes from 0 to pi
   mi = mu(im)
   mi2 = mi*mi
   do ax = 1, ndirx ! goes from 0 to 2*pi
   !alph = az(ax)
    wgt = wtdir(im)*waz(ax)/(4.d0*pi)
!---------------------------------------------------------------
! Start loop over heights. Range will depend on the polar angle;
! i.e. whether we are going inward or outward
!---------------------------------------------------------------
    if(mi.le.0) then !INWARD
     kO = 1
     k1 = nz
     kdel = 1
    else              !OUTWARD
     kO = nz
     k1 = 1
     kdel = -1
    endif
    !First height point
    do k = kO + kdel, k1, kdel
    ! UPWIND
     kU = k-kdel
     dlU = DABS((ht(kU) - ht(k))/mi)
     dtU = 0.5d0*(etaI(k,i,im,ax) + etaI(kU,i,im,ax))*dlU !Change in optical depth
    ! DOWNWIND
     if(k.ne.k1) then
      kD = k+kdel
      dlD = DABS((ht(k) - ht(kD))/mi)
      dtD = 0.5d0*(etaI(k,i,im,ax) + etaI(kD,i,im,ax))*dlD !Change in optical depth
     else !If we are on the surface kD is undefined (out of the array)
      dlD = dlU
      dtD = dtU
     endif
    ! Calculate e^-{delta tau_M}
     if(dtU.ge.trick) then
      exu = DEXP(-dtU)
     else
      exu =  1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
       (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
       (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0
     endif

     IF(k.ne.k1) THEN
     !ecP = etaI_c(k+kdel,i)
      eP = etaI(k+kdel,i,im,ax)
      IF(k.eq.kO+kdel) THEN !First point
       ! ecO = etaI_c(k,i)
       !ecM = etaI_c(k-kdel,i)
       eO = etaI(k,i,im,ax)
       eM = etaI(k-kdel,i,im,ax)
       do in = 0,3
        a = 0.d0; b = 0.d0; c = 0.d0
        do ik = 0, 2
         do iq = 0, ik
!------------------------------------------------------------------------------------------------
! Calculating total S_i source functions at upwind (M), current (O) and downwind (P) points. 
! It is insider a loop over Stokes vector components 0 through 3
!------------------------------------------------------------------------------------------------
         if(iq.eq.0) then
          ctre = 1.d0 !ONLY REAL PART IF Q=0
         else
   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
          ctre = 2.d0
         endif
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
         a = a + ctre*(lre + cre)/eP!P point is downwind
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k,i),sline(2,ik,iq,k,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k,i),scont(2,ik,iq,k,i),cre,cim)
         b =  b + ctre*(lre + cre)/eO!current point
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k-kdel,i),sline(2,ik,iq,k-kdel,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k-kdel,i),scont(2,ik,iq,k-kdel,i),cre,cim)
         c = c + ctre*(lre+cre)/eM !M is the upwind point
         enddo !iq
        enddo!ik
        SP(in) = a
        SO(in) = b
        SM(in) = c
      enddo!in
      sto_prev = 0.d0 !Boundary point
      if(mi.ge.0) sto_prev(0) = Bpb(i)!Bottom point. Assumed thermalized.
       do ik = 0,2
        do iq = 0,ik!At boundary conditions we assume no (incoming) polarization, neither from the bottom nor top. Therefore
         jrad(1,ik,iq,kO,i) = jrad(1,ik,iq,kO,i) + wgt*TTR(ik,iq,0,im,ax)*sto_prev(0)!we only need 0 Stokes index
         if(iq.gt.0) jrad(2,ik,iq,kO,i) = jrad(2,ik,iq,kO,i) + wgt*TTI(ik,iq,0,im,ax)*sto_prev(0)
        enddo ! iq
       enddo ! ik
      ELSE ! Middle points
       !HERE we only need to look for SP, since SO and SM have already been found in previous k points
       do in = 0, 3
        a = 0.d0!; b = 0.d0; c = 0.d0
        do ik = 0, 2
         do iq = 0, ik
          if(iq.eq.0) then
           ctre = 1.d0 !ONLY REAL PART IF Q=0
          else
! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
           ctre = 2.d0
          endif
          call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
          sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
          call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
          scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
          a = a + ctre*(lre+cre)/eP
         enddo !iq
        enddo!ik
        SP(in) = a
       enddo!in
      ENDIF
!----------------------------------------------------------------------------------------------
! BESSER routine itself
!----------------------------------------------------------------------------------------------
!DELOPAR
!      call PARAGUASB(dtU,dtD,CM,CO,CP)!For tests
!BESSER
      call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
      do in = 0,3
!DELOPAR
!       svector(in) = CM*SM(in) + CO*SO(in) + CP*SP(in)
!BESSER
       call Beziercm(dtU,dtD,SM(in),SO(in),SP(in),cpm,flov)
       svector(in) = omega_m*SM(in) + omega_o*SO(in) + omega_c*cpm
!       if(in.eq.0.and.im.eq.3.and.ax.eq.3.and.i.eq.54) then
!        write(*,*) k
!!        write(*,*) SM(0),SO(0),SP(0)
!        write(*,*) omega_m,SM(0),omega_o,SO(0),omega_c,cpm
!       endif
      enddo!in
      amat(:,:) = mmat(:,:,k,i,im,ax)
      bmat(:,:) = invk(:,:,k,i,im,ax)
      call MatVec(amat,sto_prev,avec) ![1 + lPhi_O K'_O]^-1 [e^-dtM -  lPhi_M K'_M] I_M
      call MatVec(bmat,svector,bvec) ![1 + lPhi_O K'_O]^-1 (Phi_M S_M + Phi_O S_O + Phi_P S_P)
! TEST: RT control
!      if(i.eq.99.and.im.eq.13.and.ax.eq.3) then
!       write(*,*) "Height:",k
!       write(*,*) CM,CO,CP
!       write(*,*) SM(0),SO(0),SP(0)
!       write(*,*) sto_prev(0),svector(0),avec(0)+bvec(0)
!      endif
!     if(im.eq.3.and.ax.eq.3.and.i.eq.54) then
!      write(*,*) sto_prev(0),svector(0),avec(0)+bvec(0)
!     endif
     do in = 0,3!TEST
       stokes(in) = avec(in) + bvec(in)
       sto_prev(in) = stokes(in)
       SM(in) = SO(in)
       SO(in) = SP(in)
       !ecM = ecO
       !ecO = ecP
      ! eM = eO
      ! eO = eP
      enddo !in
! ----------- TEST ----------------
!      if(im.eq.8.and.ax.eq.4.and.i.eq.78) then
!       write(*,*) k,stokes(0),stokes(1),stokes(2),stokes(3)
!      endif
! ----------- TEST ----------------
      eM = eO
      eO = eP
     ELSE !k=k1. LINEAR INTERPOLATION INSTEAD OF PARABOLIC
      !We already have SM and SO, and we don't need SP
      call LINEAR(exu,dtU,ClO,ClM)
      do in=0,3
       svector(in) = ClM*SM(in) + ClO*SO(in)
      enddo!in
      amat(:,:) = mmat(:,:,k,i,im,ax)
      bmat(:,:) = invk(:,:,k,i,im,ax)
      call MatVec(amat,sto_prev,avec)
      call MatVec(bmat,svector,bvec)
      do in = 0,3
       stokes(in) = avec(in) + bvec(in)
      enddo
! ----------- TEST ----------------
!       if(im.eq.8.and.ax.eq.4.and.i.eq.78) then
!       write(*,*) k,stokes(0),stokes(1),stokes(2),stokes(3)
!      endif
! ----------- TEST ----------------
! TEST: RT control (Last point)
!      if(i.eq.99.and.im.eq.13.and.ax.eq.3) then
!       write(*,*) "Height:",k
!       write(*,*) ClM,ClO
!       write(*,*) SM(0),SO(0)
!       write(*,*) sto_prev(0),svector(0),avec(0)+bvec(0) 
!      endif
! TEST: RT control (last point)
     ENDIF
!----------------------------------------------------------------------
! Control
!----------------------------------------------------------------------
     IF(stokes(0).lt.0.d0) THEN
       print*,'DISASTER! I .lt. 0 at k=', k, stokes(0), sline(1,0,0,k,i), sline(1,1,0,k,i) ,sline(1,2,0,k,i)
      df = 1!flag that indicates something going wrong in the rt_solver
     ELSE
     ENDIF
!----------------------------------------------------------------------
! Calculate the radiation tensors (multipolar components)
!----------------------------------------------------------------------
     do ik = 0,2!
      do iq = 0,ik
       ar = 0.d0
       ai = 0.d0
       do in = 0,3 !summing over the stokes parameters
        ar = ar + TTR(ik,iq,in,im,ax)*stokes(in)
        !IF(ik.eq.0) THEN
        ! if(k.eq.60) print*, i,im,ax,stokes(0),stokes(1),stokes(2),stokes(3)
        !ENDIF
        if(iq.gt.0) ai = ai + TTI(ik,iq,in,im,ax)*stokes(in)
       enddo!in
       jrad(1,ik,iq,k,i) = jrad(1,ik,iq,k,i) + wgt*ar
       if(iq.gt.0) jrad(2,ik,iq,k,i) = jrad(2,ik,iq,k,i) + wgt*ai
    !   IF(k.eq.60) THEN
    !    print*, jrad(1,0,0,k,i)
    !   ENDIF
      enddo!iq
     enddo!ik
!---RADIATION TENSOR FOUND
    enddo!k End height grid
   enddo!ax End azimuthal angle grid
  enddo!im End polar angle grid
 enddo!i End frequency grid
close(unit=19)
!stop
!----------------------------------------------------------------------
 return
!----------------------------------------------------------------------
END SUBROUTINE rt_solver
!===========================================================================


!===========================================================================
SUBROUTINE rt_emer(mi,chi,TR,TI,Bpb,etaI_c,sline,scont,eta_ref,rho_ref,emer)
DOUBLE PRECISION, INTENT(IN) :: mi, chi
DOUBLE PRECISION, INTENT(IN) :: TR(0:2,-2:2,0:3), TI(0:2,-2:2,0:3)
DOUBLE PRECISION, INTENT(IN) :: Bpb(nwl)
!DOUBLE PRECISION, INTENT(IN) :: sr(nz,nwl), setaI(nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: etaI_c(nz,nwl)!Continuum part doesn't depend on angle
!DOUBLE PRECISION, INTENT(IN) :: sinvk(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1
!DOUBLE PRECISION, INTENT(IN) :: smmat(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1 [1e^-dt - Phu_M*K']
DOUBLE PRECISION, INTENT(IN) :: eta_ref(0:3,nz,nwl), rho_ref(0:3,nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: sline(2,0:2,0:2,nz,nwl), scont(2,0:2,0:2,nz,nwl)
DOUBLE PRECISION, INTENT(OUT) :: emer(0:3,nwl)
INTEGER :: k, k1, kO, kdel, i, iz, im, ax, kU, kD
INTEGER :: ik, iq, in
DOUBLE PRECISION :: dlU, dtU, dlD, dtD, exu
DOUBLE PRECISION :: ClM, ClO, CM, CO, CP
DOUBLE PRECISION :: a,b,c, ecM, ecO, ecP, elM, elO, elP
DOUBLE PRECISION :: cre, cim, lre, lim
DOUBLE PRECISION :: ctre, wgt, ar, ai
DOUBLE PRECISION :: rP, rM, rO
DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
DOUBLE PRECISION, DIMENSION(0:3) :: SM, SO, SP,cptst
DOUBLE PRECISION, DIMENSION(0:3) :: stokes, sto_prev, svector
DOUBLE PRECISION, DIMENSION(0:3) :: avec, bvec! not allocatable anymore
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: bmat !not allocatable anymore
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: amat, omat, tmat, smat
DOUBLE PRECISION, PARAMETER :: trick = 1.d-3
LOGICAL :: flov
!Initializing the Stokes parameters as zero:
!if(mi.eq.0.16d0) open(unit=31,file='itef_last.res')
stokes = 0.d0
sto_prev = 0.d0
svector = 0.d0
emer = 0.d0
! Initialize constant
ctre = 0.d0
df = 0!disaster flag
!---------------------------------------------------------------
! Begin frequency loop
!---------------------------------------------------------------
do i = 1, nwl
!---------------------------------------------------------------
! Start loops over polar angle (using Gauss-Legendre quadrature)
! and azimuthal angle (trapezoidal quadrature). The weights are
! not considered here because we aren't performing a numerical
! integration here. For the polar angle, the mu values are found
! so they go from -1 to 1
!---------------------------------------------------------------
if (mi .le. 0) then !INWARD
print*, 'We should be considering an emergent ray'
stop
kO = 1
k1 = nz
kdel = 1
else              !OUTWARD
kO = nz
k1 = 1
kdel = -1
endif
!First height point
do k = kO + kdel, k1, kdel
! UPWIND
kU = k-kdel
dlU = DABS((ht(kU) - ht(k))/mi)
dtU = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kU,i) + eta_ref(0,kU,i))*dlU !Change in optical depth
! DOWNWIND
if(k.ne.k1) then
kD = k+kdel
dlD = DABS((ht(k) - ht(kD))/mi)
dtD = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kD,i) + eta_ref(0,kD,i))*dlD !Change in optical depth
else !If we are on the surface kD is undefined (out of the array)
dlD = dlU
dtD = dtU
endif
! Calculate e^-{delta tau_M}
if(dtU.ge.trick) then
exu = DEXP(-dtU)
else
exu = 1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
  (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
  (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0
endif
!------------------------------------------------------------------------------
! Calculating the K_0 matrices and K_M at the first point
!------------------------------------------------------------------------------
if(k.eq.kO+kdel) then
smat = 0.d0
smat(0,1) = eta_ref(1,kO,i)! Line contribution is the only one there is
smat(0,2) = eta_ref(2,kO,i)! since continuum part is zero 
smat(0,3) = eta_ref(3,kO,i)
smat(1,0) = eta_ref(1,kO,i)
smat(1,2) = rho_ref(3,kO,i)
smat(1,3) = -rho_ref(2,kO,i)
smat(2,0) = eta_ref(2,kO,i)
smat(2,1) = -rho_ref(3,kO,i)
smat(2,3) = rho_ref(1,kO,i)
smat(3,0) = eta_ref(3,kO,i)
smat(3,1) = rho_ref(2,kO,i)
smat(3,2) = -rho_ref(1,kO,i)
smat = smat/(etaI_c(kO,i) + eta_ref(0,kO,i))
else
smat = tmat!K_M
endif
tmat = 0.d0
tmat(0,1) = eta_ref(1,k,i)! Line contribution is the only one there is
tmat(0,2) = eta_ref(2,k,i)! since continuum part is zero 
tmat(0,3) = eta_ref(3,k,i)
tmat(1,0) = eta_ref(1,k,i)
tmat(1,2) = rho_ref(3,k,i)
tmat(1,3) = -rho_ref(2,k,i)
tmat(2,0) = eta_ref(2,k,i)
tmat(2,1) = -rho_ref(3,k,i)
tmat(2,3) = rho_ref(1,k,i)
tmat(3,0) = eta_ref(3,k,i)
tmat(3,1) = rho_ref(2,k,i)
tmat(3,2) = -rho_ref(1,k,i)
tmat = tmat/(etaI_c(k,i)+eta_ref(0,k,i))!K_O
call LINEAR(exu,dtU,ClO,ClM)
bmat = iden + ClO*tmat
call MatInv(bmat)
omat = exu*iden - ClM*smat
call MatMat(bmat,omat,amat)

IF(k.ne.k1) THEN
ecP = etaI_c(k+kdel,i)
elP = eta_ref(0,k+kdel,i)! - ecP
rP = elP/(ecP+elP)
  IF(k.eq.kO+kdel) THEN !First point
  ecO = etaI_c(k,i)
  ecM = etaI_c(k-kdel,i)
  elO = eta_ref(0,k,i)
  elM = eta_ref(0,k-kdel,i)
  rO = elO/(ecO + elO)
  rM = elM/(ecM + elM)
  do in = 0,3
  a = 0.d0; b = 0.d0; c = 0.d0
  do ik = 0, 2
  do iq = 0, ik
  if(iq.eq.0) then
   ctre = 1.d0 !ONLY REAL PART IF Q=0
  else
   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART, due to +Q and -Q terms
   ctre = 2.d0
  endif
  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
   a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k,i),sline(2,ik,iq,k,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k,i),scont(2,ik,iq,k,i),cre,cim)
   b =  b + ctre*rO*lre/elO + ctre*(1.d0 - rO)*cre/ecO
  call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k-kdel,i),sline(2,ik,iq,k-kdel,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k-kdel,i),scont(2,ik,iq,k-kdel,i),cre,cim)
   c = c + ctre*rM*lre/elM + ctre*(1.d0 - rM)*cre/ecM
  enddo !iq
  enddo!ik
  SP(in) = a
  SO(in) = b
  SM(in) = c
  enddo!in
  sto_prev = 0.d0 !Boundary point
  sto_prev(0) = Bpb(i)!Bottom point. Assumed thermalized.
  ELSE ! Middle points
   !HERE we only need to look for SP, since SO and SM have already been found in previous k points
  do in = 0,3
  a = 0.d0; b = 0.d0; c = 0.d0
  do ik = 0, 2
   do iq = 0, ik
    if(iq.eq.0) then
     ctre = 1.d0 !ONLY REAL PART IF Q=0
    else
   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
     ctre = 2.d0
    endif
   call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
   call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
    a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
   enddo !iq
  enddo!ik
  SP(in) = a
  enddo!in
  ENDIF
!DELOPAR
!call PARAGUASB(dtU,dtD,CM,CO,CP)
!BESSER
call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
  svector = CM*SM + CO*SO + CP*SP
 do in = 0,3
!DELOPAR
!  svector(in) = CM*SM(in) + CO*SO(in) + CP*SP(in)
!BESSER
  call Beziercm(dtU,dtD,SM(in),SO(in),SP(in),cpm,flov)
!  cptst(in) = cpm
  svector(in) = omega_m*SM(in) + omega_o*SO(in) + omega_c*cpm
 enddo
 ! write(34,*) k,i,CM,CO,CP,omega_m,omega_o,omega_c,SM(0),SO(0),SP(0),cptst(0),&
 !            SM(1),SO(1),SP(1),cptst(1)
  call MatVec(amat,sto_prev,avec)
  call MatVec(bmat,svector,bvec)
 stokes = avec + bvec
!  deallocate(avec)
!  deallocate(bvec)
 sto_prev = stokes
!-------------TEST--------------------------------------
SM = SO
SO = SP
ecM = ecO
ecO = ecP
elM = elO
elO = elP
rM = rO
rO = rP
ELSE !k=k1. LINEAR INTERPOLATION INSTEAD OF PARABOLIC
!We already have SM and SO, and we don't need SP
 svector = ClM*SM + ClO*SO
!  amat(:,:) = smmat(:,:,k,i)
!  bmat(:,:) = sinvk(:,:,k,i)
call MatVec(amat,sto_prev,avec)
call MatVec(bmat,svector,bvec)
 stokes = avec + bvec
 do in = 0,3
 emer(in,i) = stokes(in)
 enddo!in
ENDIF
!----------------------------------------------------------------------
! Control
!----------------------------------------------------------------------
enddo!k grid
enddo !i
!----------------------------------------------------------------------
 return
!----------------------------------------------------------------------
END SUBROUTINE rt_emer
!===========================================================================

!===========================================================================
SUBROUTINE LINEAR(exu,dtM,CO,CM)
!------------------------------------------------------------------------------
! Based on routines by Javier Trujillo Bueno.
! Used to obtain CO (Phi'_O) and CM (Phi'_M): the factors for linear interpolation
!------------------------------------------------------------------------------
  DOUBLE PRECISION :: PLINEA
  DOUBLE PRECISION, INTENT(IN) :: dtM,exu
  DOUBLE PRECISION, INTENT(OUT) :: CO, CM
  DOUBLE PRECISION :: U0, U1, D2
  DOUBLE PRECISION, PARAMETER :: trick=0.11d0
!------------------------------------------------------------------------------
if(dtM .gt. trick) then
 CM = ((1.d0-exu*(1.d0+dtM))/dtM)
 CO = ((exu+dtM-1.d0)/(dtM))
else
 CM = ((dtM*(dtM*(dtM*(dtM*(dtM*(dtM*((63d0-8d0*dtM)*dtM-432d0)+2520d0)-&
       12096d0)+45360d0)-120960d0)+181440d0))/362880d0)
 CO = ((dtM*(dtM*(dtM*(dtM*(dtM*(dtM*((9d0-dtM)*dtM-72d0)+504d0)-&
       3024d0)+15120d0)-60480d0)+181440d0))/362880d0)
endif
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END SUBROUTINE LINEAR
!===========================================================================


!===========================================================================
SUBROUTINE OMEGA(exu,dtM,omega_m,omega_o,omega_c)
!------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: exu,dtM
DOUBLE PRECISION, INTENT(OUT) :: omega_m,omega_o,omega_c
DOUBLE PRECISION, PARAMETER :: trick1 = 0.14d0,trick2 = 0.18d0
!DOUBLE PRECISION, PARAMETER :: trick1 = 1.d-3,trick2 = 1.d-3
!------------------------------------------------------------------------------
if(dtM.ge.trick1) then
 omega_m = ((2.d0-exu*(dtM*dtM+2.d0*dtM+2.d0))/(dtM*dtM))
else
 omega_m = ((dtM*(dtM*(dtM*(dtM*(dtM*(dtM*((140d0-18d0*dtM)*dtM-945d0)+5400d0)-&
  25200d0)+90720d0)-226800d0)+302400d0))/907200d0)
endif
if(dtM.ge.trick2) then
 omega_o = (1d0-2d0*(dtM+exu-1d0)/(dtM*dtM))
 omega_c = (2d0*(dtM-2d0+exu*(dtM+2d0))/(dtM*dtM))
else
 omega_o = ((dtM*(dtM*(dtM*(dtM*(dtM*(dtM*((10d0-dtM)*dtM-90d0)+720d0)-5040d0)+&
  30240d0)-151200d0)+604800d0))/1814400d0)
 omega_c = ((dtM*(dtM*(dtM*(dtM*(dtM*(dtM*((35d0-4d0*dtM)*dtM-270d0)+1800d0)-&
   10080d0)+45360d0)-151200d0)+302400d0))/907200d0)
endif
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END SUBROUTINE OMEGA
!===========================================================================

!===========================================================================
SUBROUTINE Beziercm(h0,h1,ym,yo,yp,cm,flov)
!------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: h0,h1,ym,yo,yp
DOUBLE PRECISION, INTENT(OUT) :: cm
LOGICAL,INTENT(OUT) :: flov!TEST
DOUBLE PRECISION :: dm,dp,yder,c0,c1,c2
LOGICAL :: cond0,cond1
!------------------------------------------------------------------------------
if(h0.gt.0d0 .and. h1.gt.0.d0) then!optical distances larger than zero
 dm = (yo-ym)/h0
 dp = (yp-yo)/h1
else
 cm = yo
 return
endif
!******************************************************************************
if(dm*dp .le. 0.d0) then !Not monotonic
 cm = yo
 return
endif
!******************************************************************************
yder = (h0*dp+h1*dm)/(h0+h1)!initial estimate for derivative at O
c0 = yo - 0.5d0*h0*yder! Bezier control point M
c1 = yo + 0.5d0*h1*yder! Bezier control point P
cond0 = interv(c0,ym,yo)!TRUE if c0 is in the [ym,yo] interval
cond1 = interv(c1,yo,yp)!TRUE if c1 is in the [yo,yp] interval
!******************************************************************************
flov = .true.!TEST
if(cond0.and.cond1) then !No corrections needed
 cm = c0
 flov = .false.!TEST
 return
elseif(.not.cond0) then
 cm = correctyab(c0,ym,yo) !Correct overshoot at cm, then exit
elseif(.not.cond1) then
 c2 = correctyab(c1,yo,yp) !Correct overshoot at cp
 yder = 2.d0*(c2-yo)/h1!estimate new derivative
 c0 = yo - 0.5d0*h0*yder!calculate cm such that spline is smooth at yo
 cond0 = interv(c0,ym,yo)!Check for overshoots at this new cm
 if(.not.cond0) then
  cm = correctyab(c0,ym,yo)
 else
  cm = c0
 endif
endif
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END SUBROUTINE Beziercm
!==============================================================================


!==============================================================================
LOGICAL FUNCTION interv(y,a,b)
!------------------------------------------------------------------------------
DOUBLE PRECISION,INTENT(IN) :: y,a,b
!------------------------------------------------------------------------------
if((a.le.b .and. y.ge.a .and. y.le.b) .or.&
  (a.gt.b .and. y.le.a .and. y.ge.b)) then
 interv = .true.
else
 interv = .false.
endif
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END FUNCTION interv
!==============================================================================

!==============================================================================
DOUBLE PRECISION FUNCTION correctyab(y,a,b)
!------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: y,a,b
DOUBLE PRECISION :: maxi,mini
!------------------------------------------------------------------------------
if(a.le.b) then
 mini = a
 maxi = b
else
 mini = b
 maxi = a
endif
!******************************************************************************
if(y.ge.mini .and. y .le. maxi) then
 correctyab = y
elseif(y.lt.mini) then
 correctyab = mini
else
 correctyab = maxi
endif
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END FUNCTION
!==============================================================================

!===========================================================================
SUBROUTINE PARAGUASB(DTM,DTP,CM,CO,CP)
!----------------------------------------------------------------------
! Subroutine for the calculation of the integral which appears in the
! formal solution of the transport equation between point M (up-wind) and
! point O, through the Short Characteristic Method by parabolic  
! interopolation of the source function
! Modified on 9/9/14 (Ernest Alsina) to return Psi_M, Psi_0 and Psi_P,
! from the difference of optical depths with the upwind and downwind 
! points
!----------------------------------------------------------------------
! INPUT
! DTM: Optical depth diffence between point 0 and point M
! DTP: Optical depth diffence between point 0 and point M
!--
! OUTPUT
! PARA: Psi_M*SM + Psi_O*SO + Psi_P*SP 
! CM: Psi_M
! CO: Psi_O
! CP: Psi_P
!----------------------------------------------------------------------
! Original subroutine written by Javier Trujillo Bueno
!----------------------------------------------------------------------
  REAL*8, INTENT(IN) :: DTM, DTP
  REAL*8, INTENT(OUT) :: CM, CO, CP
  REAL*8 :: EXU, U0, U1, U2, D2, D3, D4
  REAL*8, PARAMETER :: trick=1.d-3
!----------------------------------------------------------------------
! If the exponetial's argument's value is smaller than trick, 
! the exponenetial is calculated by its Taylor series (order 2).
!----------------------------------------------------------------------
  IF(DTM.GE.trick) THEN
    EXU=DEXP(-DTM)
    U0=1.d0-EXU
    U1=DTM-U0
    U2=DTM*DTM-2.d0*U1
  ELSE
    D2=DTM*DTM
    D3=DTM*D2
    D4=DTM*D3
    U0=DTM-(D2/2.D0)
    U1=(D2/2.D0)-(D3/6.D0)
    U2=(D3/3.D0)-(D4/12.D0)
  END IF
!----------------------------------------------------------------------
! Calcolo dei coefficienti Psi_M, Psi_O e Psi_P
!----------------------------------------------------------------------
  CM=(U2-U1*(DTP+2.d0*DTM))/(DTM*(DTM+DTP))+U0
  CO=(U1*(DTM+DTP)-U2)/(DTM*DTP)
  CP=(U2-DTM*U1)/(DTP*(DTM+DTP))

return
!----------------------------------------------------------------------
END SUBROUTINE PARAGUASB
!===========================================================================


!===========================================================================
SUBROUTINE MatVec(a,b,c)
DOUBLE PRECISION, DIMENSION(0:3,0:3), INTENT(IN) :: a
DOUBLE PRECISION, DIMENSION(0:3), INTENT(IN) :: b
DOUBLE PRECISION, DIMENSION(0:3), INTENT(OUT) :: c
DOUBLE PRECISION, DIMENSION(0:3) :: c2
INTEGER :: Nx1, Ny1, ix, iin
DOUBLE PRECISION, PARAMETER :: snap = 1.d-35

!Nx1 = size(a,1)
!Ny1 = size(a,2)
!ix = size(b,1)
!if(Ny1 .ne. ix) then
!print*, 'Elements of vector not equal to columns of matrix'
!stop
!else
!endif

!allocate(c(Nx1))
c2 = 0.d0

do ix = 0,3
 do iin = 0, 3
 c2(ix) = c2(ix) + a(ix,iin)*b(iin)
 enddo
 if(DABS(c2(ix)).le.snap) c2(ix)=0.d0
 c(ix) = c2(ix)
enddo
!----------------------------------------------------------------------
return
!----------------------------------------------------------------------
END SUBROUTINE MatVec
!===========================================================================

!===========================================================================
SUBROUTINE MatMat(a,b,c) 
!Product of two matrices
DOUBLE PRECISION, DIMENSION(0:3,0:3), INTENT(IN) :: a
DOUBLE PRECISION, DIMENSION(0:3,0:3), INTENT(IN) :: b
DOUBLE PRECISION, DIMENSION(0:3,0:3), INTENT(OUT) :: c
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: c2
INTEGER :: Nx1, Nx2, Ny1, Ny2, iin, ix, iy
DOUBLE PRECISION, PARAMETER :: snap = 1.d-35

!Nx1 = size(a,1)
!Ny1 = size(a,2)
!Nx2 = size(b,1)
!Ny2 = size(b,2)
!if(Ny1 .ne. Nx2) then
!print*, 'These matrices cannot be multiplied'
!stop
!else 
!endif

!allocate(c(Nx1,Ny2))
c2 = 0.d0

do ix = 0, 3
 do iy = 0, 3
 do iin = 0, 3
 c2(ix,iy) = c2(ix,iy) + a(ix,iin)*b(iin,iy)
 enddo!iin
if(DABS(c2(ix,iy)).le.snap) c2(ix,iy) = 0.d0
 c(ix,iy) = c2(ix,iy)
 enddo
enddo
!----------------------------------------------------------------------
return
!----------------------------------------------------------------------
END SUBROUTINE MatMat
!===========================================================================

!===========================================================================
SUBROUTINE rt_unpol(iter,Bpb,etaI,r,invk,mmat,sline,scont)
! Formal solution routine using the DELOPAR method.
INTEGER :: iter
DOUBLE PRECISION, INTENT(IN) :: Bpb(nwl)
DOUBLE PRECISION, INTENT(IN) ::  etaI(nz,nwl,ndir,ndirx),r(nz,nwl,ndir,ndirx)!Keep the r for a test. 
                      !Later we can take it out again. 
!DOUBLE PRECISION, INTENT(IN) :: etaI_c(nz,nwl)!Continuum part doesn't depend on angle
DOUBLE PRECISION, INTENT(IN) :: invk(0:3,0:3,nz,nwl,ndir,ndirx) ![1 + Phi_0*K']^-1
DOUBLE PRECISION, INTENT(IN) :: mmat(0:3,0:3,nz,nwl,ndir,ndirx) ![1 + Phi_0*K']^-1 [1e^-dt - Phu_M*K']
DOUBLE PRECISION, INTENT(IN) :: sline(2,0:2,0:2,nz,nwl), scont(2,0:2,0:2,nz,nwl)
!DOUBLE PRECISION, INTENT(INOUT) :: jrad(2,0:2,0:2,nz,nwl)
!radiation tensors: Don't consider -Q terms since J^K_-Q = (-1)^Q J^K_Q^ast
INTEGER :: k, k1, kO, kdel, i, iz, im, ax, kU, kD
INTEGER :: ik, iq, in
DOUBLE PRECISION :: dlU, dtU, dlD, dtD, exu
DOUBLE PRECISION :: ClM, ClO, CM, CO, CP
DOUBLE PRECISION :: a,b,c, eM, eO, eP!, elM, elO, elP
DOUBLE PRECISION :: mi, mi2
DOUBLE PRECISION :: cre,cim,lre,lim
DOUBLE PRECISION :: ctre, wgt, ar, ai
DOUBLE PRECISION, DIMENSION(0:3) :: SM, SO, SP
DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
DOUBLE PRECISION, DIMENSION(0:3) :: stokes, sto_prev, svector
DOUBLE PRECISION, DIMENSION(0:3) :: avec, bvec !not an allocatable
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: amat, bmat
DOUBLE PRECISION, PARAMETER :: trick = 1.d-3
CHARACTER(len=1024) :: filop,format_op
LOGICAL :: flov
!Initializing the Stokes parameters as zero:
do in=0,3!TEST
 stokes(in) = 0.d0
 sto_prev(in) = 0.d0
 svector(in) = 0.d0
enddo!in
! Set the radiation tensor initially to zero
jrout = 0.d0
! Initialize constant
ctre = 0.d0
df = 0!disaster flag
if(iter.eq.1) open(unit=29, file='rt_test.dat')
!TEST
!if(MOD(iter,50).eq.1) then
! if(iter .lt. 10) then
!  format_op = "(A5,I1)"
! elseif(iter.lt.100) then
!  format_op = "(A5,I2)"
! else
!  format_op = "(A5,I3)"
! endif
!write(filop,format_op) "itef_",iter
! filop = trim(filop)//".res"
! open(unit=31,file=trim(filop))
!endif
!---------------------------------------------------------------
! Begin frequency loop
!---------------------------------------------------------------
 do i = 1, nwl
!---------------------------------------------------------------
! Start loops over polar angle (using Gauss-Legendre quadrature)
! and azimuthal angle (trapezoidal quadrature). The weights are
! not considered here because we aren't performing a numerical
! integration here. For the polar angle, the mu values are found
! so they go from -1 to 1
!---------------------------------------------------------------
  do im = 1, ndir !polar angle goes from 0 to pi
   mi = mu(im)
   mi2 = mi*mi
   do ax = 1, ndirx ! goes from 0 to 2*pi
   !alph = az(ax)
    wgt = wtdir(im)*waz(ax)/(4.d0*pi)
!---------------------------------------------------------------
! Start loop over heights. Range will depend on the polar angle;
! i.e. whether we are going inward or outward
!---------------------------------------------------------------
    if (mi .le. 0) then !INWARD
     kO = 1
     k1 = nz
     kdel = 1
    else              !OUTWARD
     kO = nz
     k1 = 1
     kdel = -1
    endif
    !First height point
    do k = kO + kdel, k1, kdel
    ! UPWIND
     kU = k-kdel
     dlU = DABS((ht(kU) - ht(k))/mi)
     dtU = 0.5d0*(etaI(k,i,im,ax) + etaI(kU,i,im,ax))*dlU !Change in optical depth
    ! DOWNWIND
     if(k.ne.k1) then
      kD = k+kdel
      dlD = DABS((ht(k) - ht(kD))/mi)
      dtD = 0.5d0*(etaI(k,i,im,ax) + etaI(kD,i,im,ax))*dlD !Change in optical depth
     else !If we are on the surface kD is undefined (out of the array)
      dlD = dlU
      dtD = dtU
     endif
    ! Calculate e^-{delta tau_M}
     if(dtU.ge.trick) then
      exu = DEXP(-dtU)
     else
      exu =  1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
       (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
       (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
       (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0
     endif

     IF(k.ne.k1) THEN
     !ecP = etaI_c(k+kdel,i)
      eP = etaI(k+kdel,i,im,ax)
      IF(k.eq.kO+kdel) THEN !First point
       ! ecO = etaI_c(k,i)
       !ecM = etaI_c(k-kdel,i)
       eO = etaI(k,i,im,ax)
       eM = etaI(k-kdel,i,im,ax)
!       do in = 0,3
        in = 0
        a = 0.d0; b = 0.d0; c = 0.d0
        do ik = 0, 2
         do iq = 0, ik
!------------------------------------------------------------------------------------------------
! Calculating total S_i source functions at upwind (M), current (O) and downwind (P) points. 
! It is insider a loop over Stokes vector components 0 through 3
!------------------------------------------------------------------------------------------------
         if(iq.eq.0) then
          ctre = 1.d0 !ONLY REAL PART IF Q=0
         else
   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
          ctre = 2.d0
         endif
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
         a = a + ctre*(lre + cre)/eP!P point is downwind
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k,i),sline(2,ik,iq,k,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k,i),scont(2,ik,iq,k,i),cre,cim)
         b =  b + ctre*(lre + cre)/eO!current point
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         sline(1,ik,iq,k-kdel,i),sline(2,ik,iq,k-kdel,i),lre,lim)
         call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
         scont(1,ik,iq,k-kdel,i),scont(2,ik,iq,k-kdel,i),cre,cim)
         c = c + ctre*(lre+cre)/eM !M is the upwind point
         enddo !iq
        enddo!ik
        SP(in) = a
        SO(in) = b
        SM(in) = c
!      enddo!in
      do in = 1,3
       SP(in) = a
       SO(in) = b
       SM(in) = c
      enddo!in
      sto_prev = 0.d0 !Boundary point
      if(mi.ge.0) sto_prev(0) = Bpb(i)!Bottom point. Assumed thermalized.
       do ik = 0,2
        do iq = 0,ik!At boundary conditions we assume no (incoming) polarization, neither from the bottom nor top. Therefore
         jrout(1,ik,iq,kO,i) = jrout(1,ik,iq,kO,i) + wgt*TTR(ik,iq,0,im,ax)*sto_prev(0)!we only need 0 Stokes index
         if(iq.gt.0) jrout(2,ik,iq,kO,i) = jrout(2,ik,iq,kO,i) + wgt*TTI(ik,iq,0,im,ax)*sto_prev(0)
        enddo ! iq
       enddo ! ik
      ELSE ! Middle points
       !HERE we only need to look for SP, since SO and SM have already been found in previous k points
!       do in = 0, 3
        in = 0
        a = 0.d0!; b = 0.d0; c = 0.d0
        do ik = 0, 2
         do iq = 0, ik
          if(iq.eq.0) then
           ctre = 1.d0 !ONLY REAL PART IF Q=0
          else
! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
           ctre = 2.d0
          endif
          call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
          sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
          call comprod(TTR(ik,iq,in,im,ax),TTI(ik,iq,in,im,ax),&
          scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
          a = a + ctre*(lre+cre)/eP
         enddo !iq
        enddo!ik
        SP(in) = a
!       enddo!in
       do in = 1,3
        SP(in) = 0.d0
       enddo
      ENDIF
!----------------------------------------------------------------------------------------------
! BESSER routine itself
!----------------------------------------------------------------------------------------------
      call PARAGUASB(dtU,dtD,CM,CO,CP)!For tests
      call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
      in = 0
!      do in = 0,3
       call Beziercm(dtU,dtD,SM(in),SO(in),SP(in),cpm,flov)
!       svector(in) = CM*SM(in) + CO*SO(in) + CP*SP(in)
       svector(in) = omega_m*SM(in) + omega_o*SO(in) + omega_c*cpm
!      enddo!in
      do in = 1,3
       svector(in) = 0.d0
      enddo
      amat(:,:) = mmat(:,:,k,i,im,ax)
      bmat(:,:) = invk(:,:,k,i,im,ax)
      call MatVec(amat,sto_prev,avec) ![1 + lPhi_O K'_O]^-1 [e^-dtM -  lPhi_M K'_M] I_M
      call MatVec(bmat,svector,bvec) ![1 + lPhi_O K'_O]^-1 (Phi_M S_M + Phi_O S_O + Phi_P S_P)
      in = 0
!      do in = 0,3!TEST
       stokes(in) = avec(in) + bvec(in)
       sto_prev(in) = stokes(in)
       SM(in) = SO(in)
       SO(in) = SP(in)
       !ecM = ecO
       !ecO = ecP
      ! eM = eO
      ! eO = eP
!      enddo !in
      do in = 1,3
       SM(in) = 0.d0
       SO(in) = 0.d0
      enddo
      eM = eO
      eO = eP
     ELSE !k=k1. LINEAR INTERPOLATION INSTEAD OF PARABOLIC
      call LINEAR(exu,dtU,ClO,ClM)
       in = 0
!      do in=0,3
       svector(in) = ClM*SM(in) + ClO*SO(in)
!      enddo!in
      do in = 1,3
       svector(in) = 0.d0
      enddo
      amat(:,:) = mmat(:,:,k,i,im,ax)
      bmat(:,:) = invk(:,:,k,i,im,ax)
      call MatVec(amat,sto_prev,avec)
      call MatVec(bmat,svector,bvec)
      do in = 0,3
       stokes(in) = avec(in) + bvec(in)
      enddo
     ENDIF
!----------------------------------------------------------------------
! Control
!----------------------------------------------------------------------
     IF(stokes(0).lt.0.d0) THEN
       print*,'DISASTER! I .lt. 0 at k=', k, stokes(0), sline(1,0,0,k,i), sline(1,1,0,k,i) ,sline(1,2,0,k,i)
      df = 1!flag that indicates something going wrong in the rt_solver
     ELSE
     ENDIF
!----------------------------------------------------------------------
! Calculate the radiation tensors (multipolar components)
!----------------------------------------------------------------------
     do ik = 0,2!
      do iq = 0,ik
       ar = 0.d0
       ai = 0.d0
       do in = 0,3 !summing over the stokes parameters
        ar = ar + TTR(ik,iq,in,im,ax)*stokes(in)
        if(iq.gt.0) ai = ai + TTI(ik,iq,in,im,ax)*stokes(in)
       enddo!in
       jrout(1,ik,iq,k,i) = jrout(1,ik,iq,k,i) + wgt*ar
       if(iq.gt.0) jrout(2,ik,iq,k,i) = jrout(2,ik,iq,k,i) + wgt*ai
      enddo!iq
     enddo!ik
!---RADIATION TENSOR FOUND
    enddo!k End height grid
   enddo!ax End azimuthal angle grid
  enddo!im End polar angle grid
 enddo!i End frequency grid
close(unit=19)
!TEST
!close(unit=31)
!TEST
!----------------------------------------------------------------------
 return
!----------------------------------------------------------------------
END SUBROUTINE rt_unpol
!===========================================================================

!===========================================================================
SUBROUTINE rt_emer_par(mi,chi,TR,TI,Bpb,etaI_c,sline,scont,eta_ref,rho_ref,emer)
DOUBLE PRECISION, INTENT(IN) :: mi, chi
DOUBLE PRECISION, INTENT(IN) :: TR(0:2,-2:2,0:3), TI(0:2,-2:2,0:3)
DOUBLE PRECISION, INTENT(IN) :: Bpb(nwl)
!DOUBLE PRECISION, INTENT(IN) :: sr(nz,nwl), setaI(nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: etaI_c(nz,nwl)!Continuum part doesn't depend on angle
!DOUBLE PRECISION, INTENT(IN) :: sinvk(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1
!DOUBLE PRECISION, INTENT(IN) :: smmat(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1 [1e^-dt - Phu_M*K']
DOUBLE PRECISION, INTENT(IN) :: eta_ref(0:3,nz,nwl), rho_ref(0:3,nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: sline(2,0:2,0:2,nz,nwl), scont(2,0:2,0:2,nz,nwl)
DOUBLE PRECISION, INTENT(OUT) :: emer(0:3,nwl)
INTEGER :: k, k1, kO, kdel, i, iz, im, ax, kU, kD
INTEGER :: ik, iq, in
DOUBLE PRECISION :: dlU, dtU, dlD, dtD, exu
DOUBLE PRECISION :: ClM, ClO, CM, CO, CP
DOUBLE PRECISION :: a,b,c, ecM, ecO, ecP, elM, elO, elP
DOUBLE PRECISION :: cre, cim, lre, lim
DOUBLE PRECISION :: ctre, wgt, ar, ai
DOUBLE PRECISION :: rP, rM, rO
DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
DOUBLE PRECISION, DIMENSION(0:3) :: SM, SO, SP,cptst
DOUBLE PRECISION, DIMENSION(0:3) :: stokes, sto_prev, svector
DOUBLE PRECISION, DIMENSION(0:3) :: avec, bvec! not allocatable anymore
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: bmat !not allocatable anymore
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: amat, omat, tmat, smat
DOUBLE PRECISION, PARAMETER :: trick = 1.d-3
LOGICAL :: flov
!Initializing the Stokes parameters as zero:
!if(mi.eq.0.16d0) open(unit=31,file='itef_last.res')
stokes = 0.d0
sto_prev = 0.d0
svector = 0.d0
emer = 0.d0
! Initialize constant
ctre = 0.d0
df = 0!disaster flag
!---------------------------------------------------------------
! Begin frequency loop
!---------------------------------------------------------------
do i = 1, nwl
!---------------------------------------------------------------
! Start loops over polar angle (using Gauss-Legendre quadrature)
! and azimuthal angle (trapezoidal quadrature). The weights are
! not considered here because we aren't performing a numerical
! integration here. For the polar angle, the mu values are found
! so they go from -1 to 1
!---------------------------------------------------------------
if (mi .le. 0) then !INWARD
print*, 'We should be considering an emergent ray'
stop
kO = 1
k1 = nz
kdel = 1
else              !OUTWARD
kO = nz
k1 = 1
kdel = -1
endif
!First height point
do k = kO + kdel, k1, kdel
! UPWIND
kU = k-kdel
dlU = DABS((ht(kU) - ht(k))/mi)
dtU = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kU,i) + eta_ref(0,kU,i))*dlU !Change in optical depth
! DOWNWIND
if(k.ne.k1) then
kD = k+kdel
dlD = DABS((ht(k) - ht(kD))/mi)
dtD = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kD,i) + eta_ref(0,kD,i))*dlD !Change in optical depth
else !If we are on the surface kD is undefined (out of the array)
dlD = dlU
dtD = dtU
endif
! Calculate e^-{delta tau_M}
if(dtU.ge.trick) then
exu = DEXP(-dtU)
else
exu = 1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
  (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
  (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0

endif
!------------------------------------------------------------------------------
! Calculating the K_0 matrices and K_M at the first point
!------------------------------------------------------------------------------
if(k.eq.kO+kdel) then
smat = 0.d0
smat(0,1) = eta_ref(1,kO,i)! Line contribution is the only one there is
smat(0,2) = eta_ref(2,kO,i)! since continuum part is zero 
smat(0,3) = eta_ref(3,kO,i)
smat(1,0) = eta_ref(1,kO,i)
smat(1,2) = rho_ref(3,kO,i)
smat(1,3) = -rho_ref(2,kO,i)
smat(2,0) = eta_ref(2,kO,i)
smat(2,1) = -rho_ref(3,kO,i)
smat(2,3) = rho_ref(1,kO,i)
smat(3,0) = eta_ref(3,kO,i)
smat(3,1) = rho_ref(2,kO,i)
smat(3,2) = -rho_ref(1,kO,i)
smat = smat/(etaI_c(kO,i) + eta_ref(0,kO,i))
else
smat = tmat!K_M
endif
tmat = 0.d0
tmat(0,1) = eta_ref(1,k,i)! Line contribution is the only one there is
tmat(0,2) = eta_ref(2,k,i)! since continuum part is zero 
tmat(0,3) = eta_ref(3,k,i)
tmat(1,0) = eta_ref(1,k,i)
tmat(1,2) = rho_ref(3,k,i)
tmat(1,3) = -rho_ref(2,k,i)
tmat(2,0) = eta_ref(2,k,i)
tmat(2,1) = -rho_ref(3,k,i)
tmat(2,3) = rho_ref(1,k,i)
tmat(3,0) = eta_ref(3,k,i)
tmat(3,1) = rho_ref(2,k,i)
tmat(3,2) = -rho_ref(1,k,i)
tmat = tmat/(etaI_c(k,i)+eta_ref(0,k,i))!K_O
call LINEAR(exu,dtU,ClO,ClM)
bmat = iden + ClO*tmat
call MatInv(bmat)
omat = exu*iden - ClM*smat
call MatMat(bmat,omat,amat)

IF(k.ne.k1) THEN
ecP = etaI_c(k+kdel,i)
elP = eta_ref(0,k+kdel,i)! - ecP
rP = elP/(ecP+elP)
  IF(k.eq.kO+kdel) THEN !First point
  ecO = etaI_c(k,i)
  ecM = etaI_c(k-kdel,i)
  elO = eta_ref(0,k,i)
  elM = eta_ref(0,k-kdel,i)
  rO = elO/(ecO + elO)
  rM = elM/(ecM + elM)
  do in = 0,3
  a = 0.d0; b = 0.d0; c = 0.d0
  do ik = 0, 2
  do iq = 0, ik
  if(iq.eq.0) then
   ctre = 1.d0 !ONLY REAL PART IF Q=0
  else
   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART, due to +Q and -Q terms
   ctre = 2.d0
  endif
  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
   a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k,i),sline(2,ik,iq,k,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k,i),scont(2,ik,iq,k,i),cre,cim)
   b =  b + ctre*rO*lre/elO + ctre*(1.d0 - rO)*cre/ecO
  call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k-kdel,i),sline(2,ik,iq,k-kdel,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k-kdel,i),scont(2,ik,iq,k-kdel,i),cre,cim)
   c = c + ctre*rM*lre/elM + ctre*(1.d0 - rM)*cre/ecM
  enddo !iq
  enddo!ik
  SP(in) = a
  SO(in) = b
  SM(in) = c
  enddo!in
  sto_prev = 0.d0 !Boundary point
  sto_prev(0) = Bpb(i)!Bottom point. Assumed thermalized.
  ELSE ! Middle points
   !HERE we only need to look for SP, since SO and SM have already been found in previous k points
  do in = 0,3
  a = 0.d0; b = 0.d0; c = 0.d0
  do ik = 0, 2
  do iq = 0, ik
     if(iq.eq.0) then
     ctre = 1.d0 !ONLY REAL PART IF Q=0
     else
     ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
     ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
     ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
     ctre = 2.d0
     endif
  call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
  call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
   a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
  enddo !iq
  enddo!ik
  SP(in) = a
  enddo!in
  ENDIF

call PARAGUASB(dtU,dtD,CM,CO,CP)
call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
  svector = CM*SM + CO*SO + CP*SP
! do in = 0,3
!  call Beziercm(dtU,dtD,SM(in),SO(in),SP(in),cpm,flov)
!!  cptst(in) = cpm
!  svector(in) = omega_m*SM(in) + omega_o*SO(in) + omega_c*cpm
!TEST
!  if(mi.eq.0.16d0) then
!   if(flov) then
!    write(31,*) in,k,i,SM(in),SO(in),SP(in)
!    print*,in,k,i,SM(in),SO(in),SP(in)
!   endif
!  endif
!TEST
!!  svector(in) = CM*SM(in) + CO*SO(in) + CP*SP(in)
! enddo
 ! write(34,*) k,i,CM,CO,CP,omega_m,omega_o,omega_c,SM(0),SO(0),SP(0),cptst(0),&
 !            SM(1),SO(1),SP(1),cptst(1)
  call MatVec(amat,sto_prev,avec)
  call MatVec(bmat,svector,bvec)
 stokes = avec + bvec
!  deallocate(avec)
!  deallocate(bvec)
 sto_prev = stokes
!-------------TEST--------------------------------------
SM = SO
SO = SP
ecM = ecO
ecO = ecP
elM = elO
elO = elP
rM = rO
rO = rP
ELSE !k=k1. LINEAR INTERPOLATION INSTEAD OF PARABOLIC
!We already have SM and SO, and we don't need SP
!call PLINEA(dtU,ClM,ClO) Already computed in the loop
 svector = ClM*SM + ClO*SO
!  amat(:,:) = smmat(:,:,k,i)
!  bmat(:,:) = sinvk(:,:,k,i)
call MatVec(amat,sto_prev,avec)
call MatVec(bmat,svector,bvec)
 stokes = avec + bvec

 do in = 0,3
 emer(in,i) = stokes(in)
 enddo!in
ENDIF
!----------------------------------------------------------------------
! Control
!----------------------------------------------------------------------
!IF(stokes(0).lt.0.0) THEN
! print*,'DISASTER! I .lt. 0 at k=', k, sline(1,0,0,k,i), sline(1,2,0,k,i)
! df = 1!flag that indicates something going wrong in the rt_solver
!ELSE
!ENDIF
enddo!k grid
enddo !i
!TEST
!if(mi.eq.31) close(unit=31)
!----------------------------------------------------------------------
 return
!----------------------------------------------------------------------
END SUBROUTINE rt_emer_par
!===========================================================================



!!==============================================================================
!FUNCTION PARABO(dtM,dtP,SO,SM,SP)
!!----------------------------------------------------------------------
!! Funzione pari al valore dell'integrale che compare nella soluzione
!! formale dell'equazione del trasporto tra il punto M (up-wind) e il 
!! punto O, quale si ottiene applicando lo Short Characteristic Method 
!! con interpolazione parabolica della funzione sorgente.
!!----------------------------------------------------------------------
!! INPUT
!! dtM: differenza di profondita' ottica tra il punto 0 e il punto M
!! dtP: differenza di profondita' ottica tra il punto 0 e il punto P
!! SO: funzione sorgente nel punto O
!! SM: funzione sorgente nel punto M
!! SP: funzione sorgente nel punto P
!!--
!! OUTPUT
!! PARABO: Psi_M*SM + Psi_O*SO + Psi_P*SP
!!----------------------------------------------------------------------
!! Subroutine originale scritta da Javier Trujillo Bueno
!!----------------------------------------------------------------------
!  DOUBLE PRECISION :: PARABO
!  DOUBLE PRECISION, INTENT(IN) :: dtM, dtP, SO, SM, SP
!  DOUBLE PRECISION :: exu, U0, U1, U2, D2, D3, D4, CM, CO, CP
!  DOUBLE PRECISION, PARAMETER :: trick=1.d-3
!!----------------------------------------------------------------------
!  if(dtM.ge.trick) then
!    exu=DEXP(-dtM)
!    U0=1.d0-exu
!    U1=dtM-U0
!    U2=dtM*dtM-2.d0*U1
!  else
!    D2=dtM*dtM
!    D3=dtM*D2
!    D4=dtM*D3
!    U0=dtM-(D2/2.D0)
!    U1=(D2/2.D0)-(D3/6.D0)
!    U2=(D3/3.D0)-(D4/12.D0)
!  endif
!!----------------------------------------------------------------------
!  CM=(U2-U1*(dtP+2.d0*dtM))/(dtM*(dtM+dtP))+U0
!  CO=(U1*(dtM+dtP)-U2)/(dtM*dtP)
!  CP=(U2-dtM*U1)/(dtP*(dtM+dtP))
!!----------------------------------------------------------------------
!  PARABO=CM*SM+CO*SO+CP*SP
!!----------------------------------------------------------------------
!END FUNCTION PARABO

!FUNCTION PLINEA(dtM,SO,SM)
!!----------------------------------------------------------------------
!! Funzione pari al valore dell'integrale che compare nella soluzione
!! formale dell'equazione del trasporto tra il punto M (up-wind) e il 
!! punto O, quale si ottiene applicando lo Short Characteristic Method 
!! con interpolazione lineare della funzione sorgente.
!!----------------------------------------------------------------------
!! INPUT
!! dtM: differenza di profondita' ottica tra il punto 0 e il punto M
!! SO: funzione sorgente nel punto O
!! SM: funzione sorgente nel punto M
!!--
!! OUTPUT
!! PLINEA: Psi_M*SM + Psi_O*SO
!!----------------------------------------------------------------------
!! Subroutine originale scritta da Javier Trujillo Bueno
!!----------------------------------------------------------------------
!  DOUBLE PRECISION :: PLINEA
!  DOUBLE PRECISION, INTENT(IN) :: dtM, SO, SM
!  DOUBLE PRECISION :: exu, U0, U1, D2, CO, CM
!  DOUBLE PRECISION, PARAMETER :: trick=1.d-3
!!----------------------------------------------------------------------
!  if(dtM.ge.trick) then
!    exu=DEXP(-dtM)
!    U0=1.d0-exu
!    U1=dtM-U0
!    CO=U1/dtM
!    CM=U0-CO
!  else
!    D2=dtM*dtM
!    CO=(dtM/2.D0)-(D2/6.D0)
!    CM=(dtM/2.D0)-(D2/3.D0)
!  endif
!!----------------------------------------------------------------------
! PLINEA=CM*SM+CO*SO
!return
!!----------------------------------------------------------------------
!END FUNCTION PLINEA




!SUBROUTINE rt_emer_extra(mi,chi,TR,TI,Bpb,etaI_c,sline,scont,eta_ref,&
!  rho_ref,emer)
!DOUBLE PRECISION, INTENT(IN) :: mi, chi
!DOUBLE PRECISION, INTENT(IN) :: TR(0:2,-2:2,0:3), TI(0:2,-2:2,0:3)
!DOUBLE PRECISION, INTENT(IN) :: Bpb(nwl)
!!DOUBLE PRECISION, INTENT(IN) :: sr(nz,nwl), setaI(nz,nwl)
!DOUBLE PRECISION, INTENT(IN) :: etaI_c(nz,nwl)!Continuum part doesn't depend on angle
!!DOUBLE PRECISION, INTENT(IN) :: sinvk(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1
!!DOUBLE PRECISION, INTENT(IN) :: smmat(0:3,0:3,nz,nwl) ![1 + Phi_0*K']^-1 [1e^-dt - Phu_M*K']
!DOUBLE PRECISION, INTENT(IN) :: eta_ref(0:3,nz,nwl), rho_ref(0:3,nz,nwl)
!DOUBLE PRECISION, INTENT(IN) :: sline(2,0:2,0:2,nz,nwl), scont(2,0:2,0:2,nz,nwl)
!DOUBLE PRECISION, INTENT(OUT) :: emer(0:3,nwl)
!INTEGER :: k, k1, kO, kdel, i, iz, im, ax, kU, kD
!INTEGER :: ik, iq, in
!DOUBLE PRECISION :: dlU, dtU, dlD, dtD, exu
!DOUBLE PRECISION :: ClM, ClO, CM, CO, CP
!DOUBLE PRECISION :: a,b,c, ecM, ecO, ecP, elM, elO, elP
!DOUBLE PRECISION :: cre, cim, lre, lim
!DOUBLE PRECISION :: ctre, wgt, ar, ai
!DOUBLE PRECISION :: rP, rM, rO
!DOUBLE PRECISION :: omega_m,omega_o,omega_c,cpm
!DOUBLE PRECISION, DIMENSION(0:3) :: SM, SO, SP,cptst
!DOUBLE PRECISION, DIMENSION(0:3) :: stokes, sto_prev, svector
!DOUBLE PRECISION, DIMENSION(0:3) :: avec, bvec! not allocatable anymore
!DOUBLE PRECISION, DIMENSION(0:3,0:3) :: bmat !not allocatable anymore
!DOUBLE PRECISION, DIMENSION(0:3,0:3) :: amat, omat, tmat, smat
!DOUBLE PRECISION, PARAMETER :: trick = 1.d-3
!LOGICAL :: flov
!!Initializing the Stokes parameters as zero:
!!if(mi.eq.0.16d0) open(unit=31,file='itef_last.res')
!stokes = 0.d0
!sto_prev = 0.d0
!svector = 0.d0
!emer = 0.d0
!! Initialize constant
!ctre = 0.d0
!df = 0!disaster flag
!!---------------------------------------------------------------
!! Begin frequency loop
!!---------------------------------------------------------------
!do i = 1, nwl
!!---------------------------------------------------------------
!! Start loops over polar angle (using Gauss-Legendre quadrature)
!! and azimuthal angle (trapezoidal quadrature). The weights are
!! not considered here because we aren't performing a numerical
!! integration here. For the polar angle, the mu values are found
!! so they go from -1 to 1
!!---------------------------------------------------------------
!if (mi .le. 0) then !INWARD
!print*, 'We should be considering an emergent ray'
!stop
!kO = 1
!k1 = nz
!kdel = 1
!else     !OUTWARD
!kO = nz
!k1 = 1
!kdel = -1
!endif
!!First height point
!do k = kO + kdel, k1, kdel
!! UPWIND
!kU = k-kdel
!dlU = DABS((ht(kU) - ht(k))/mi)
!dtU = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kU,i) + eta_ref(0,kU,i))*dlU !Change in optical depth
!! DOWNWIND
!if(k.ne.k1) then
!kD = k+kdel
!dlD = DABS((ht(k) - ht(kD))/mi)
!dtD = 0.5d0*(eta_ref(0,k,i) + etaI_c(k,i) + etaI_c(kD,i) + eta_ref(0,kD,i))*dlD !Change in optical depth
!else !If we are on the surface kD is undefined (out of the array)
!dlD = dlU
!dtD = dtU
!endif
!! Calculate e^-{delta tau_M}
!if(dtU.ge.trick) then
!exu = DEXP(-dtU)
!else
!exu = 1.d0 - dtU + (dtU*dtU)/2.d0 - (dtU*dtU*dtU)/6.d0 + & 
!  (dtU*dtU*dtU*dtU)/24.d0 - (dtU*dtU*dtU*dtU*dtU)/120.d0 + &
!  (dtU*dtU*dtU*dtU*dtU*dtU)/720.d0 - (dtU*dtU*dtU*dtU*dtU*dtU*dtU)/5040 + &
!  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/40320.d0 - &
!  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/362880.d0 + &
!  (dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU*dtU)/3628800.d0
!endif
!!------------------------------------------------------------------------------
!! Calculating the K_0 matrices and K_M at the first point
!!------------------------------------------------------------------------------
!if(k.eq.kO+kdel) then
!smat = 0.d0
!smat(0,1) = eta_ref(1,kO,i)! Line contribution is the only one there is
!smat(0,2) = eta_ref(2,kO,i)! since continuum part is zero 
!smat(0,3) = eta_ref(3,kO,i)
!smat(1,0) = eta_ref(1,kO,i)
!smat(1,2) = rho_ref(3,kO,i)
!smat(1,3) = -rho_ref(2,kO,i)
!smat(2,0) = eta_ref(2,kO,i)
!smat(2,1) = -rho_ref(3,kO,i)
!smat(2,3) = rho_ref(1,kO,i)
!smat(3,0) = eta_ref(3,kO,i)
!smat(3,1) = rho_ref(2,kO,i)
!smat(3,2) = -rho_ref(1,kO,i)
!smat = smat/(etaI_c(kO,i) + eta_ref(0,kO,i))
!else
!smat = tmat!K_M
!endif
!tmat = 0.d0
!tmat(0,1) = eta_ref(1,k,i)! Line contribution is the only one there is
!tmat(0,2) = eta_ref(2,k,i)! since continuum part is zero 
!tmat(0,3) = eta_ref(3,k,i)
!tmat(1,0) = eta_ref(1,k,i)
!tmat(1,2) = rho_ref(3,k,i)
!tmat(1,3) = -rho_ref(2,k,i)
!tmat(2,0) = eta_ref(2,k,i)
!tmat(2,1) = -rho_ref(3,k,i)
!tmat(2,3) = rho_ref(1,k,i)
!tmat(3,0) = eta_ref(3,k,i)
!tmat(3,1) = rho_ref(2,k,i)
!tmat(3,2) = -rho_ref(1,k,i)
!tmat = tmat/(etaI_c(k,i)+eta_ref(0,k,i))!K_O
!call LINEAR(exu,dtU,ClO,ClM)
!bmat = iden + ClO*tmat
!call MatInv(bmat)
!omat = exu*iden - ClM*smat
!call MatMat(bmat,omat,amat)
!
!IF(k.ne.k1) THEN
!ecP = etaI_c(k+kdel,i)
!elP = eta_ref(0,k+kdel,i)! - ecP
!rP = elP/(ecP+elP)
!  IF(k.eq.kO+kdel) THEN !First point
!  ecO = etaI_c(k,i)
!  ecM = etaI_c(k-kdel,i)
!  elO = eta_ref(0,k,i)
!  elM = eta_ref(0,k-kdel,i)
!  rO = elO/(ecO + elO)
!  rM = elM/(ecM + elM)
!  do in = 0,3
!  a = 0.d0; b = 0.d0; c = 0.d0
!  do ik = 0, 2
!  do iq = 0, ik
!  if(iq.eq.0) then
!   ctre = 1.d0 !ONLY REAL PART IF Q=0
!  else
!   ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
!   ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
!   ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART, due to +Q and -Q terms
!   ctre = 2.d0
!  endif
!  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
!  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
!   a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
!  call comprod(TR(ik,iq,in),TI(ik,iq,in), sline(1,ik,iq,k,i),sline(2,ik,iq,k,i),lre,lim)
!  call comprod(TR(ik,iq,in),TI(ik,iq,in), scont(1,ik,iq,k,i),scont(2,ik,iq,k,i),cre,cim)
!   b =  b + ctre*rO*lre/elO + ctre*(1.d0 - rO)*cre/ecO
!  call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k-kdel,i),sline(2,ik,iq,k-kdel,i),lre,lim)
!  call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k-kdel,i),scont(2,ik,iq,k-kdel,i),cre,cim)
!   c = c + ctre*rM*lre/elM + ctre*(1.d0 - rM)*cre/ecM
!  enddo !iq
!  enddo!ik
!  SP(in) = a
!  SO(in) = b
!  SM(in) = c
!  enddo!in
!  sto_prev = 0.d0 !Boundary point
!  sto_prev(0) = Bpb(i)!Bottom point. Assumed thermalized.
!  stokh(0,kO,i) = Bpb(i)
!  do in = 1,3
!   stokh(in,kO,i) = 0.d0
!  enddo
!  ELSE ! Middle points
!   !HERE we only need to look for SP, since SO and SM have already been found in previous k points
!   do in = 0,3
!   a = 0.d0; b = 0.d0; c = 0.d0
!    do ik = 0, 2
!     do iq = 0, ik
!      if(iq.eq.0) then
!       ctre = 1.d0 !ONLY REAL PART IF Q=0
!      else
!     ! Here we are looking for |Q|>0 components. Consider positive and negative Q in the same line
!     ! WE DON'T NEED TO CONSIDER IMAG. PART SINCE IT MUST CANCEL W/ THE -Q COMPONENT
!     ! WE DO, HOWEVER, GET A FACTOR 2 FOR THE REAL PART.
!       ctre = 2.d0
!      endif
!      call comprod(TR(ik,iq,in),TI(ik,iq,in),sline(1,ik,iq,k+kdel,i),sline(2,ik,iq,k+kdel,i),lre,lim)
!      call comprod(TR(ik,iq,in),TI(ik,iq,in),scont(1,ik,iq,k+kdel,i),scont(2,ik,iq,k+kdel,i),cre,cim)
!      a = a + ctre*rP*lre/elP + ctre*(1.d0 - rP)*cre/ecP
!     enddo !iq
!    enddo!ik
!    SP(in) = a
!   enddo!in
!  ENDIF
!!call PARAGUASB(dtU,dtD,CM,CO,CP)
! call OMEGA(exu,dtU,omega_m,omega_o,omega_c)
!!  svector = CM*SM + CO*SO + CP*SP
! do in = 0,3
!  call Beziercm(dtU,dtD,SM(in),SO(in),SP(in),cpm,flov)
!  svector(in) = omega_m*SM(in) + omega_o*SO(in) + omega_c*cpm
!!  svector(in) = CM*SM(in) + CO*SO(in) + CP*SP(in)
! enddo
!  call MatVec(amat,sto_prev,avec)
!  call MatVec(bmat,svector,bvec)
! stokes = avec + bvec
! do in = 0,3
!  stokh(in,k,i) = stokes(in)
! enddo
! sto_prev = stokes
!!-------------TEST--------------------------------------
!SM = SO
!SO = SP
!ecM = ecO
!ecO = ecP
!elM = elO
!elO = elP
!rM = rO
!rO = rP
!ELSE !k=k1. LINEAR INTERPOLATION INSTEAD OF PARABOLIC
!!We already have SM and SO, and we don't need SP
!!call PLINEA(dtU,ClM,ClO) Already computed in the loop
! svector = ClM*SM + ClO*SO
!!  amat(:,:) = smmat(:,:,k,i)
!!  bmat(:,:) = sinvk(:,:,k,i)
!call MatVec(amat,sto_prev,avec)
!call MatVec(bmat,svector,bvec)
! stokes = avec + bvec
! do in = 0,3
!  emer(in,i) = stokes(in)
!  stokh(in,k,i) = stokes(in)
! enddo!in
!ENDIF
!!----------------------------------------------------------------------
!! Control
!!----------------------------------------------------------------------
!!IF(stokes(0).lt.0.0) THEN
!! print*,'DISASTER! I .lt. 0 at k=', k, sline(1,0,0,k,i), sline(1,2,0,k,i)
!! df = 1!flag that indicates something going wrong in the rt_solver
!!ELSE
!!ENDIF
!enddo!k grid
!enddo !i
!!TEST
!!if(mi.eq.31) close(unit=31)
!!----------------------------------------------------------------------
! return
!!----------------------------------------------------------------------
!END SUBROUTINE rt_emer_extra

END MODULE rt_formal

