MODULE redis
!------------------------------------------------------------------------------
! Ernest Alsina 5/02/2015
! IAC
! Finds the frequency-dependent part of the redistribution
! matrices in order to find the frequency-redistributed 
! radiation tensors. They will carry real and imaginary parts.
! The subroutines calculate K, K' and Q components
! ------------------------------------------------------------
USE global_vars
USE parameters
USE math
USE profile!
USE abs_prof!
! ------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------


!=========================================================================
SUBROUTINE crd_red(epsprm,D2,Qel,Hup,Hup2,gprof_re,gprof_im,r3)
!------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: epsprm(1:nz), D2(0:ju2,1:nz),Qel(1:nz),Hup(0:ju2,1:nz),Hup2(1:nz)
DOUBLE PRECISION, INTENT(IN) :: gprof_re(0:ju2,0:2,-2:2,1:nz,1:nwl)
DOUBLE PRECISION, INTENT(IN) :: gprof_im(0:ju2,0:2,-2:2,1:nz,1:nwl)
DOUBLE PRECISION, INTENT(OUT) :: r3(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz)
INTEGER :: ku, ik, ikp, iq, kmin, iz,i,ip
DOUBLE PRECISION :: c, d, e, f,gre,gim, t1,t2,a1,b1
DOUBLE PRECISION,DIMENSION(-2:2) :: denom2
DOUBLE PRECISION,DIMENSION(0:ju2,-2:2) :: denom1
!------------------------------------------------------------------------
call cpu_time(t1)
print*,'Enters R3'
do iz = 1, nz
 do ik = 0,ju2
  do iq = -2,2 !we'll won't need more than 2, since iq is also limited by ikp
   denom1(ik,iq) = (1.d0 + epsprm(iz) + D2(ik,iz)/Aul)*(1.d0 + (Hup(ik,iz)*dble(iq))**2)
   denom2(iq) = (1.d0 + epsprm(iz) + Qel(iz)/Aul)*(1.d0 + (Hup2(iz)*dble(iq))**2)
  enddo!ik
 enddo!iq
 do i = 1, nwl
  do ip = 1, nwl
   do ik = 0,2
    do ikp = 0,2
     kmin = MIN(ik,ikp)
     do iq = -kmin, kmin 
      gre = 0.d0; gim = 0.d0
       do ku = 0,ju2
        a1 = 1.d0/denom1(ku,iq) - 1.d0/denom2(iq)
        b1 = Hup2(iz)*dble(iq)/denom2(iq) - Hup(ku,iz)*dble(iq)/denom1(ku,iq) 
        call comprod(a1,b1,gprof_re(ku,ikp,iq,iz,i),gprof_im(ku,ikp,iq,iz,i),c,d)
        call comprod(c,d,gprof_re(ku,ik,iq,iz,ip),gprof_im(ku,ik,iq,iz,ip),e,f)
        gre = gre + e
        gim = gim + f
       enddo !ku
      r3(1,ik,ikp,iq,i,ip,iz) = gre !real part
      r3(2,ik,ikp,iq,i,ip,iz) = gim !imaginary part
     enddo !iq
    enddo !ikp
   enddo !ik
  enddo !ip
 enddo !i
enddo!iz
call cpu_time(t2)
print*, 'R3 time', t2 - t1
!-------------------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE crd_red
!=========================================================================


!=========================================================================
SUBROUTINE coh_red(epsprm,Qel,Hup2,gprof_re,gprof_im,r2)
!-------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: epsprm(1:nz),Qel(1:nz),Hup2(1:nz)
DOUBLE PRECISION, INTENT(IN) :: gprof_re(0:ju2,0:2,-2:2,1:nz,1:nwl)
DOUBLE PRECISION, INTENT(IN) :: gprof_im(0:ju2,0:2,-2:2,1:nz,1:nwl)
DOUBLE PRECISION, INTENT(OUT) :: r2(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz)
INTEGER :: iz,i,j,iq,kp,kpp,qmax,nt,it,mu2,mup2,ml2,mlp2,ip
INTEGER :: is, dj,ja,jb,wnorm,hnorm,fl1,npl,nph,np,testip,testi,testiz!TESTS
INTEGER :: nlt,st,ov,ib,ib2,flt,ij,nwf,dp,pa,pb,jp,j1,j2
DOUBLE PRECISION :: ac,bc,vs,trll,fn
DOUBLE PRECISION, DIMENSION(30) :: tllp
DOUBLE PRECISION :: qden,vta,vtb,delnu,delx,xa,xb,vph,vpl,Gmax,fact1
DOUBLE PRECISION :: dvl,dvh,vul,vpulp,nllp,ai,ar,z1,z2,gr,gi,acr,aci
DOUBLE PRECISION :: norm,dnorm,mnorm,tst,sr,si,rr,ri,ttr,tti
DOUBLE PRECISION ::f1,f2,c1,c2,fac,arg,phi,psi,val,step
DOUBLE PRECISION :: normt,dnormt,mnormt,thres,ctf,dv,dwf,vp,awr1,awr2,awi1,awi2
DOUBLE PRECISION, DIMENSION(0:2) :: wfac
DOUBLE PRECISION, DIMENSION(-2:2,1:2) :: qfac
DOUBLE PRECISION, DIMENSION(50) :: nut,vt
DOUBLE PRECISION, DIMENSION(:),allocatable :: f0,wf,mkr,mki
DOUBLE PRECISION, DIMENSION(:), allocatable :: v,sp,sp2
DOUBLE PRECISION, DIMENSION(:), allocatable :: splsum
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), allocatable :: Aquan,Awe
DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), allocatable :: Agrid,rfin
DOUBLE PRECISION :: tima,timb,testmax,testp
testmax = 0.d0
!-------------------------------------------------------------------------
call cpu_time(tima)
allocate(v(nwl),f0(1:nwl))
allocate(Aquan(1:nwl,-ju2:ju2,-jl2:jl2,-jl2:jl2,1:2))
allocate(Awe(1:nwl,-ju2:ju2,-jl2:jl2,-jl2:jl2,1:2))
allocate(rfin(0:2,0:2,-2:2,1:nwl,1:nwl,1:2))
mnormt =0.d0
nt=0
do mu2 = -ju2,ju2,2
 do ml2 = -jl2,jl2,2
 IF(abs(mu2-ml2).gt.2) CYCLE
 nt = nt+1
 nut(nt) = nu0 + nular*0.5d0*(gu*dble(mu2)-gl*dble(ml2)) !For every transition
 enddo!ml2
enddo!mu2
!-------------------------------------------------------------------------------------------
! Initial
!------------------------------------------------------------------------------------------
thres = 0.01d0
ctf = 1.d-8
do kp=0,2
 fac = 1.d0
 if(mod(2+ju2+jl2,4).ne.0) fac = -1.d0
 val = w6js(2,2,2*kp,ju2,ju2,jl2)
 wfac(kp) = fac*dsqrt(3.d0*(dble(ju2)+1.d0))*val
enddo
mnorm = 0.d0
wnorm = 0.d0
hnorm = 0.d0
rfin = 0.d0
!------------------------------------------------------------------------------------------------
! Begin height loop 
!------------------------------------------------------------------------------------------------
do iz=1,nz
 IF(MOD(iz,5).eq.0) print*, 'height',iz,ht(iz)
 do iq=-2,2
  qden = (1.d0 + (dble(iq)*Hup2(iz))**2)*(1.d0 + epsprm(iz) + Qel(iz)/Aul)
  qfac(iq,1) = 1.d0/qden
  qfac(iq,2) = -dble(iq)*Hup2(iz)/qden
 enddo!iq
!------------------------------------------------------------------------------------------------
!Calculating the reduced frequencies at every height (change w/ Doppler width)
!------------------------------------------------------------------------------------------------
 do i=1,nwl
  v(i) = (nu0 - nu(i))/dnd(iz)
 enddo
 delx = dabs(delnu)/dnd(iz)!For the energy difference of the most separate lower
 if(delx.gt.12.d0) then
  print*, 'Energy level difference between Zeeman components too large'
  stop
 endif
!-------------------------------------------------------------------------------------------------
! Compute the A-weights. First we need to define them for each set
! of quantum numbers, taking into account the selection rules. 
!--------------------------------------------------------------------------------------------------
 do i = 2,nwl-1
 Aquan = 0.d0
 Awe = 0.d0
 f0(i) = gprof_re(0,0,0,iz,i)
  do mu2 = -ju2,ju2,2
   do ml2 = -jl2,jl2,2
    IF(ABS(mu2-ml2).gt.2) CYCLE
    vul = (nu0 + nular*0.5d0*(gu*dble(mu2)-gl*dble(ml2)) - nu(i))/dnd(iz)
    do mlp2 = -jl2,jl2,2
     IF(ABS(mu2-mlp2).gt.2) CYCLE
     nllp = nular*0.5d0*(gl*dble(ml2)-gl*dble(mlp2))/dnd(iz)
!--------------------------------------------------------------------------------------------------
!    In order to find the cutoff, finds grid point for which x_j minimizes the
!    argument of the Gaussian in Aquan (accounting for Zeeman splitting and Raman scattering)
!--------------------------------------------------------------------------------------------------
     if(ml2.eq.mlp2.or.nllp.eq.0.d0) then
      jp = i
     else
      call get_center(v,i,nllp,jp)
     endif
!--------------------------------------------------------------------------------------------------
!    For every j incoming frequency the value of the real and imaginary part of
!    A must be computed. (Calculated separately from the weighs because we may want to put
!    a cuttoff in the calculation later)
!--------------------------------------------------------------------------------------------------
     vpulp = v(jp) + (nular*0.5d0*(gu*dble(mu2)-gl*dble(mlp2)))/dnd(iz)
     call red_aquan(iz,v(i),v(jp),nllp,vul,vpulp,ar,ai)
     acr = ar
     aci = ai
     do is = 1,2
      if(is.eq.1) then
       ja = jp
       jb = nwl
       dj = 1
      else
       ja = jp-1
       jb = 1
       dj = -1
      endif
      do j=ja,jb,dj
       vpulp = v(j) + (nular*0.5d0*(gu*dble(mu2)-gl*dble(mlp2)))/dnd(iz)
       call red_aquan(iz,v(i),v(j),nllp,vul,vpulp,ar,ai)
       Aquan(j,mu2,ml2,mlp2,1) = ar
       Aquan(j,mu2,ml2,mlp2,2) = ai
       IF(is.eq.1) j2 = j
       IF(is.eq.2) j1 = j
       IF(dabs(ar/acr).lt.ctf .and. dabs(ai/aci).lt.ctf) EXIT
      enddo!j
     enddo!is
!-------------------------------------------------------------------------------------------------
!     Now we calculate the weights. Here we could include the ranges considering the cutoff
!     The weight now consider the trapezoidal rule. Here we start from the lowest j value and go up
!--------------------------------------------------------------------------------------------------
     do j=j1,j2
      dv = v(j+1)-v(j)
      if(dv.lt.thres) then
       Awe(j,mu2,ml2,mlp2,1) = Awe(j,mu2,ml2,mlp2,1) + &
        Aquan(j,mu2,ml2,mlp2,1)*dv/2.d0
       Awe(j,mu2,ml2,mlp2,2) = Awe(j,mu2,ml2,mlp2,2) + &
        Aquan(j,mu2,ml2,mlp2,2)*dv/2.d0
       Awe(j+1,mu2,ml2,mlp2,1) = Awe(j+1,mu2,ml2,mlp2,1) + &
        Aquan(j+1,mu2,ml2,mlp2,1)*dv/2.d0
       Awe(j+1,mu2,ml2,mlp2,2) = Awe(j+1,mu2,ml2,mlp2,2) + &
        Aquan(j+1,mu2,ml2,mlp2,2)*dv/2.d0
      else
       !Define range of smaller grid
       nwf = NINT(dv/thres)
       allocate(mkr(nwf+1),mki(nwf+1))
       dwf = dv/nwf
       mkr = 0.d0; mki = 0.d0
       mkr(1) = Aquan(j,mu2,ml2,mlp2,1)
       mkr(nwf+1) = Aquan(j+1,mu2,ml2,mlp2,1)
       mki(1) = Aquan(j,mu2,ml2,mlp2,2)
       mki(nwf+1) = Aquan(j+1,mu2,ml2,mlp2,2)
       !Reorder the limits of the loop so a cutoff can be introduced if desired
       if(mkr(1).ge.mkr(nwf+1)) then
        dp = 1
        pa = 2
        pb = nwf
       else
        dp = -1
        pa = nwf
        pb = 2
       endif
       !Calculate Voigt and Faraday-Voigt functions at intermediate points of subgrid
       do ip=pa,pb,dp
        vp = v(j)+(ip-1)*dwf
        vpulp = vp + (nular*0.5d0*(gu*dble(mu2)-gl*dble(mlp2)))/dnd(iz)
        call red_aquan(iz,v(i),vp,nllp,vul,vpulp,ar,ai)
        mkr(ip) = ar
        mki(ip) = ai
        if(dabs(ar/acr).le.ctf.and.dabs(ai/aci).le.ctf) EXIT
       enddo!ip 
       awr1 = mkr(1)*dv
       awr2 = mkr(nwf+1)*dv
       awi1 = mki(1)*dv
       awi2 = mki(nwf+1)*dv
       do ip=pa,pb,dp
        awr1 = awr1 + 2.d0*mkr(ip)*(nwf-ip+1)*dwf!
        awr2 = awr2 + 2.d0*mkr(ip)*(ip-1)*dwf!
        awi1 = awi1 + 2.d0*mki(ip)*(nwf-ip+1)*dwf!
        awi2 = awi2 + 2.d0*mki(ip)*(ip-1)*dwf!
       enddo
!
       Awe(j,mu2,ml2,mlp2,1) = Awe(j,mu2,ml2,mlp2,1)+awr1/2.d0/nwf
       Awe(j,mu2,ml2,mlp2,2) = Awe(j,mu2,ml2,mlp2,2)+awi1/2.d0/nwf
       Awe(j+1,mu2,ml2,mlp2,1) = Awe(j+1,mu2,ml2,mlp2,1)+awr2/2.d0/nwf
       Awe(j+1,mu2,ml2,mlp2,2) = Awe(j+1,mu2,ml2,mlp2,2)+awi2/2.d0/nwf
       deallocate(mkr,mki)
      endif
     enddo!j
    enddo!mlp2
   enddo!ml2
  enddo!mu2
!-------------------------------------------------------------------------------------------------------
! A quantities have now been calculated for all sets of allowed transitions, for all incoming
! frequencies.  Now we can obtain the r^{K' K''}_{Q}(x_i,xp_j) for all indices.
!-------------------------------------------------------------------------------------------------------
  Gmax = 0.d0
  do is = 1,2
   if(is.eq.1) then
    ja = i
    jb = nwl
    dj = 1
   else
    ja = i-1
    jb = 1
    dj = -1
   endif
   do j=ja,jb,dj
    do kp = 0,2
     do kpp = 0,2
      qmax = min(kp,kpp)
      do iq = -qmax,qmax
       gr = 0.d0
       gi = 0.d0
       do mu2 = -ju2,ju2,2
        do mup2 = -ju2,ju2,2
         do ml2 = -jl2,jl2,2
          do mlp2 = -jl2,jl2,2
           tst = cfac(kp,kpp,iq,mup2,mu2,ml2,mlp2)
           if(tst.eq.0) CYCLE
           rr = tst*0.5d0*(Awe(j,mu2,ml2,mlp2,1) + Awe(j,mup2,ml2,mlp2,1))
           ri = tst*0.5d0*(Awe(j,mu2,ml2,mlp2,2) - Awe(j,mup2,ml2,mlp2,2))
           sr = qfac(iq,1)
           si = qfac(iq,2)
           call comprod(sr,si,rr,ri,ttr,tti)
           gr = gr + ttr
           gi = gi + tti
          enddo!mlp2
         enddo!ml2
        enddo !mup2
       enddo !mu2
       rfin(kp,kpp,iq,i,j,1) = gr
       rfin(kp,kpp,iq,i,j,2) = gi
!
       if(Gmax.lt.gr) Gmax = gr
      enddo!iq
     enddo!kpp
    enddo!kp
   enddo !j index for v'
  enddo!is
  do ip=1,np
   testp = dabs(1.d0-splsum(ip))
   if(testp.gt.testmax) then
    testmax = testp
    testiz = iz
    testi = i
    testip = ip
   endif
  enddo
!------------------------------------------------------------------------------------------------
! Changing from reduced frequency to frequency differentials
!------------------------------------------------------------------------------------------------
  do j = 1,nwl
   do kp = 0,2
    do kpp = 0,2
     qmax = min(kp,kpp)
     do iq = -qmax,qmax
      rfin(kp,kpp,iq,i,j,1) = rfin(kp,kpp,iq,i,j,1)/dnd(iz)
      rfin(kp,kpp,iq,i,j,2) = rfin(kp,kpp,iq,i,j,2)/dnd(iz)
     enddo
    enddo
   enddo
  enddo
!------------------------------------------------------------------------------------------------
! Normalization (in terms of the 000 component). 
!------------------------------------------------------------------------------------------------
  fact1 = (1.d0 + epsprm(iz) + Qel(iz)/Aul)/(1.d0 + epsprm(iz))
  norm=0.d0
  do j=1,nwl
   norm = norm + fact1*(1.d0+epsprm(iz))*rfin(0,0,0,i,j,1)/f0(i)
  enddo!j
  dnorm = dabs(1.d0-norm)
  if(dnorm.gt.4.d-2) print*,'!!',norm,i,iz
  if(dnorm.gt.mnorm) then
   mnorm = dnorm
   wnorm = i
   hnorm = iz
  endif
  do j=1,nwl
   do kp=0,2
    do kpp=0,2
     qmax = min(kp,kpp)
     do iq=-qmax,qmax
      r2(1,kp,kpp,iq,i,j,iz) = rfin(kp,kpp,iq,i,j,1)/norm!output
      r2(2,kp,kpp,iq,i,j,iz) = rfin(kp,kpp,iq,i,j,2)/norm
     enddo!iq
    enddo!kpp
   enddo!kp
  enddo !j
 enddo !i
!
 print*, 'Normalization:', mnorm,'Height: ', ht(hnorm), 'Red. frequency: ',v(wnorm)
!------------------------------------------------------------------------------------------------
! R2 values at first and last frequency grid point
!------------------------------------------------------------------------------------------------
 do kp = 0,2
  do kpp = 0,2
   qmax = min(kp,kpp)
   do iq = -qmax,qmax
    call comprod(qfac(iq,1),qfac(iq,2),gprof_re(kp,kpp,iq,iz,1),&
         gprof_im(kp,kpp,iq,iz,1),c1,c2)
    r2(1,kp,kpp,iq,1,1,iz) = c1
    r2(2,kp,kpp,iq,1,1,iz) = c2
    do j=2,nwl
     r2(1,kp,kpp,iq,1,j,iz) = 0.d0
     r2(2,kp,kpp,iq,1,j,iz) = 0.d0
    enddo
   enddo
  enddo
 enddo
!
 do kp = 0,2
  do kpp = 0,2
   qmax = min(kp,kpp)
   do iq = -qmax,qmax
    call comprod(qfac(iq,1),qfac(iq,2),gprof_re(kp,kpp,iq,iz,nwl),&
         gprof_im(kp,kpp,iq,iz,nwl),c1,c2)
    r2(1,kp,kpp,iq,nwl,nwl,iz) = c1
    r2(2,kp,kpp,iq,nwl,nwl,iz) = c2
    do j=1,nwl-1
     r2(1,kp,kpp,iq,nwl,j,iz) = 0.d0
     r2(2,kp,kpp,iq,nwl,j,iz) = 0.d0
    enddo
   enddo
  enddo
 enddo
enddo!iz
call cpu_time(timb)
!-------------------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE coh_red
!=========================================================================

!=========================================================================
SUBROUTINE get_center(v,i,nllp,jp)
!------------------------------------------------------------------------------
! We want to find the jp that minimizes xj-xi+nllp 
! i.e. the v(j) closest to v(i)-nllp
!------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: v(nwl)
DOUBLE PRECISION, INTENT(IN) :: nllp
INTEGER, INTENT(IN) :: i
INTEGER, INTENT(OUT) :: jp
DOUBLE PRECISION :: xpr
INTEGER :: j,dx,jf
!------------------------------------------------------------------------------
 xpr = v(i)-nllp
 if(nllp.gt.0d0) then
  dx = -1
  jf = 1
 else
  dx = 1
  jf = nwl
 endif
do j = i,jf,dx
 if(DABS(v(j)-xpr).le.DABS(v(j+dx)-xpr)) EXIT
enddo
jp = j
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END SUBROUTINE get_center
!=========================================================================


!=========================================================================
SUBROUTINE red_aquan(iz,v,vp,nllp,vul,vpulp,ar,ai)
INTEGER, INTENT(IN) :: iz
DOUBLE PRECISION, INTENT(IN) :: v,vp,nllp,vul,vpulp
DOUBLE PRECISION, INTENT(OUT) :: ai, ar
INTEGER :: i
DOUBLE PRECISION :: arg,y0,b,b2, quantr,quanti,sen,sen2,dif,ag,ap
DOUBLE PRECISION :: Hp,Lp,expon
!-------------------------------------------------------------------------
arg = 0.5d0*(vul+vpulp)
y0 = 1/(2.d0*pi)!1/2 because of angle-averaging
b = 0.5d0*(v-vp+nllp) 
b2 = b*b
quantr = 0.d0
quanti = 0.d0
!----------------------------------------------------------------------------------------
! Quadrature
!----------------------------------------------------------------------------------------
if(dabs(b).le.1.d-20) then !Core
 do i=1,nqglc
  sen2 = (1.d0-yglc(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  ag = arg/yglc(i)
  ap = ad(iz)/yglc(i)
  call profilarch(ap,ag,Hp,Lp)
  quantr = quantr + y0*dif*Hp*wglc(i)!/pnorm
  quanti = quanti + y0*dif*Lp*wglc(i)!/pnorm
 enddo
!---------------------------------------------------------------------------------------
elseif(dabs(b).le.1.5d0) then !Intermediate interval with higher number of angular points
 do i=1,nqglf !Less points: function is smoother
  sen2 = (1.d0 - yglf(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  expon = dexp(-b2/sen2)
  ag = arg/yglf(i)
  ap = ad(iz)/yglf(i)
  call profilarch(ap,ag,Hp,Lp)!No sqrt(pi) because we find it by integr. in freq. directly
  quantr = quantr + y0*dif*expon*Hp*wglf(i)!/pnorm
  quanti = quanti + y0*dif*expon*Lp*wglf(i)!/pnorm
 enddo!i
!---------------------------------------------------------------------------------------
else
 do i=1,nqgl !Less points: function is smoother
  sen2 = (1.d0 - ygl(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  expon = dexp(-b2/sen2)
  ag = arg/ygl(i)
  ap = ad(iz)/ygl(i)
  call profilarch(ap,ag,Hp,Lp)
  quantr = quantr + y0*dif*expon*Hp*wgl(i)!/pnorm
  quanti = quanti + y0*dif*expon*Lp*wgl(i)!/pnorm
 enddo!i
endif!core or otherwise
!----------------------------------------------------------------------------------------
ar = quantr
ai = quanti
!---------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE red_aquan
!=========================================================================


!=========================================================================
SUBROUTINE fr2_angav(dmp,x,xp,vul,nllp,ar,ai)
DOUBLE PRECISION, INTENT(IN) :: dmp,x,xp,vul,nllp
DOUBLE PRECISION, INTENT(OUT) :: ai, ar
INTEGER :: i
DOUBLE PRECISION :: arg,y0,b,b2, quantr,quanti,sen,sen2,dif,ag,ap
DOUBLE PRECISION :: Hp,Lp,expon
!-------------------------------------------------------------------------
y0 = 1/(2.d0*pi)!1/2 due to angle average
arg = 0.5d0*(x+xp+2.d0*vul+nllp)
b = 0.5d0*(x-xp+nllp)
b2 = b*b
quantr = 0.d0
quanti = 0.d0
!----------------------------------------------------------------------------------------
! Quadrature
!----------------------------------------------------------------------------------------
if(dabs(b).le.1.d-20) then !Core
 do i=1,nqglc
  sen2 = (1.d0-yglc(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  ag = arg/yglc(i)
  ap = dmp/yglc(i)
  call profilarch(ap,ag,Hp,Lp)
  quantr = quantr + y0*dif*Hp*wglc(i)!/pnorm
  quanti = quanti + y0*dif*Lp*wglc(i)!/pnorm
 enddo
!---------------------------------------------------------------------------------------
elseif(dabs(b).le.1.5d0) then !Intermediate interval with higher number of angular points
 do i=1,nqglf !Less points: function is smoother
  sen2 = (1.d0 - yglf(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  expon = dexp(-b2/sen2)
  ag = arg/yglf(i)
  ap = dmp/yglf(i)
  call profilarch(ap,ag,Hp,Lp)!No sqrt(pi) because we find it by integr. in freq. directly
  quantr = quantr + y0*dif*expon*Hp*wglf(i)!/pnorm
  quanti = quanti + y0*dif*expon*Lp*wglf(i)!/pnorm
 enddo!i
!---------------------------------------------------------------------------------------
else
 do i=1,nqgl !Less points: function is smoother
  sen2 = (1.d0 - ygl(i)**2)
  sen = dsqrt(sen2)
  dif = 2.d0/sen
  expon = dexp(-b2/sen2)
  ag = arg/ygl(i)
  ap = dmp/ygl(i)
  call profilarch(ap,ag,Hp,Lp)
  quantr = quantr + y0*dif*expon*Hp*wgl(i)!/pnorm
  quanti = quanti + y0*dif*expon*Lp*wgl(i)!/pnorm
 enddo!i
endif!core or otherwise
!----------------------------------------------------------------------------------------
ar = quantr
ai = quanti
!----------------------------------------------------------------------------------------
return
!----------------------------------------------------------------------------------------
END SUBROUTINE fr2_angav
!========================================================================================

!========================================================================================
SUBROUTINE red_rad(k_l,jrad,r2,r3,jint)
!----------------------------------------------------------------------------------------
! Subroutine that calculates the frequency-integrated radiation tensors. 
!----------------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: k_l(1:nz)
DOUBLE PRECISION, INTENT(IN) :: jrad(1:2,0:2,0:2,1:nz,1:nwl)!If its a global variable, can be taken out
DOUBLE PRECISION, INTENT(IN) :: r3(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz)
DOUBLE PRECISION, INTENT(IN) :: r2(1:2,0:2,0:2,-2:2,1:nwl,1:nwl,1:nz)
DOUBLE PRECISION, INTENT(OUT) :: jint(1:2,0:2,0:2,1:nz,1:nwl)
INTEGER :: i,ipr,iz, ik,iq,iqp, kmain,qmain, qmin,mu2,ml2
DOUBLE PRECISION :: rl1,rl2,rl3,im1,im2,im3,rint,iint,arg,f1,f2,phi,psi
DOUBLE PRECISION :: rsum, isum, rneg,ineg,fac,ikl,rar,rai,f0,norm,dnorm,mnorm
!-------------------------------------------------------------------------
mnorm = 0.d0
do iz = 1, nz
 ikl = k_l(iz)
 do i = 1, nwl
 norm = 0.d0
  do kmain = 0, 2
   do qmain = 0, kmain!Though it probably doesn't matter
    rsum = 0.d0
    isum = 0.d0
    do ik = 0, 2 !J^K_Q radiation tensor cannot have a K index larger than 2
     do iqp = 0, ik! The radiation tensor array has only been defined for positive values of q
      qmin = MIN(ik,kmain)
      do iq = -qmin,qmin
       rint = 0.d0
       iint = 0.d0
       do ipr = 1,nwl !Integral over frequencies
        rar = r3(1,ik,kmain,iq,i,ipr,iz)*uj(ipr) + r2(1,ik,kmain,iq,i,ipr,iz)
        rai = r3(2,ik,kmain,iq,i,ipr,iz)*uj(ipr) + r2(2,ik,kmain,iq,i,ipr,iz)
        call comprod(jrad(1,ik,iqp,iz,ipr), -jrad(2,ik,iqp,iz,ipr),rar,rai,rl1,im1)
        call comprod(rl1,im1,DR(ik,iq,iqp),DI(ik,iq,iqp),rl2,im2)
        call comprod(rl2,im2,DR(kmain,iq,qmain),-DI(kmain,iq,qmain),rl3,im3)
        rint = rint + rl3
        iint = iint + im3
       enddo!ipr
       rneg=0.d0; ineg = 0.d0
       IF(iqp.ne.0) THEN !Finds contribution from negative Q' terms
        do ipr = 1,nwl !Integral over frequencies 
         rar = r3(1,ik,kmain,iq,i,ipr,iz)*uj(ipr) + r2(1,ik,kmain,iq,i,ipr,iz)
         rai = r3(2,ik,kmain,iq,i,ipr,iz)*uj(ipr) + r2(2,ik,kmain,iq,i,ipr,iz)
         call comprod(jrad(1,ik,iqp,iz,ipr),jrad(2,ik,iqp,iz,ipr),rar,rai,rl1,im1)
         call comprod(rl1,im1,DR(ik,iq,-iqp),DI(ik,iq,-iqp),rl2,im2)
         call comprod(rl2,im2,DR(kmain,iq,qmain),-DI(kmain,iq,qmain),rl3,im3)
         IF(MOD(iqp,2).NE.0) THEN ! (-1)^Q'
          fac=-1.d0
         ELSE
          fac = 1.d0
         ENDIF! (-1)^Q'
         rneg = rneg + fac*rl3
         ineg = ineg + fac*im3
        enddo!ipr
       ENDIF! Negative Q'
       rsum = rsum + ikl*rint + ikl*rneg
       isum = isum + ikl*iint + ikl*ineg
      enddo !iq
     enddo !iqp
    enddo!ik
    jint(1,kmain,qmain,iz,i) = rsum
    jint(2,kmain,qmain,iz,i) = isum
   enddo !qmain
  enddo !kmain
 enddo !i
 do i = 1,nwl
  norm = 0.d0
  f0 = gprof_re(0,0,0,iz,i)
  do ipr=1,nwl
   rar = r3(1,0,0,0,i,ipr,iz)*uj(ipr) + r2(1,0,0,0,i,ipr,iz)
   norm = norm + (1.d0+epsprm(iz))*rar/f0
  enddo
  dnorm = dabs(1.d0-norm)
  if(mnorm.lt.dnorm) mnorm = dnorm
 enddo
!-------------------------------------------------------------------------
! if(MOD(iz,10).eq.0) print*, 'Height', iz
enddo !iz
!-------------------------------------------------------------------------
!print*, 'Normalization:', mnorm
!-------------------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE red_rad
!=========================================================================

!=========================================================================
SUBROUTINE lambd_line(k_l,epsprm,Wm,gprof_re,gprof_im,jint,sline)
!-------------------------------------------------------------------------
! Calculates the multipolar components of the line source function in terms
! of the frequency-redistributed source functions
!-------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: k_l(nz),epsprm(nz),Wm(nz)
DOUBLE PRECISION, INTENT(IN) :: gprof_re(0:ju2,0:2,-2:2,nz,nwl), gprof_im(0:ju2,0:2,-2:2,nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: jint(2,0:2,0:2,nz,nwl)
DOUBLE PRECISION, INTENT(OUT) :: sline(2,0:2,0:2,nz,nwl)
DOUBLE PRECISION :: re1, im1, fac
INTEGER :: ik, iq, i, iz
!-------------------------------------------------------------------------
do iz = 1, nz
 fac = k_l(iz)*Wm(iz)*epsprm(iz)/(1.d0 + epsprm(iz))
 do i = 1, nwl
  do ik = 1,2
   do iq = 0, ik
    call comprod(gprof_re(0,ik,0,iz,i), gprof_im(0,ik,0,iz,i), DR(ik,0,iq),-DI(ik,0,iq),re1,im1)
    sline(1,ik,iq,iz,i) = fac*re1 + jint(1,ik,iq,iz,i)
    if(iq.ne.0) then
    sline(2,ik,iq,iz,i) = fac*im1 + jint(2,ik,iq,iz,i)
    else
     sline(2,ik,iq,iz,i) = 0.d0!Since S^K_-Q = (-1)^Q S^K_Q(nu)^ast
    endif
   enddo !iq
  enddo !ik
 enddo !i
enddo !iz
!-------------------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE lambd_line
!=========================================================================


!=========================================================================
SUBROUTINE lambd_cont(s,eps_c,etaI_c,jrad,scont)
!--------------------------------------------------------------------------
! Calculates the multipolar components of the line source function in terms
! of the frequency-redistributed source functions
!--------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: s(nz,nwl), etaI_c(nz,nwl), eps_c(nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: jrad(2,0:2,0:2,nz,nwl)
DOUBLE PRECISION, INTENT(OUT) :: scont(2,0:2,0:2,nz,nwl)
INTEGER :: i,iz,ik,iq
!--------------------------------------------------------------------------
do iz = 1,nz
 do i = 1, nwl
  do ik = 0, 2
   do iq = 0, ik
    IF(iq.EQ.0) THEN
     IF(ik.EQ.0) THEN
      scont(1,ik,iq,iz,i) = etaI_c(iz,i)*(1.d0-s(iz,i))*jrad(1,ik,iq,iz,i) + eps_c(iz,i)
      scont(2,ik,iq,iz,i) = 0.d0                                          
     ELSE
      scont(1,ik,iq,iz,i) = etaI_c(iz,i)*(1.d0-s(iz,i))*jrad(1,ik,iq,iz,i)
      scont(2,ik,iq,iz,i) = 0.d0
     ENDIF
    ELSE
     scont(1,ik,iq,iz,i) = etaI_c(iz,i)*(1.d0-s(iz,i))*jrad(1,ik,iq,iz,i)
     scont(2,ik,iq,iz,i) = -etaI_c(iz,i)*(1.d0-s(iz,i))*jrad(2,ik,iq,iz,i)
    ENDIF
   enddo !iq
  enddo!ik
 enddo!i
enddo !iz
!-------------------------------------------------------------------------
return
!-------------------------------------------------------------------------
END SUBROUTINE lambd_cont
!=========================================================================

END MODULE redis
