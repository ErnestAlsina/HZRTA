MODULE matrmod
!------------------------------------------------------------------------------
! Module only contains the subroutine for the M matrix which contains height, 
!and two frequency indices. It carries a differential for the nu' frequency 
!(since the expression for Delta S^0_0 for the Jacobi iteration method carries)
!------------------------------------------------------------------------------
USE global_vars
USE parameters
USE math
!------------------------------------------------------------------------------
IMPLICIT NONE
!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------


!==============================================================================
SUBROUTINE m_mat(iz,k_l,r,r2,r3,lstar,eta_re,m)!,kmat,m)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: iz 
DOUBLE PRECISION, INTENT(IN) :: k_l(nz)
DOUBLE PRECISION, INTENT(IN) :: r(nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: r2(2,0:2,0:2,-2:2,nwl,nwl,nz)
DOUBLE PRECISION, INTENT(IN) :: r3(2,0:2,0:2,-2:2,nwl,nwl,nz)
DOUBLE PRECISION, INTENT(IN) :: lstar(0:3,0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: eta_re(0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(OUT) :: m(nwl,nwl)
INTEGER :: ip,ipr,im,ax,ik,iq, in,in1
DOUBLE PRECISION :: x0,fac, term_re,term_im, a_re,a_im, val
DOUBLE PRECISION :: wgt, mi,rai,rar
DOUBLE PRECISION :: angint_re, angint_im, res_re,res_im, b_re, b_im
DOUBLE PRECISION :: absim, absm!TEST
DOUBLE PRECISION, DIMENSION(0:3,0:3) :: minv
DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: tstvec
INTEGER :: un
!------------------------------------------------------------------------------
absim = 0.d0!TEST
absm = 0.d0!TEST2
do ip = 1,nwl !i frequency: Outgoing. No differential
 do ipr = 1,nwl!j frequency: Incoming. Integrated over in S^0_0 
  res_re = 0.d0; res_im = 0.d0
  do ik = 0,2
   do iq = -ik,ik
    if(MOD(IQ,2).ne.0) then
     fac = -1.d0
    else
     fac = 1.d0
    endif
    rar = fac*r3(1,ik,0,0,ip,ipr,iz)*uj(ipr) + fac*r2(1,ik,0,0,ip,ipr,iz)!Should I account for the reduced freq. differential?
    rai = fac*r3(2,ik,0,0,ip,ipr,iz)*uj(ipr) + fac*r2(2,ik,0,0,ip,ipr,iz)
    call comprod(rar,rai,DR(ik,0,iq),DI(ik,0,iq),a_re,a_im)
    angint_re = 0.d0; angint_im = 0.d0 !angular integration for every freq. and k', q' indices
    do im = 1, ndir
     mi = mu(im)
     if((mi.le.0 .and. iz.eq.1) .or. (mi.gt.0 .and. iz.eq.nz)) then
     else
      do ax = 1, ndirx
       wgt = wtdir(im)*waz(ax)/(4.d0*pi)! angular weight. Same as 'integrating' over all azimuths and 
!      then multiplying by the polar angle differential. 1/2 due to the mu weights.
       val = wgt*r(iz,ipr,im,ax)/eta_re(0,iz,ipr,im,ax)
       b_re = 0.d0; b_im = 0.d0
       do in = 0,3!has angular dependence
        b_re = b_re + TTR(ik,-iq,in,im,ax)*val*lstar(in,0,iz,ipr,im,ax)
        b_im = b_im + TTI(ik,-iq,in,im,ax)*val*lstar(in,0,iz,ipr,im,ax)
       enddo!in
       angint_re = angint_re + b_re
       angint_im = angint_im + b_im
      enddo!ax
     endif
    enddo !im 
    call comprod(a_re,a_im,angint_re,angint_im, term_re,term_im)
    res_re = res_re + term_re
    res_im = res_im + term_im
   enddo! iq
  enddo!ik
  if(DABS(res_im).gt.absim) absim = DABS(res_im)!TEST
  if(ipr.eq.ip) then
   x0 = 1.d0
   else
   x0 = 0.d0
  endif
  m(ip,ipr) = x0 - k_l(iz)*res_re!
  if(DABS(m(ip,ipr)).gt.absm) absm = DABS(m(ip,ipr))!TEST2
 enddo!ipr
enddo !ip
!print*,'Imaginary part:', absim
!print*, 'Largest m value (absolute):', absm 
!------------------------------------------------------------------------------
return
!------------------------------------------------------------------------------
END SUBROUTINE m_mat
!==============================================================================


!==============================================================================
SUBROUTINE find_rho(k_l,etaI,lstar,dline,add,r2,r3,rhol)
!------------------------------------------------------------------------------
DOUBLE PRECISION, INTENT(IN) :: k_l(nz)
DOUBLE PRECISION, INTENT(IN) :: etaI(nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: lstar(0:3,0:3,nz,nwl,ndir,ndirx)
DOUBLE PRECISION, INTENT(IN) :: dline(1:2,0:2,0:2,nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: add(nz,nwl)
DOUBLE PRECISION, INTENT(IN) :: r2(2,0:2,0:2,-2:2,nwl,nwl,nz)
DOUBLE PRECISION, INTENT(IN) :: r3(2,0:2,0:2,-2:2,nwl,nwl,nz)
DOUBLE PRECISION, INTENT(OUT) :: rhol(nz,nwl)
INTEGER :: iz,ip,ipr,si,sj,im,ax,ik,iq,ikp,iqp
DOUBLE PRECISION :: tone,ttwo,prt,cntrl,ddr,ddi,rr,ri,roint,wgt,val,qfac
DOUBLE PRECISION :: brt_re,brt_im,stprod,xr,xi,vr,vi,f1,f2,f3,f4,a_re,a_im
DOUBLE PRECISION :: pkq_re(0:2,-2:2,nwl,ndir,ndirx),pkq_im(0:2,-2:2,nwl,ndir,ndirx)
DOUBLE PRECISION,allocatable :: psum(:,:,:,:)
DOUBLE PRECISION :: rar(0:2,-2:2),rai(0:2,-2:2)
!----------------------------------------------------------------------
call cpu_time(tone)
allocate(psum(0:3,1:nwl,1:ndir,1:ndirx))
do iz = 1,nz
 pkq_re = 0.d0
 pkq_im = 0.d0
 psum = 0.d0
 do ipr = 1,nwl
  do im = 1,ndir
  IF((mu(im).le.0 .and. iz.eq.1) .or. (mu(im).gt.0 .and. iz.eq.nz)) THEN
     !nothing. There is no lstar at these points/directions
  ELSE
   do ax = 1, ndirx
    do si = 0,3
     prt = 0.d0
     do ikp = 1,2
      do iqp = 0,ikp
       cntrl = 2.d0
       if(iqp.eq.0) cntrl = 1.d0
       ddr = dline(1,ikp,iqp,iz,ipr)
       ddi = dline(2,ikp,iqp,iz,ipr)
       call comprod(TTR(ikp,iqp,si,im,ax),TTI(ikp,iqp,si,im,ax),ddr,ddi,rr,ri)
       prt = prt + cntrl*rr
      enddo!iqp
     enddo!ikp
     psum(si,ipr,im,ax) = prt
    enddo!si
    brt_re = 0.d0; brt_im = 0.d0
    do si = 0,3
     do sj = 0,3
      stprod = lstar(si,sj,iz,ipr,im,ax)*psum(sj,ipr,im,ax)
      do ik = 0,2
       do iq = -ik,ik
        pkq_re(ik,iq,ipr,im,ax) = pkq_re(ik,iq,ipr,im,ax) + &
        TTR(ik,-iq,si,im,ax)*stprod
        pkq_im(ik,iq,ipr,im,ax) = pkq_im(ik,iq,ipr,im,ax) + &
        TTI(ik,-iq,si,im,ax)*stprod
       enddo!iq
      enddo!ik
     enddo!sj
    enddo!si
   enddo!ax
  ENDIF
  enddo!im
 enddo!ipr
!------------------------------------------------------------------------------
 do ip = 1,nwl
  roint = 0.d0
  do ipr = 1,nwl
   do ik = 0,2
    do iq = -ik,ik
     f1 = r3(1,ik,0,0,ip,ipr,iz)*uj(ipr) + r2(1,ik,0,0,ip,ipr,iz)
     f2 = r3(2,ik,0,0,ip,ipr,iz)*uj(ipr) + r2(2,ik,0,0,ip,ipr,iz)
     f3 = DR(ik,0,iq)
     f4 = DI(ik,0,iq)
     call comprod(f1,f2,f3,f4,a_re,a_im)
     rar(ik,iq) = a_re
     rai(ik,iq) = a_im
    enddo
   enddo
   do im = 1,ndir
    IF((mu(im).le.0 .and. iz.eq.1) .or. (mu(im).gt.0 .and. iz.eq.nz)) THEN
     !nothing. There is no lstar at these points/directions
    ELSE
     do ax = 1,ndirx
      wgt = wtdir(im)*waz(ax)/(4.d0*pi)
      val = wgt/etaI(iz,ipr,im,ax)
      do ik = 0,2
       do iq = -ik,ik
       qfac = 1.d0
       if(MOD(iq,2).ne.0) qfac = -1.d0
       xr = pkq_re(ik,iq,ipr,im,ax)
       xi = pkq_im(ik,iq,ipr,im,ax)
       call comprod(rar(ik,iq),rai(ik,iq),xr,xi,vr,vi)
       roint = roint + k_l(iz)*qfac*val*vr
       enddo!iq
      enddo!ik
     enddo!ax
    ENDIF
   enddo!im
  enddo !ipr
 rhol(iz,ip) = roint + add(iz,ip)!
 enddo!ip
if(MOD(iz,10).eq.0) print*, 'Height:', iz
enddo!iz
call cpu_time(ttwo)
print*, 'Time', ttwo-tone
!----------------------------------------------------------------------
return
!----------------------------------------------------------------------
END SUBROUTINE find_rho
!==============================================================================

END MODULE matrmod
