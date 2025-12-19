!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Developd by Neil Ribe (Laboratoire FAST, Universit\'e Paris-Saclay, CNRS, Orsay, 91405, France)
! Minor modifications by Manuele Faccenda (Dipartimento di Geoscienze, Universita' di Padova)   
!
! Related article: 
! Ribe, N., Faccenda, M., VanderBeek, B.P. ODFTEX: A continuum model for texture evolution with
! dynamic recrystallization. Submitted to GJI, 2025-2026
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE ODFDRX(odfmode,tid,m,fractdisl,f1_ol,f1_opx,acs_ol,acs_opx)

   USE comvar

!  Calculates evolution of ODF with active drx (new one-parameter model)
!  Accepts an arbitrary velocity gradient tensor (assumed constant in time)
!  Operative slip system assumed (only one for now):
!     s=1: (010)[100] (olivine)
!  The limitation to one slip system means that the only nonzero spin is \dot\psi
!  Time is scaled by 1/Sqrt[E_{ij} E_{ij}/2]

!  Input parameters:
!     nbox = number of intervals into which the ranges of ph \in [0,pi], cos(th) \in[-1,1]
!       and ps \in [0,pi] are divided. The total number of boxes in the discretized Euler 
!       space is nbox**3.
!  Output: 
!     f_ol = ODF for olivine at the centers of each  box in the Euler space
!     f_opx = ODF for opx
!
   IMPLICIT NONE

   INTEGER, intent(in) :: odfmode,tid,m
   DOUBLE PRECISION, intent(in) :: fractdisl
   DOUBLE PRECISION, DIMENSION(nboxnum), intent(out) ::  f1_ol,f1_opx
   DOUBLE PRECISION, DIMENSION(nboxnum,3,3), intent(out) ::  acs_ol,acs_opx

   INTEGER :: iph,ith,ips,i,j,k,ijk,n,i1,i2,i3,nnum,ntint,s_ol,ideg,ilimit,isteptype
   DOUBLE PRECISION :: zero,one,tiny,cf
   DOUBLE PRECISION :: rlam,eii,sumf
   DOUBLE PRECISION :: dtau,dtau2,taumax,taunew,taustart
   DOUBLE PRECISION :: phdot,thdot,psdot,divpsdot
   DOUBLE PRECISION :: c1,drx,frac_opx
   DOUBLE PRECISION :: rmax
   DOUBLE PRECISION :: cph,sph,phdot_ext,thdot_ext,psdot_ext,psdot_int
   DOUBLE PRECISION :: ph,th,ps,dph,dth,dcosth,dps,costh
   DOUBLE PRECISION :: r12,r23,fintegral,lnf
   DOUBLE PRECISION, DIMENSION(3) ::  fse_axis,omdot,evals
   DOUBLE PRECISION, DIMENSION(3,3) ::  o,vg,sr,evects
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: pha,tha,psa,phm,thm,psm        
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: f,f_nodrx,f_ol,f_opx,fnew
   CHARACTER(LEN=60):: fileout_ol, fileout_opx

   zero = 0.0d0
   one = 1.0d0
   tiny = 1.0d-6

   s_ol = 1
   frac_opx = Xol(rocktype(m))
   rlam = lambda(rocktype(m))
   rmax = Mob(rocktype(m))
   ilimit = 0
   ideg = 1

   !Strain rate tensor
   sr = e(tid,:,:)
   !Second invariant of strain rate
   CALL DSYEVQ3(sr,evects,evals)
   eii = sqrt(0.5*(evals(1)**2.0 + evals(2)**2.0 + evals(3)**2.0))

   !Scale strain rate tensor
   sr = sr/eii

   !Scaled velocity gradient tensor
   vg = l(tid,:,:)/eii

! ideg determines whether output Eulerian angles are in degrees or radians
   if (ideg.eq.1) then
     cf = 180.0d0/pi
   else
     cf = one
   endif

! Extrinsic spin (one-half the vorticity)
   omdot(1) = 0.5*(vg(3,2) - vg(2,3))
   omdot(2) = 0.5*(vg(1,3) - vg(3,1))
   omdot(3) = 0.5*(vg(2,1) - vg(1,2))

   !Compute FSE
   CALL eigen(m,Fij(:,:,m),evects,evals,1)
   
   fse_axis(1) = evals(3)
   fse_axis(2) = evals(2)
   fse_axis(3) = evals(1)
   r12 = dlog(fse_axis(1)/fse_axis(2))
   r23 = dlog(fse_axis(2)/fse_axis(3))
   !write(6,*) 'r23, r12 =', r23, r12

   if (r12.gt.rmax.or.r23.gt.rmax) ilimit = 1  ! limit FSE to (r12,r23) < rmax

   !Dimensionless timestep and number of steps
   dtau = MIN(dt,1e-2/eii) 
   ntint = NINT(dt/dtau)
   if(ntint > 100) then
      write(6,*) dt,eii,dtau,ntint
      stop
   end if
   dtau = dt/REAL(ntint)
   dtau = dtau*eii !non-dimensionalize
   dtau = dtau*fractdisl
   dtau2  = dtau/2.0

   c1 = - 2.d0*(sr(1,1)**2 + sr(1,2)**2 + sr(1,1)*sr(2,2) + sr(2,2)**2 + sr(2,3)**2 + sr(3,1)**2)/5.d0

   !Dimensionalize rlam
   !rlam = rlam *eii     

   pha = pham(m,:,:,:)
   tha = tham(m,:,:,:)
   psa = psam(m,:,:,:)
   f   =   fm(m,:,:,:)
   f_nodrx   =   fm_nodrx(m,:,:,:)

!Compute odf evolution
   IF(odfmode == 0) THEN

      DO nnum = 1, ntint  ! start of time-stepping loop

   !taustart = (nt - 1)*dtau
   
! Half step to midpoint tau + dtau/2
         CALL STEP(1,pha,tha,psa,pha,tha,psa,phm,thm,psm,f,f_nodrx,nbox,omdot,sr,rlam,c1,ilimit,dtau2)

! Full step to tau + dtau
         CALL STEP(2,pha,tha,psa,phm,thm,psm,pha,tha,psa,f,f_nodrx,nbox,omdot,sr,rlam,c1,ilimit,dtau)

!  Normalize integral of ODF to unity
         CALL INTEGRAL_DRX(f,f_nodrx,fintegral)
         if(fintegral < 0 .OR. fintegral > 2.0) write(6,*) 'fintegral =', fintegral, m
         DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
           f(i,j,k) = f(i,j,k)/fintegral
         ENDDO; ENDDO; ENDDO

      ENDDO  ! end of time-stepping loop

      !Save the new ODF
      pham(m,:,:,:) = pha
      tham(m,:,:,:) = tha
      psam(m,:,:,:) = psa
      fm(m,:,:,:) = f   
      fm_nodrx(m,:,:,:) = f_nodrx   

!Save ODF to 1D array for Ol and Opx
   ELSE

!Weight the odf to get sum = 1
      CALL SCALE_ODF(f,f_nodrx,fnew)

! the subroutine EULERIAN_TRANSFORM is particularly fast if 
! the only active olivine slip system is s_ol = 1, because this slip
! system defines the reference Eulerian angles.
!  Transform Eulerian angles for ol and opx 
      CALL EULERIAN_TRANSFORM(s_ol,fnew,pha,tha,psa,f1_ol,acs_ol)
      if (frac_opx > 0) CALL EULERIAN_TRANSFORM(4,fnew,pha,tha,psa,f1_opx,acs_opx)

      !Check ODFs               
      DO iph = 1, nbox
      DO ith = 1, nbox
      DO ips = 1, nbox

         i = ips + (ith-1)*nbox + (iph-1)*nbox*nbox

         IF(isnan(f1_ol(i))) then
            print *,"f1_ol Nan at",i,f1_ol(i),m,rocktype(m),mx1(m),mx2(m),mx3(m)
            stop
         END IF
         IF(isnan(f1_opx(i))) then
            print *,"f1_opx Nan at",i,f1_opx(i),m,rocktype(m),mx1(m),mx2(m),mx3(m)
            stop
         END IF

      ENDDO; ENDDO; ENDDO

   END IF

   if(0==1) THEN
      OPEN (31,file='ol.out')
      !OPEN (32,file=fileout_opx)
      DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
         ijk = i + (j-1)*nbox + (k-1)*nbox*nbox
         write(31,1002) cf*pha(i,j,k),cf*tha(i,j,k),cf*psa(i,j,k),f1_ol(ijk)
         !write(32,1002) cf*ph_opx(i,j,k),cf*th_opx(i,j,k),cf*ps_opx(i,j,k),f_opx(i,j,k)
      ENDDO; ENDDO; ENDDO
      CLOSE (31)
      !CLOSE (32)

     1000 format(3i5,4(1pe13.5))
     1001 format(4(1pe14.6))
     1002 format(4(1pe11.3))
   end if

   RETURN

   END SUBROUTINE ODFDRX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE STEP(isteptype,ph1,th1,ps1,phm,thm,psm,ph2,th2,ps2,f,f_nodrx,nbox,omdot,sr,rlam,c1,ilimit,dt)
!  Advances solution by a time dt
!  (ph1,th1,ps1): initial Eulerian angles
!  (phm,thm,psm): Eulerian angles for which RHS is  calculated
!  (ph2,th2,ps2): final Eulerian angles
!  isteptype = 1: Euler half-step by dtau/2 to advance Eulerian angles along characteristics
!  isteptype = 2: full midpoint time step by dtau; also advances the ODF
   IMPLICIT NONE
   INTEGER :: i,j,k,nbox,ilimit,isteptype
   DOUBLE PRECISION :: rlam,c1,dt,cph,sph,phdot_ext,thdot_ext,psdot_ext,drx,psdot_int,phdot, &
    thdot,psdot,divpsdot,zero,lnf,lnf_nodrx
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: ph1,th1,ps1,ph2,th2,ps2,phm,thm,psm,f,f_nodrx        
   DOUBLE PRECISION, DIMENSION(3,3) :: sr
   DOUBLE PRECISION, DIMENSION(3) :: omdot

   zero = 0.0d0

   DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
     cph = dcos(phm(i,j,k)) 
     sph = dsin(phm(i,j,k))
     phdot_ext = omdot(3) + (omdot(2)*cph - omdot(1)*sph)/dtan(thm(i,j,k))
     thdot_ext = omdot(1)*cph + omdot(2)*sph
     psdot_ext = (omdot(1)*sph - omdot(2)*cph)/dsin(thm(i,j,k))
     if (ilimit.eq.1) then
       drx = zero
       psdot_int = zero
       divpsdot = zero
       phdot = phdot_ext
       thdot = thdot_ext
       psdot = psdot_ext
     else
        CALL PSDOT_CALC_GENERAL(phm(i,j,k),thm(i,j,k),psm(i,j,k),sr,psdot_int,divpsdot)
        phdot = phdot_ext
        thdot = thdot_ext
        psdot = psdot_ext + psdot_int
        IF (isteptype.eq.2) then
           drx = rlam*(c1 + 2.d0*psdot_int**2)
           !lnf = dlog(f(i,j,k))
           !lnf = lnf + dt*(drx - divpsdot)
           !f(i,j,k) = dexp(lnf)
           lnf = dt*(drx - divpsdot)
           f(i,j,k) = f(i,j,k)*dexp(lnf)
           !lnf_nodrx = dlog(f_nodrx(i,j,k))
           !lnf_nodrx = lnf_nodrx - dt*divpsdot
           !f_nodrx(i,j,k) = dexp(lnf_nodrx)
           lnf = -dt*divpsdot
           f_nodrx(i,j,k) = f_nodrx(i,j,k)*dexp(lnf)
        ENDIF
     endif
     ph2(i,j,k) = ph1(i,j,k) + dt*phdot
     th2(i,j,k) = th1(i,j,k) + dt*thdot
     ps2(i,j,k) = ps1(i,j,k) + dt*psdot
   ENDDO; ENDDO; ENDDO

   END SUBROUTINE STEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE EULERIAN_TRANSFORM(s,f1,ph1,th1,ps1,f3,acs3)

   USE comvar

   IMPLICIT NONE
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: f1,ph1,th1,ps1        
   DOUBLE PRECISION, DIMENSION(3,3) :: bmat,amat,bamat
   DOUBLE PRECISION, DIMENSION(nboxnum) :: f3
   DOUBLE PRECISION, DIMENSION(nboxnum,3,3) :: acs3
   DOUBLE PRECISION :: zero,one,volfract
   INTEGER :: s,i,j,k,ijk
   zero = 0.0d0
   one = 1.0d0

   IF (s.eq.1) THEN
     DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
       ijk = k + (j-1)*nbox + (i-1)*nbox*nbox
       CALL EULER_TO_DIRCOS(ph1(i,j,k),th1(i,j,k),ps1(i,j,k),acs3(ijk,:,:))
       f3(ijk) = f1(i,j,k)
     ENDDO; ENDDO; ENDDO
     RETURN   
   ENDIF

   DO i=1,3; DO j=1,3
     bmat(i,j) = zero 
   ENDDO; ENDDO
   IF (s.eq.2) THEN
     bmat(1,1) = one
     bmat(2,3) = one
     bmat(3,2) = -one
   ENDIF
   IF (s.eq.3) THEN
     bmat(1,3) = 1.0d0
     bmat(2,2) = -1.0d0
     bmat(3,1) = 1.0d0
   ENDIF
   IF (s.eq.4) THEN
     bmat(1,2) = one
     bmat(2,3) = one
     bmat(3,1) = one
   ENDIF

   DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
     ijk = k + (j-1)*nbox + (i-1)*nbox*nbox
     CALL EULER_TO_DIRCOS(ph1(i,j,k),th1(i,j,k),ps1(i,j,k),amat)
     acs3(ijk,:,:) = MATMUL(bmat(:,:),amat(:,:)) 
     
     f3(ijk) = f1(i,j,k)
   ENDDO; ENDDO; ENDDO

   RETURN
   END SUBROUTINE EULERIAN_TRANSFORM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PSDOT_CALC_GENERAL(ph,th,ps,sr,psdot,divpsdot)
! Calculates intrinsic spin gdot and its derivatives
   IMPLICIT NONE
   DOUBLE PRECISION :: ph,th,ps,psdot,divpsdot
   DOUBLE PRECISION :: cph,sph,cth,sth,cps,sps
   DOUBLE PRECISION, DIMENSION (3,3) :: sr,a,daps
   
   CALL EULER_TO_DIRCOS(ph,th,ps,a)

   cph = dcos(ph)
   sph = dsin(ph)
   cth = dcos(th)
   sth = dsin(th)
   cps = dcos(ps)
   sps = dsin(ps)
! derivatives of a(i,j) w/r to ps
   daps(1,1) = -(cps*cth*sph) - cph*sps 
   daps(1,2) = cph*cps*cth - sph*sps
   daps(1,3) = cps*sth
   daps(2,1) = -(cph*cps) + cth*sph*sps
   daps(2,2) = -(cps*sph) - cph*cth*sps
   daps(2,3) = -(sps*sth)
   daps(3,1) = 0.d0
   daps(3,2) = 0.d0
   daps(3,3) = 0.d0

   psdot = a(1,2)*(a(2,1)*sr(1,2) + a(2,2)*sr(2,2) + a(2,3)*sr(2,3))  &
         + a(1,1)*(a(2,1)*sr(1,1) + a(2,2)*sr(1,2) + a(2,3)*sr(3,1))  &
         + a(1,3)*(a(2,2)*sr(2,3) + a(2,1)*sr(3,1) + a(2,3)*sr(3,3)) 
   divpsdot = daps(1,2)*(a(2,1)*sr(1,2) + a(2,2)*sr(2,2) + a(2,3)*sr(2,3)) +   &
         a(1,2)*(daps(2,1)*sr(1,2) + daps(2,2)*sr(2,2) + daps(2,3)*sr(2,3)) +      &
         daps(1,1)*(a(2,1)*sr(1,1) + a(2,2)*sr(1,2) + a(2,3)*sr(3,1)) +        &
         a(1,1)*(daps(2,1)*sr(1,1) + daps(2,2)*sr(1,2) + daps(2,3)*sr(3,1)) +      &
         daps(1,3)*(a(2,2)*sr(2,3) + a(2,1)*sr(3,1) + a(2,3)*sr(3,3)) +        &
         a(1,3)*(daps(2,2)*sr(2,3) + daps(2,1)*sr(3,1) + daps(2,3)*sr(3,3))

   RETURN
   END SUBROUTINE PSDOT_CALC_GENERAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PSDOT_CALC(ph,th,ps,e1,e2,psdot,divpsdot)
! Calculates psdot and its derivative divpsdot w/r to ps
   IMPLICIT NONE
   DOUBLE PRECISION :: ph,th,ps,e1,e2,e3,psdot,divpsdot 
   DOUBLE PRECISION :: cph,sph,c2ph,s2ph
   DOUBLE PRECISION :: cth,sth,c2th,s2th
   DOUBLE PRECISION :: c2ps,s2ps
   DOUBLE PRECISION :: e1m2,e2m3
   DOUBLE PRECISION :: bigf,bigg,bigh,bigp,bigq

   e3 = - e1 - e2
   cph = dcos(ph)
   sph = dsin(ph)
   c2ph = dcos(2*ph)
   s2ph = dsin(2*ph)
   cth = dcos(th)
   sth = dsin(th)
   c2th = dcos(2*th)
   s2th = dsin(2*th)
   c2ps = dcos(2*ps)
   s2ps = dsin(2*ps)
   e1m2 = e1 - e2
   e2m3 = e2 - e3
   bigf = - s2ph*cth
   bigg = sph**2*cth**2 - cph**2
   bigh = - sth**2
   bigp = e1m2*bigf
   bigq = e1m2*bigg + e2m3*bigh

   psdot = 0.5*(bigp*c2ps + bigq*s2ps)
   divpsdot = -bigp*s2ps + bigq*c2ps

   RETURN
   END SUBROUTINE PSDOT_CALC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FSE_LIMIT(rmax,r12,r23)
!  Limits the FSE to prevent excessively concentrated textures

   IMPLICIT NONE

   DOUBLE PRECISION :: rmax,r12,r23,varphi,pi

   pi = 4.0d0*datan(1.0d0)

   varphi = datan2(r12,r23)
   if (r12.le.r23) then 
     r23 = rmax
     r12 = rmax*dtan(varphi)
   else
     r12 = rmax
     r23 = rmax*dtan(pi/2.d0 - varphi)
   endif

   END SUBROUTINE FSE_LIMIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE EULER_TO_DIRCOS(ph,th,ps,o)
!  Calculate direction cosines corresponding to given Eulerian angles

   IMPLICIT NONE

   DOUBLE PRECISION :: ph,th,ps
   DOUBLE PRECISION :: cph,cth,cps,sph,sth,sps
   DOUBLE PRECISION, DIMENSION(3,3) :: o

   cph = dcos(ph)
   cth = dcos(th)
   cps = dcos(ps)
   sph = dsin(ph)
   sth = dsin(th)
   sps = dsin(ps)

   o(1,1) =   cps*cph  - cth*sph*sps 
   o(1,2) =   cps*sph  + cth*cph*sps
   o(1,3) =   sps*sth 
   o(2,1) = - sps*cph  - cth*sph*cps 
   o(2,2) = - sps*sph  + cth*cph*cps 
   o(2,3) =   cps*sth 
   o(3,1) =   sph*sth 
   o(3,2) = - cph*sth 
   o(3,3) =   cth 

   return

   END SUBROUTINE EULER_TO_DIRCOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DIRCOS_TO_EULER(o,ph,th,ps)
!  Calculates Eulerian angles corresponding to a given matrix of direction cosines

   IMPLICIT NONE

   DOUBLE PRECISION :: ph,th,ps,pi
   DOUBLE PRECISION, DIMENSION(3,3) :: o

   pi = 4.0d0*datan(1.0d0)

   th = dacos(o(3,3))
   ph = datan2(o(3,1),-o(3,2))
   ps = datan2(o(1,3),o(2,3))
 
   return

   END SUBROUTINE DIRCOS_TO_EULER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INTEGRAL_DRX(odff,odff_nodrx,integral)
!  Approximate calculation of (1/(2 Pi**2)) Integral(f dg)
!  Uses dg/odf_nodrx as an estimate of the elemental volume of the Euler
!    space for point (i,j,k)
!  Result = 1 exactly if drx is turned off

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,iphm1,ithm1,ipsm1
   DOUBLE PRECISION :: integral,df,ph_side,x_side,ps_side,xp,xm,box
   DOUBLE PRECISION :: dcosth,dph,dps,dg
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: odff,odff_nodrx

   !pi = 4.0d0*datan(1.0d0)

   dcosth = 2.0d0/nbox
   dph = pi/nbox
   dps = pi/nbox
   dg = dcosth*dph*dps

   integral = 0.0d0
   DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
      integral = integral + dg*odff(i,j,k)/odff_nodrx(i,j,k)
   ENDDO; ENDDO; ENDDO
   integral = integral/(2.0*pi**2)

   END SUBROUTINE INTEGRAL_DRX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SCALE_ODF(odff,odff_nodrx,odffscaled)
!  Approximate calculation of (1/(2 Pi**2)) Integral(f dg)
!  Uses dg/odf_nodrx as an estimate of the elemental volume of the Euler
!    space for point (i,j,k)
!  Result = 1 exactly if drx is turned off

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k
   DOUBLE PRECISION :: dcosth,dph,dps,dg
   DOUBLE PRECISION, DIMENSION(nbox,nbox,nbox) :: odff,odff_nodrx,odffscaled

   !pi = 4.0d0*datan(1.0d0)

   dcosth = 2.0d0/nbox
   dph = pi/nbox
   dps = pi/nbox
   dg = dcosth*dph*dps/(2.0*pi**2)

   DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
      odffscaled(i,j,k) =  dg*odff(i,j,k)/odff_nodrx(i,j,k)
   ENDDO; ENDDO; ENDDO

   END SUBROUTINE SCALE_ODF    
