 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2023, Universita' di Padova, Manuele Faccenda
 !!    All rights reserved.
 !!
 !!    This software package was developed at:
 !!
 !!         Dipartimento di Geoscienze
 !!         Universita' di Padova, Padova         
 !!         via Gradenigo 6,            
 !!         35131 Padova, Italy 
 !!
 !!    project:    ECOMAN
 !!    funded by:  ERC StG 758199 - NEWTON
 !!
 !!    ECOMAN is free software package: you can redistribute it and/or modify
 !!    it under the terms of the GNU General Public License as published
 !!    by the Free Software Foundation, version 3 of the License.
 !!
 !!    ECOMAN is distributed WITHOUT ANY WARRANTY; without even the implied
 !!    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 !!    See the GNU General Public License for more details.
 !!
 !!    You should have received a copy of the GNU General Public License
 !!    along with ECOMAN. If not, see <http://www.gnu.org/licenses/>.
 !!
 !!
 !!    Contact:
 !!        Manuele Faccenda    [manuele.faccenda@unipd.it]
 !!        Brandon VanderBeek  [brandon.vanderbeek@unipd.it]
 !!
 !!
 !!    Main development team:
 !!        Manuele Faccenda    [manuele.faccenda@unipd.it]
 !!        Brandon VanderBeek  [brandon.vanderbeek@unipd.it]
 !!        Albert de Montserrat Navarro
 !!        Jianfeng Yang   
 !!
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------

   PROGRAM DREX_S  

   USE comvar
   USE hdf5    

   IMPLICIT NONE

   INTEGER :: i,j,k,ijk,m,n,nx(3),ti(1),anisnum,nsave

   DOUBLE PRECISION, DIMENSION(3) :: evals,c2
   DOUBLE PRECISION, DIMENSION (3,3) ::LSij,evects,fseacs,acs1,acs2,Rotm,ee,Rotm1,lrot
   DOUBLE PRECISION :: fractvoigt0,fractdisl,phi1,theta,phi2,a0,w3,cf
   DOUBLE PRECISION :: pphi1,ttheta,pphi2
   DOUBLE PRECISION, DIMENSION(6,6) :: Voigt,Reuss,Mixed,Cstilwe
   DOUBLE PRECISION, DIMENSION(1000,6,6) :: Cdem
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_a
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f1_ol,f1_opx
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs_ol,acs_opx
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: ph_ol,th_ol,ps_ol
   ! average orientation of a-axis
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: perc_a,perc_anis,perc_hexa,perc_tetra,perc_ortho,perc_mono,perc_tri
   ! percentage of S wave anisotropy

   CHARACTER (3) :: dt_str3
   CHARACTER (4) :: dt_str4
   CHARACTER (100) :: fname,str

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag

!!! initialization

   CALL init0

   fractvoigt0 = fractvoigt

   anisnum = 1000
   ALLOCATE(phi_a(anisnum,3))
   ALLOCATE(perc_a(anisnum))    ; perc_a     = 0d0 ! norm of hexa / norm tensor (%)
   ALLOCATE(perc_anis(anisnum)) ; perc_anis  = 0d0 ! norm of anis / norm tensor (%)
   ALLOCATE(perc_hexa(anisnum)) ; perc_hexa  = 0d0 ! norm of hexa / norm of anis(%)
   ALLOCATE(perc_tetra(anisnum)); perc_tetra = 0d0 ! norm of tetra/ norm of anis(%)
   ALLOCATE(perc_ortho(anisnum)); perc_ortho = 0d0 ! norm of ortho/ norm of anis(%)
   ALLOCATE(perc_mono(anisnum)) ; perc_mono  = 0d0 ! norm of mono / norm of anis(%)
   ALLOCATE(perc_tri(anisnum))  ; perc_tri   = 0d0 ! norm of tri  / norm of anis(%)
   IF(texmod > 0) THEN
      ALLOCATE(f1_ol(nboxnum),f1_opx(nboxnum))
      ALLOCATE(acs_ol(nboxnum,3,3),acs_opx(nboxnum,3,3))
      ALLOCATE(ph_ol(nbox,nbox,nbox),th_ol(nbox,nbox,nbox),ps_ol(nbox,nbox,nbox))
   END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  LPO  and FSE calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Constant timestep   
   dt = 1d-2/epsnot(1)

   write(*,"(a,f8.2)"),' Timestep = ',dt
   write(*,*)

   fractdisl = fractdislrock(rocktype(1))

   n = 0
   nsave = 0

   DO WHILE (max_strain <= strain_max)

      n = n + 1

!!! FSE and LPO calculation

      IF(strain_max > 0) THEN

         CALL strain(1,1,fractdisl)

         max_strain = max_strain + dt*epsnot(1)!*fractdisl

         timesum = timesum + dt

      END IF

      IF(n==10 .OR. max_strain == strain_max) THEN

      write(*,"(a,f8.2)"),' Finite strain:',max_strain

      nsave = nsave + 1

!!! Cijkl tensor (using Voigt average)
      !CALL stifftenz(1)
      fractvoigt = fractvoigt0 ; CALL tensorscalc(1,mtk0,mpgpa0,Mixed)
      fractvoigt = 1d0 ; CALL tensorscalc(1,mtk0,mpgpa0,Voigt)
      fractvoigt = 0d0 ; CALL tensorscalc(1,mtk0,mpgpa0,Reuss)

!!! Percentage of anisotropy and orientation of axis of hexagonal symmetry
      CALL DECSYM(Mixed,rocktype(1),perc_a(nsave),phi_a(nsave,:),1,perc_anis(nsave),perc_hexa(nsave),perc_tetra(nsave),perc_ortho(nsave),perc_mono(nsave),perc_tri(nsave))

      n = 0

      !!! Write infos in hdf5 format

      !Initialize FORTRAN interface.

      CALL H5open_f (error)
   
      ! Create a new file using default properties.
      IF(max_strain < 9.999) THEN 
         write(dt_str3,'(1f3.1)') max_strain
         fname = trim(output_name)//dt_str3//'.h5'
      END IF
      IF(max_strain>= 9.999 .AND. max_strain < 100) THEN 
         write(dt_str4,'(1f4.1)') max_strain
         fname = trim(output_name)//dt_str4//'.h5'
      END IF

      print *
      print *,' Save to ',trim(fname)
      print *

      CALL H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

      IF(texmod == 0) THEN

      nx(1)=size
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs(:,:,:,1),'acs',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf(1,:),'odf',1)
      IF(rocktype(1) == 1) THEN
         CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens(:,:,:,1),'acs_ens',1)
         CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens(1,:),'odf_ens',1)
      END IF

      ELSE

      CALL ODFDRX(2,0,1,1,f1_ol,f1_opx,acs_ol,acs_opx)

      nx(1)=nboxnum
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ol(:,:,:),'acs',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,f1_ol(:),'odf',1)

      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_opx(:,:,:),'acs_ens',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,f1_opx(:),'odf_ens',1)



      END IF

      nx(1)=3
      nx(2)=3
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)
      nx(1)=6
      nx(2)=6
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Voigt,'Voigt',1)
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Reuss,'Reuss',1)
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Mixed,'Mixed',1)
      !Density
      CALL loadsave_double(0,1,file_id,1,H5T_NATIVE_DOUBLE,rho,'Density',1)
      !Rocktype
      CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,rocktype,'Rocktype',1)

      !Terminate access to the file.
      CALL H5Fclose_f(file_id, error)

      !Close FORTRAN interface.
      CALL H5close_f(error)

      IF(strain_max == 0) max_strain = 1.0

      END IF
      !END OUTPUT

   END DO

   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   print *,'Voigt average'
   write(*,'(6f10.3)') Voigt
   print *
   print *,'Reuss average'
   write(*,'(6f10.3)') Reuss
   print *
   print *,'Mixed Voigt/Ruess average'
   write(*,'(6f10.3)') Mixed
   print *
   print *,'Density'
   write(*,'(1f10.3)') rho   
   print *
   write(*,"(a)"),'--------------------------------------------------------'
   print *

   if(1==0 .AND. texmod > 0) then

      cf = 180.0d0/pi
      OPEN (31,file='ol.out')
      OPEN (32,file='opx.out')
      DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
         ijk = k + (j-1)*nbox + (i-1)*nbox*nbox
         write(31,1002) cf*ph_ol(i,j,k),cf*th_ol(i,j,k),cf*ps_ol(i,j,k),f1_ol(ijk)
         write(32,1002) cf*pha0(i,j,k),cf*tha0(i,j,k),cf*psa0(i,j,k),f1_opx(ijk)
      ENDDO; ENDDO; ENDDO
      CLOSE (31)
      CLOSE (32)

   end if

   ! Save tensor aniisotrorpic components
   !!! Write infos in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)
   
   CALL H5Fcreate_f('anisdec.h5', H5F_ACC_TRUNC_F, file_id, error)

   nx(1)=nsave
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_anis(1:nsave) ,'perc_anis' ,1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_hexa(1:nsave) ,'perc_hexa' ,1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_tetra(1:nsave),'perc_tetra',1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_ortho(1:nsave),'perc_ortho',1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_mono(1:nsave) ,'perc_mono' ,1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,perc_tri(1:nsave)  ,'perc_tri'  ,1)
   CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,nsave,'nsave',1)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  SPO calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF(spomod > 0) THEN
  
      !!! Read input file infos for DEM model and allocate TD database matrices
      CALL initspo

      i = 1       

      !Find principal axes of FSE
      IF(spomod == 1) THEN

         LSij = MATMUL(Fij(:,:,1),TRANSPOSE(Fij(:,:,1)))
         CALL DSYEVQ3(LSij,evects,evals)

         !Order from smallest to largest semiaxis
         DO j = 1,3
            ti = MINLOC(evals) ;fseacs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1))**0.5 ; evals(ti(1))= 1d60
         END DO

         !Save semiaxes orientation
         !1st column: orientaton mininum axis
         !2nd column: orientaton medium  axis
         !3rd column: orientaton maximum axis
         Rotm1 = fseacs
         IF(phi_spo /=0) CALL rot3Dspo(fseacs(:,:),Rotm1(:,:),fseacs(:,2),phi_spo)

         !Rotate tensor parallel to max semiaxis of FSE
         phi1  = atan2(Rotm1(1,1),-Rotm1(2,1)) !Angle of FSE max semiaxis with x1 axis 
         theta =  acos(Rotm1(3,1)) !Angle of FSE min semiaxis with x3 axis
         phi2  = atan2(Rotm1(3,3),Rotm1(3,2))
         CALL rotmatrixZXZ(phi1,theta,phi2,acs2)

         !STILWE (Backus, 1962, JGR)
         CALL stilwe(1,mtk0,mpgpa0,Sav(:,:,i)) 
         
         IF(ptmod == 0) rho(1) = ro_back*(1.0-Vmax) + ro_incl*Vmax

      END IF
    
      IF(spomod > 1) THEN

         ee = e(1,:,:)
         CALL DSYEVQ3(ee,evects,evals)
         
         !Order from smallest to largest semiaxis
         DO j = 1,3
            ti = MINLOC(evals) ;fseacs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1)) ; evals(ti(1))= 1d60
         END DO

         !Save semiaxes orientation
         !1st column: orientaton compressive  axis
         !2nd column: orientaton intermediate axis
         !3rd column: orientaton extensionsal axis
         Rotm1 = fseacs

         !Rotate tensor parallel to compressive (minimum) strain rate axis
         phi1  = atan2(Rotm1(1,3),-Rotm1(2,3)) !Angle of compressive strain rate axis with x1 axis 
         theta =  acos(Rotm1(3,3)) !Angle of extensional strain rate axis with x3 axis
         phi2  = atan2(Rotm1(3,1),Rotm1(3,2))

         a0 = phi_spo
         IF(phi_spo /=0) THEN
            !Rotate velocity gradient tensor parallel to intermediate strain rate axis
            pphi1  = atan2(Rotm1(1,2),-Rotm1(2,2)) 
            ttheta =  acos(Rotm1(3,2))
            pphi2  = 0d0
            CALL rotmatrixZXZ(pphi1,ttheta,pphi2,acs1)
            !L' = R'*L*R
            lrot = MATMUL(TRANSPOSE(acs1),l(1,:,:))
            lrot = MATMUL(lrot,acs1)
            !Vorticity
            w3 = 0.5*(lrot(2,1) - lrot(1,2))
            IF(w3 < 0) a0 = -phi_spo
            !CALL rot3Dspo(fseacs(:,:),Rotm1(:,:),fseacs(:,2),a0)

         END IF

         IF(spomod == 3) Sav(:,:,i) = Mixed

         CALL rotmatrixZXZ(phi1,theta,phi2,acs2)
         CALL tensorrot_aggr(i,TRANSPOSE(acs2)) 
   
         !Now rotate around Y axis if phi_spo different from 0
         IF(phi_spo /= 0) THEN
            CALL rotmatrixZYZ(0d0,a0,0d0,acs1)
            CALL tensorrot_aggr(i,TRANSPOSE(acs1)) !Reverse rotation with transpose rot. matrix
         END IF

         !DEM (Mainprice, 1997, EPSL)
         CALL DEM(1,Vmax,Cdem) 

         !Find DEM elastic tensor according to porosity of the geodynamic model
         j = INT(Vmax/Vstp) + 1
         Sav(:,:,i) = Cdem(j,:,:)

         IF(ptmod == 0 .OR. spomod == 2) rho(1) = ro_back*(1.0-Vmax) + ro_incl*Vmax
         IF(ptmod >  0) rho(1) = rho(1)*(1.0-Vmax)  + ro_incl*Vmax

         !Now rotate around Y axis if phi_spo different from 0
         IF(phi_spo /= 0) THEN
            CALL tensorrot_aggr(i,acs1) !Reverse rotation with transpose rot. matrix
         END IF

      END IF

      !Rotate tensor parallel to:
      !spomod = 1: max semiaxis of FSE 
      !spomod > 1: compressive axis
      CALL tensorrot_aggr(i,acs2) 
    
      Mixed = Sav(:,:,i)

      fname = trim(output_name)//'SPO.h5'
      print *
      print *,'Save to ',trim(fname)
      print *

      !!! Write infos in hdf5 format

      !Initialize FORTRAN interface.

      CALL H5open_f (error)
   
      CALL H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

      nx(1)=3
      nx(2)=3
      IF(spomod == 1) CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)
      IF(spomod >  1) CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,ee,'Eij',1)
      nx(1)=6
      nx(2)=6
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Mixed,'Mixed',1)
      !Density
      CALL loadsave_double(0,1,file_id,1,H5T_NATIVE_DOUBLE,rho,'Density',1)
      !spomod  
      CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,spomod,'spomod',1)

      !Terminate access to the file.
      CALL H5Fclose_f(file_id, error)

      !Close FORTRAN interface.
      CALL H5close_f(error)

      print *,'Mixed Voigt/Ruess average + SPO'
      IF(spomod == 1) print *,'Tensor rotated parallel to FSE max semiaxis'
      IF(spomod >  1) print *,'Tensor rotated parallel to compressive axis'
      write(*,'(6f10.3)') Mixed 
      print *
      write(*,"(a)"),'--------------------------------------------------------'
      print *
 
   END IF

   str='mkdir '//trim(output_name)
   CALL system(str)
   str='mv '//trim(output_name)//'*.h5 '//trim(output_name)
   CALL system(str)
   str='mv anisdec.h5 '//trim(output_name)
   CALL system(str)

   !!! Write infos in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   CALL H5Fcreate_f('fossilfabric.h5', H5F_ACC_TRUNC_F, file_id, error)

   IF(texmod == 0) THEN

      nx(1)=size
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs(:,:,:,1),'acs',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf(1,:),'odf',1)
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens(:,:,:,1),'acs_ens',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens(1,:),'odf_ens',1)
   
   ELSE
   
      nx(1)=nboxnum
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ol(:,:,:),'acs',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,f1_ol(:),'odf',1)
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_opx(:,:,:),'acs_ens',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,f1_opx(:),'odf_ens',1)



   END IF

   nx(1)=3
   nx(2)=3
   CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)

   !spomod and texmod  
   CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,spomod,'spomod',1)
   CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,texmod,'texmod',1)
   CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,size,'size',1)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   1002 format(4(1pe11.3))

   STOP

   END PROGRAM DREX_S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0

   USE comvar
   USE hdf5  

   IMPLICIT NONE

   INTEGER :: gi,j,i1,i3,i,j1,j2,j3,nrot,k ! loop counters
   INTEGER :: iph,ith,ips

   DOUBLE PRECISION :: dph,dcosth,dps,ph,costh,ps,th
   ! matrix of random numbers used to generate initial random LPO

   DOUBLE PRECISION :: phi1,theta,phi2
   ! eulerian angles

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects,ee
   ! eigen values and vectors in jacobi

   DOUBLE PRECISION phi
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phix0,phiy0,phiz0,phiz1
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Rhodum,Vpdum,Vsdum

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag

   ! to fix random distribution
   INTEGER, ALLOCATABLE :: seed(:)

   l = 0d0 ; e = 0d0

   pi = 3.141592653589793238462643383279
   deg2rad = pi/180.0

!!! Name of input files from file input_fabric0.dat

   CALL read_input_file

!!! initial size
   size = size3**3
   nboxnum = size3**3

   Xol = Xol /100.0

!!! strain rate tensor
   e(1,1,1) = l(1,1,1) ; e(1,3,3) = l(1,3,3) ; e(1,2,2) = l(1,2,2)
   e(1,3,1) = (l(1,1,3)+l(1,3,1))/2d0 ; e(1,1,3) = e(1,3,1)
   e(1,1,2) = (l(1,2,1)+l(1,1,2))/2d0 ; e(1,2,1) = e(1,1,2)
   e(1,2,3) = (l(1,3,2)+l(1,2,3))/2d0 ; e(1,3,2) = e(1,2,3)

!! reference strain rate
   ee = e(1,:,:)
   CALL DSYEVQ3(ee,evects,evals)
   epsnot(1) = MAXVAL(ABS(evals))

   write(*,"(a)"),' VELOCITY GRADIENT TENSOR Dij'
   write(*,*)
   write(*,"(3f8.2)"),(l(1,1,:))
   write(*,"(3f8.2)"),(l(1,2,:))
   write(*,"(3f8.2)"),(l(1,3,:))
   write(*,*)
   write(*,"(a)"),' STRAIN RATE tensor'
   write(*,*)
   write(*,"(3f8.2)"),(e(1,1,:))
   write(*,"(3f8.2)"),(e(1,2,:))
   write(*,"(3f8.2)"),(e(1,3,:))
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   max_strain = 0d0
        
!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:) = 1d0 ; Fij(2,2,:) = 1d0 ; Fij(3,3,:) = 1d0

!!! tensor \epsilon_{ijk}

   alt=0d0
   alt(1,2,3) = 1d0 ; alt(2,3,1) = 1d0 ; alt(3,1,2) = 1d0
   alt(1,3,2) = -1d0 ; alt(2,1,3) = -1d0 ; alt(3,2,1) = -1d0

!!! tensor \delta_{ij}

   del=0d0
   del(1,1) = 1d0 ; del(2,2) = 1d0 ; del(3,3) = 1d0

!!! tensors of indices

   ijkl(1,1) = 1 ; ijkl(1,2) = 6 ; ijkl(1,3) = 5
   ijkl(2,1) = 6 ; ijkl(2,2) = 2 ; ijkl(2,3) = 4
   ijkl(3,1) = 5 ; ijkl(3,2) = 4 ; ijkl(3,3) = 3

   l1(1) = 1 ; l1(2) = 2 ; l1(3) = 3
   l1(4) = 2 ; l1(5) = 3 ; l1(6) = 1
   l2(1) = 1 ; l2(2) = 2 ; l2(3) = 3
   l2(4) = 3 ; l2(5) = 1 ; l2(6) = 2

   mandel_scale = 1d0

   DO j1=1,6 ; DO j2=1,6
      IF(j1 .GT. 3 .AND. j2 .GT. 3) THEN
         mandel_scale(j1,j2) = 2
         !mandel_scale(j1,j2) = 4
      ELSE IF(j1 .GT. 3 .OR. j2 .GT. 3) THEN
         mandel_scale(j1,j2) = 2**0.5d0
         !mandel_scale(j1,j2) = 2
      END IF
   END DO ; END DO

!!! Loading stiffness tensors (GPa)

   CALL elastic_database(S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5)

   
   ALLOCATE(acs0(size,3,3))

   nbox = size3
   dph = pi/REAL(nbox)
   dcosth = 2.0d0/REAL(nbox)
   dps = pi/REAL(nbox)

   !CALL rotmatrixZXZ(0d0,pi/2d0,0d0,acsrot)

   !  Grid points at centers of boxes
   DO iph = 1, nbox
      ph = (iph - 0.5)*dph
   DO ith = 1, nbox
      costh = -1.0 + (ith - 0.5)*dcosth
      th = dacos(costh)
   DO ips = 1, nbox
      ps = (ips - 0.5)*dps

      i = ips + (ith-1)*nbox + (iph-1)*nbox*nbox

!!! Direction cosine matrix
!!! acs(k,j) = cosine of the angle between the kth crystallographic axes and the
!jth external axes 

      CALL EULER_TO_DIRCOS(ph,th,ps,acs0(i,:,:))  ! direction cosines g_(ij) of crystal axes

   ENDDO; ENDDO; ENDDO

   !If only 1 crystal, aling crystallographic axes with external reference frame
   IF(size == 1) THEN
      acs0 = 0d0
      acs0(1,1,1) = 1d0
      acs0(1,2,2) = 1d0
      acs0(1,3,3) = 1d0
   END IF

   IF(texmod == 0) THEN

!!! allocation of the dimensions of the arrays

   ALLOCATE(odf(1,size))
   ALLOCATE(odf_ens(1,size))
   ALLOCATE(acs(size,3,3,1),acs_ens(size,3,3,1))

!!! Set random initial LPO and same grain size

   DO i = 1 , 1        

      odf(i,:) = 1d0/REAL(size3**3)
      odf_ens(i,:) = odf(i,:)
      acs(:,:,:,i) = acs0
      acs_ens(:,:,:,i) = acs0

   END DO

   ELSE

   ALLOCATE(pham(1,nbox,nbox,nbox),pha0(nbox,nbox,nbox))
   ALLOCATE(tham(1,nbox,nbox,nbox),tha0(nbox,nbox,nbox))
   ALLOCATE(psam(1,nbox,nbox,nbox),psa0(nbox,nbox,nbox))
   ALLOCATE(fm(1,nbox,nbox,nbox),fm_nodrx(1,nbox,nbox,nbox),f0(nbox,nbox,nbox))
   rocktype(1) = 1

   dph = pi/REAL(nbox)
   dcosth = 2.0/REAL(nbox)
   dps = pi/REAL(nbox)

!  Initial uniform grid in the Euler space; grid points at centers of boxes
!  to avoid singularities
   DO i = 1, nbox; DO j = 1, nbox; DO k = 1, nbox
     pha0(i,j,k) = (i - 0.5)*dph
     costh = -1.0 + (j - 0.5)*dcosth
     tha0(i,j,k) = dacos(costh)
     psa0(i,j,k) = (k - 0.5)*dps
     f0(i,j,k) = 1.d0
   ENDDO; ENDDO; ENDDO

   pham(1,:,:,:) = pha0
   tham(1,:,:,:) = tha0
   psam(1,:,:,:) = psa0
   fm(1,:,:,:) = f0
   fm_nodrx(1,:,:,:) = f0

   END IF

!!! Cartesian coordinates of a spherical object    

   degstp=15 !!! step in degrees
   nxy=360/degstp
   nz=90/degstp + 1
    
   ALLOCATE(phix0(nxy),phiy0(nxy),phiz0(nz),phiz1(nz))
   ALLOCATE(phix(nz,nxy),phiy(nz,nxy),phiz(nz,nxy))

   DO i = 1 , nz 
      phi = degstp*(i-1)*acos(-1.0)/180
      phiz0(i) = cos(phi); 
      phiz1(i) = sin(phi); 
   END DO 

   DO i = 1 , nxy
      phi = degstp*(i-1)*acos(-1.0)/180
      phix0(i) = cos(phi); 
      phiy0(i) = sin(phi); 
   END DO 
  
   DO i = 1 , nz
      DO j = 1 , nxy
      phix(i,j) = phix0(j) * phiz1(i)
      phiy(i,j) = phiy0(j) * phiz1(i)
      phiz(i,j) = phiz0(i)           
      END DO
   END DO

   DEALLOCATE(phix0,phiy0,phiz0,phiz1)

   !!! Load mantle density database made with PERPLE_X
   rho(1) = 3353d0
  
   ALLOCATE(mdb(1))

   IF(ptmod > 0) THEN

      !P-T domain grid
      tknum = 85 ; tkmin = 300; tkstp = 50
      pbnum = 1401 ; pbmin = 0d0; pbstp = 1d3
      ptnum = tknum*pbnum

      ALLOCATE(Rhodum(ptnum),Vpdum(ptnum),Vsdum(ptnum))

      CALL H5open_f (error)

      IF(eosmod == 1) CALL H5Fopen_f("../DATABASES/MMA-EoS/dunite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 2) CALL H5Fopen_f("../DATABASES/MMA-EoS/hartzburgite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 3) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyrolite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 4) CALL H5Fopen_f("../DATABASES/MMA-EoS/morb.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 5) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyroxenite.h5", H5F_ACC_RDONLY_F, file_id, error)

      CALL H5Gopen_f(file_id, "/rock", group_id, error)

      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum,'Rho',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum,'Vp',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum,'Vs',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)

      ALLOCATE(td_rho(1,tknum,pbnum),td_vp(1,tknum,pbnum),td_vs(1,tknum,pbnum))

      DO j = 1 , tknum ; DO i = 1 , pbnum
 
         gi = j + (i-1)*tknum
         td_rho(1,j,i) = Rhodum(gi) !kg/m^3
         td_vp(1,j,i)  = Vpdum(gi)  !km/s
         td_vs(1,j,i)  = Vsdum(gi)  !km/s

      END DO ; END DO

      DEALLOCATE(Rhodum,Vpdum,Vsdum)

      CALL rhopt(1,mtk0,mpgpa0*1d9)

   END IF

   RETURN

   END SUBROUTINE init0
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Rotate tensor parallel to shortes axis of FSE  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE tensorrot_aggr(m,acs1)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,p,q,r,ss,m
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,Cav
   DOUBLE PRECISION, DIMENSION(3,3) :: acs1

   !Angles are commonly defined according to the right-hand rule. Namely, they
   !have positive values when they represent a rotation that appears clockwise
   !when looking in the positive direction of the rotating axis, and negative
   !values when
   !the rotation appears counter-clockwise. 

   Cav = 0d0 ; C0 = 0d0

   !Convert to Cijkl 4th order tensor
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      C0(i,j,k,ll) = Sav(ijkl(i,j),ijkl(k,ll),m)
   END DO ; END DO ; END DO ; END DO

   !Rotate 4th order tensor
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
          Cav(i,j,k,ll) = Cav(i,j,k,ll) + acs1(i,p)*acs1(j,q)*acs1(k,r)*acs1(ll,ss)*C0(p,q,r,ss)
      END DO ; END DO ; END DO ; END DO
   END DO ; END DO ; END DO ; END DO

   !Convert to Voigt notation
   DO i = 1 , 6 ; DO j = 1 , 6
      Sav(i,j,m) = Cav(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   RETURN

   END SUBROUTINE tensorrot_aggr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rhopt(m,mtk,mpa)

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION mtk,mpa,mpgpa,mpb,ee,n,R0,R1,R2,R3
   INTEGER m,n1,n2,db

   IF(mtk<300d0) mtk=300d0
   IF(mpa<1d5) mpa=1d5

   ! Transform pressure in GPa
   mpgpa = mpa / 1d9

   !!! Density

   ! Transform pressure from Pascal to bar
   mpb = mpa / 1d5

   ! ABCD-4Cell Number 
   ee = (mtk-tkmin)/tkstp
   IF(ee < 0) ee=0
   IF(ee > REAL(tknum)) ee=REAL(tknum)
   n=(mpb-pbmin)/pbstp 
   IF(n < 0) n=0
   IF(n > REAL(pbnum)) n=REAL(pbnum)
   n1=FLOOR(ee) + 1
   IF(n1 < 1) n1=1
   IF(n1 > tknum-1) n1=tknum-1
   n2=FLOOR(n)+1
   IF(n2 < 1) n2=1
   IF(n2 > pbnum-1) n2=pbnum-1
   ! Calc normalized distances 
   ee=(ee-REAL(n1-1))
   n=(n-REAL(n2-1))

   ! Choose database
   db = mdb(m)

   ! Ro values
   ! 0 2
   ! 1 3 
   R0=td_rho(db,n1  ,n2  )
   R1=td_rho(db,n1  ,n2+1)
   R2=td_rho(db,n1+1,n2  )
   R3=td_rho(db,n1+1,n2+1)
   
   rho(m)=((R0*(1.0-n)+R1*n)*(1.0-ee)+(R2*(1.0-n)+R3*n)*ee)

   END SUBROUTINE rhopt           

