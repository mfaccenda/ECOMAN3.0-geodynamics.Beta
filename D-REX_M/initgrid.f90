 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2026, Universita' di Padova, Manuele Faccenda
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0(input_dir,output_dir)

   USE comvar
   USE hdf5
   USE omp_lib

   !M3E!!!!!!!!!!!!!!
   use class_DistGrid
   !M3E!!!!!!!!!!!!!!

   IMPLICIT NONE

   INTEGER :: i1,i2,i3,i,j,k,gi,j1,j2,j3,m,ptmod0,ptnum,yy,marknum00,t,tdmax,db ! loop counters
   INTEGER :: iph,ith,ips

   ! eulerian angles
   DOUBLE PRECISION :: phi1,theta,phi2,memnodes,memmark,dummy,lr,cr
   DOUBLE PRECISION :: dph,dcosth,dps,ph,costh,ps,th

   ! matrixes of initial random eulerian angles
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Rhodum,Vpdum,Vsdum

   ! input and output directories
   CHARACTER (len=500) :: input_dir,output_dir

   CHARACTER(500) :: str

   INTEGER(HID_T)  :: file_id, group_id !Handles
   INTEGER     ::   error  ! Error flag

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(kind=double) :: DistTime
   integer           :: errMPI,rankMPI,sizeMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeMPI,errMPI)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !$omp parallel    
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel

   if ( rankMPI .eq. 1 ) then
      print *,' '
      write(*,'(a,i0)') ' Number of threads = ',nt
      print *,' '
   endif

   ALLOCATE(l(nt,3,3),e(nt,3,3),epsnot(nt))

   pi = 3.141592653589793238462643383279d0
   deg2rad = pi/180.0d0

   tau = 1d60

   !Read input file infos
   CALL read_input_file(input_dir,output_dir)

   !Number of files to be processed
   tnum=(Tend-Tinit)/Tstep+1

   !Scale coordinates 
   IF(cartspher>1) THEN

      !Set min/max grid coordinates for Yin-Yang grids
      IF(yinyang == 2) THEN
         x1stp = 270d0/REAL(nx1-1-2)
         x3stp =  90d0/REAL(nx3-1-2)
         x1min = -135d0 - x1stp; x1max = 135d0 + x1stp
         x3min =   45d0 - x3stp; x3max = 135d0 + x3stp
         !mx1min = x1min; mx1max = x1max
         !mx3min = x3min; mx3max = x3max
      END IF

      ! Convert longitude/colatitude to radians
      x1min     = x1min*deg2rad;     x1max     = x1max*deg2rad
      mx1min    = mx1min*deg2rad;    mx1max    = mx1max*deg2rad
      mx1minfab = mx1minfab*deg2rad; mx1maxfab = mx1maxfab*deg2rad
      x3min     = x3min*deg2rad;     x3max     = x3max*deg2rad
      mx3min    = mx3min*deg2rad;    mx3max    = mx3max*deg2rad
      mx3minfab = mx3minfab*deg2rad; mx3maxfab = mx3maxfab*deg2rad
   END IF

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   DistTime = MPI_Wtime()

   call DistGrid(nproc1,nproc2,nproc3,                                                    &
   &             dimensions,cartspher,                                                    &
   &             x1min,x1max,nx1,x1periodic,                                              &
   &             x2min,x2max,nx2,x2periodic,                                              &
   &             x3min,x3max,nx3,x3periodic,                                              &
   &             mx1min,mx1max,mx1stp,                                                    &
   &             mx2min,mx2max,mx2stp,                                                    &
   &             mx3min,mx3max,mx3stp)

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   DistTime = MPI_Wtime() - DistTime

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,f10.5,a)') ' M3E | TIME DISTGRID : ',DistTime,' s'
      write(*,*)
   endif 
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if ( rankMPI .eq. 1 ) then
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
      write(*,'(a)') ' INITIALIZE EULERIAN GRID, AGGREGATES, ETC.'     
      write(*,*)
   endif 

   !Check if domain is defined correctly
   IF(x1max < x1min) THEN ; print *,' Wrong Eulerian domain because x1max < x1min '; stop; END IF
   IF(x2max < x2min) THEN ; print *,' Wrong Eulerian domain because x2max < x2min '; stop; END IF
   IF(x3max < x3min) THEN ; print *,' Wrong Eulerian domain because x3max < x3min '; stop; END IF
   IF(mx1max < mx1min) THEN ; print *,' Wrong Lagrangian domain because mx1max < mx1min '; stop; END IF
   IF(mx2max < mx2min) THEN ; print *,' Wrong Lagrangian domain because mx2max < mx2min '; stop; END IF
   IF(mx3max < mx3min) THEN ; print *,' Wrong Lagrangian domain because mx3max < mx3min '; stop; END IF
   IF(mx1min < x1min) THEN ; print *,' Markers defined outside the domain because mx1min < x1min '; stop; END IF
   IF(mx2min < x2min) THEN ; print *,' Markers defined outside the domain because mx2min < x2min '; stop; END IF
   IF(mx3min < x3min) THEN ; print *,' Markers defined outside the domain because mx3min < x3min '; stop; END IF
   IF(mx1max > x1max) THEN ; print *,' Markers defined outside the domain because mx1max > x1max '; stop; END IF
   IF(mx2max > x2max) THEN ; print *,' Markers defined outside the domain because mx2max > x2max '; stop; END IF
   IF(mx3max > x3max) THEN ; print *,' Markers defined outside the domain because mx3max > x3max '; stop; END IF

   IF(cartspher==2) THEN
      IF(x1max <  0d0) THEN; print *,'Wrong Eulerian domain because Longitude  <    0'; stop; END IF
      IF(x1max > 2*pi) THEN; print *,'Wrong Eulerian domain because Longitude  > 2*pi'; stop; END IF
!     IF(mx1max <  0d0) THEN; print *,'Wrong Lagrangian domain because Longitude  <    0'; stop; END IF ! wrong error
      IF(mx1max > 2*pi) THEN; print *,'Wrong Lagrangian domain because Longitude  > 2*pi'; stop; END IF
      IF(x3max <  0d0) THEN; print *,'Wrong Eulerian domain because Colatitude  <    0'; stop; END IF
      IF(x3max > pi) THEN; print *,'Wrong Eulerian domain because Colatitude  > pi'; stop; END IF
!     IF(mx3max <  0d0) THEN; print *,'Wrong Lagrangian domain because Colatitude  <    0'; stop; END IF ! wrong error
      IF(mx3max > pi) THEN; print *,'Wrong Lagrangian domain because Colatitude  > pi'; stop; END IF
   END IF

   IF(fossilfabric > 0) THEN
      IF(mx1maxfab < mx1minfab) THEN ; print *,' Wrong fossil fabric domain because mx1maxfab < mx1minfab '; stop; END IF
      IF(mx2maxfab < mx2minfab) THEN ; print *,' Wrong fossil fabric domain because mx2maxfab < mx2minfab '; stop; END IF
      IF(mx3maxfab < mx3minfab) THEN ; print *,' Wrong fossil fabric domain because mx3maxfab < mx3minfab '; stop; END IF
      IF(mx1minfab < x1min) THEN ; print *,' Fossil fabric defined outside the domain because mx1minfab < x1min '; stop; END IF
      IF(mx2minfab < x2min) THEN ; print *,' Fossil fabric defined outside the domain because mx2minfab < x2min '; stop; END IF
      IF(mx3minfab < x3min) THEN ; print *,' Fossil fabric defined outside the domain because mx3minfab < x3min '; stop; END IF
      IF(mx1maxfab > x1max) THEN ; print *,' Fossil fabric defined outside the domain because mx1maxfab > x1max '; stop; END IF
      IF(mx2maxfab > x2max) THEN ; print *,' Fossil fabric defined outside the domain because mx2maxfab > x2max '; stop; END IF
      IF(mx3maxfab > x3max) THEN ; print *,' Fossil fabric defined outside the domain because mx3maxfab > x3max '; stop; END IF
   END IF

   !Check if LPO parameters are defined correctly
   DO i=1,5

      IF(Xol(i) < 0 .OR. Xol(i) > 100d0) THEN
         print *,'Wrong volume fraction of main phase for rocktype ',i,', it should be >= 0 and <= 1 '
         stop
      END IF 
      IF(chi(i) < 0 .OR. chi(i) > 1d0) THEN
         print *,'Wrong GBS volume threshold for rocktype ',i,', it should be >= 0 and <= 1 '
         stop
      END IF
      IF(fractdislrock(i) < 0 .OR. fractdislrock(i) > 1d0) THEN
         print *,'Wrong fraction of deformation accommodated by dislocation creep for rocktype ',i,', it should be >= 0 and <= 1 '
         stop
      END IF   
      IF(single_crystal_elastic_db(i,1) < 1 .OR.  single_crystal_elastic_db(i,1) > 20) THEN
         print *,'Wrong single crystal elastic tensor for main phase for rocktype ',i,', it should be > 0 and < 20'
         stop
      END IF   
      IF(single_crystal_elastic_db(i,2) < 1 .OR.  single_crystal_elastic_db(i,2) > 20) THEN
         print *,'Wrong single crystal elastic tensor for minor phase for rocktype ',i,', it should be > 0 and < 20'
         stop
      END IF   
      IF(stressexp(i) < 0) THEN ; print *,'Negative power-law exponent for rocktype ',i,', it should be >= 1'; stop; END IF   
      IF(Mob(i) < 0)       THEN ; print *,'Negative Mob parameter for rocktype ',i,', it should be > 0'; stop; END IF   
      IF(lambda(i) < 0)    THEN ; print *,'Negative lambda parameter for rocktype ',i,', it should be > 0'; stop; END IF    
      DO j=1,12
         IF(tau(i,j) < 1)  THEN ; print *,'nCRSS < 1 for rocktype ',i, ', it should be >= 1'; stop; END IF 
      END DO

   END DO
  
   !Create ouput directory
   if ( rankMPI .eq. 1) then 
      str='mkdir '//trim(output_dir)
      CALL SYSTEM(str)
   endif

!!! initial size
   size = size3**3
   nboxnum = size3**3

!!! Xol given as a percentage
   Xol = Xol/1d2

!!! Initialization of the Eulerian fields: velocity, velocity gradient, FSE, P, T, etc.
   yy = yinyang
   IF(dimensions == 2) nx3 = 1
   nodenum = nx1*nx2*nx3

   ALLOCATE(X1(nx1),X2(nx2),X3(nx3),Ui(yy,3,nx1,nx2,nx3),Ux(yy,3,nx1,nx2,nx3),Dij(yy,3,3,nx1,nx2,nx3))
   memnodes=nx1+nx2+nx3+nodenum*yy*17 !Including V1,V2 in loadsave.f90
   IF(dimensions == 3) memnodes = memnodes + nodenum*yy !Including V3 in loadsave.f90

   IF(ptmod > 0) THEN 
      ALLOCATE(Tk(yy,nx1,nx2,nx3),Pa(yy,nx1,nx2,nx3))
      memnodes = memnodes + nodenum*yy*4 !Including Tk0 and P0 in loadsave.f90
   END IF
   IF(fractdislmod > 0) THEN 
      ALLOCATE(Fd(yy,nx1,nx2,nx3))
      memnodes = memnodes + nodenum*yy*2 !Including Fd0 in loadsave.f90
   END IF

   x1stp = (x1max-x1min)/dble(nx1-1)
   x2stp = (x2max-x2min)/dble(nx2-1)
   x3stp = 0.d0
   if ( nx3 .gt. 1 ) x3stp = (x3max-x3min)/dble(nx3-1)

   !Define Eulerian grid coordinates
   !X
   X1(1) = x1min
   DO i1 = 2, nx1
      X1(i1) = X1(i1-1) + x1stp     
   END DO   
   X1(nx1) = x1max
   
   !Y
   X2(1) = x2min
   DO i2 = 2, nx2
      X2(i2) = X2(i2-1) + x2stp   
   END DO  
   X2(nx2) = x2max
   
   !Z
   X3(1) = x3min
   DO i3 = 2, nx3
      X3(i3) = X3(i3-1) + x3stp   
   END DO  
   X3(nx3) = x3max
   
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,i13)') ' OK'
      write(*,*)
   endif

   ncycmax = 0
   DO t = Tinit , Tend , Tstep 

!!! Load file with velocity components
      CALL load(0,t,input_dir)
      
      ncycmax = ncycmax + numdt

   END DO
   ncycmax = ncycmax + 1 !add 1 cycle to save initial setup

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,i13)') ' TOTAL NUMBER OF CYCLES FOR THE RUN + INITIAL SETUP = ',ncycmax
      write(*,*)
   endif

!!! Initialization of the Lagrangian grid with crystal aggregates

   ! Calculate number of aggregates for array/matrix memory allocation

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)

   marknum = Ranknnmrk

   ! Double number of markers when yinyang grids are active
   if ( cartspher == 2 .and. dimensions == 3 .and. yinyang == 2 ) then
      marknum = marknum * 2
   end if

   !!! Total number of markers
   marknumtot = 0
   call MPI_Allreduce(marknum, marknumtot, 1, MPI_INTEGER, MPI_SUM,  MPI_COMM_WORLD, errMPI)

   write(*,*)
   write(*,'(a,i0,a,i0)') ' Number of aggregates in the domain = ',marknum,' | rank: ',rankMPI

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)

   if ( rankMPI .eq. 1 ) then
      write(*,*)
      write(*,'(a,1i)') ' Total number of aggregates = ',marknumtot       
      write(*,*)
   endif

   ALLOCATE(mx1(marknum),mx2(marknum),mx3(marknum),mYY(marknum),rocktype(marknum),rho(marknum))
   ALLOCATE(mx123(marknum,ncycmax,3),mrtYY(marknum,ncycmax,2),mdb(marknum))
   ALLOCATE(max_strain(marknum))
   ALLOCATE(Fij(3,3,marknum))

   memmark = marknum*(16.125 + ncycmax*2.5d0) !mx123 is REAL (4 bytes), mdb and mrtYY are INT (1 bytes)
   !memmark = marknum*16 !rocktype and mYY are INT (4 bytes) and thus together amount to one 8 bytes array
   !Account for Sav and Xsave
   IF(fsemod == 0) THEN
      memmark = memmark + marknum*(36 + 21)
      IF(ptmod > 0) memmark = memmark + marknum*2 ! mtk, mpgpa 
   END IF

!!! Initialize arrays
   mx1 = 0d0 ; mx2 = 0d0 ; mx3 = 0d0 ; mx123 = 0d0
   mYY = 1 ; rocktype = 1 ; mrtYY = 1   
   rho = 0d0 ; max_strain = 0d0

!!! Set initial thermodynamic database type
   mdb = 1

!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:) = 1d0 ; Fij(2,2,:) = 1d0 ; Fij(3,3,:) = 1d0

   m = 0
  
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF(cartspher == 1) THEN

      marknum00 = marknum
      call GetCoordCart(dimensions,marknum00,mx1,mx2,mx3)

   ELSE

      marknum00 = marknum
      if(yy == 2) marknum00 = marknum / 2
      call GetCoordSphe(dimensions,marknum00,mx1,mx2,mx3)

      IF(dimensions == 3) THEN
         !Add markers for Yang grid
         IF(yy == 2) THEN
           marknum00 = marknum/2
           DO i1=marknum00+1,marknum
              mx1(i1) = mx1(i1-marknum00)
              mx2(i1) = mx2(i1-marknum00)
              mx3(i1) = mx3(i1-marknum00)
              mYY(i1) = 2
           END DO
         END IF
      END IF
   
   END IF
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   IF(fsemod == 0) THEN

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

!!! scaling factor for elastic tensor inversion

   mandel_scale = 1d0

   DO j1=1,6 ; DO j2=1,6
      IF(j1 > 3 .AND. j2 > 3) THEN
         mandel_scale(j1,j2) = 2d0
      ELSE IF(j1 > 3 .OR. j2 > 3) THEN
         mandel_scale(j1,j2) = 2**0.5d0
      END IF
   END DO ; END DO

!!! Loading stiffness tensors (GPa)

   CALL elastic_database(S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5)

!!! initialization of orientations - uniformally random distribution
!!! Rmq cos(theta) used to sample the metric Eulerian space

!!! allocation of the dimensions of the arrays

   ALLOCATE(acs0(size,3,3))

   nbox = size3
   dph = pi/REAL(nbox)
   dcosth = 2.0d0/REAL(nbox)
   dps = pi/REAL(nbox)

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
   
!!! D-REX
   IF(texmod == 0) THEN

   ALLOCATE(acs(size,3,3,marknum),acs_ens(size,3,3,marknum))

   ALLOCATE(odf(marknum,size))
   ALLOCATE(odf_ens(marknum,size))

   dummy=marknum*size
   dummy=dummy*(9 + 9 + 1 + 1) 
   memmark = memmark + dummy           

!!! Set random initial LPO and same grain size

   dph = 1d0/REAL(size3**3)
   DO i = 1 , marknum
   
      odf(i,:) = dph                  
      odf_ens(i,:) = odf(i,:)
      acs(:,:,:,i) = acs0
      acs_ens(:,:,:,i) = acs0
   
   END DO

   ELSE

!!! Define SBFTEX parameters and grid

   nbox = size3
   nboxnum = size3**3.0d0
   size = nboxnum
   ALLOCATE(pham(marknum,nbox,nbox,nbox),pha0(nbox,nbox,nbox))
   ALLOCATE(tham(marknum,nbox,nbox,nbox),tha0(nbox,nbox,nbox))
   ALLOCATE(psam(marknum,nbox,nbox,nbox),psa0(nbox,nbox,nbox))
   ALLOCATE(fm(marknum,nbox,nbox,nbox),fm_nodrx(marknum,nbox,nbox,nbox),f0(nbox,nbox,nbox))
   
   dph = pi/REAL(nbox)
   dcosth = 2.0d0/REAL(nbox)
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

!!! Set random initial LPO and same grain size

   DO i = 1 , marknum
   
      pham(i,:,:,:) = pha0
      tham(i,:,:,:) = tha0
      psam(i,:,:,:) = psa0
      fm(i,:,:,:) = f0
      fm_nodrx(i,:,:,:) = f0
   
   END DO

   END IF   

!!! End initializing crystal orientations

!!! Load mantle density database made with PERPLE_X
   rho = 3353d0
   IF(ptmod > 0) THEN

      !Maximum number of databases
      tdmax = 1

      !P-T domain grid
      tknum = 85 ; tkmin = 300; tkstp = 50
      pbnum = 1401 ; pbmin = 0d0; pbstp = 1d3
      ptnum = tknum*pbnum

      ALLOCATE(td_rho(tdmax,tknum,pbnum),td_vp(tdmax,tknum,pbnum),td_vs(tdmax,tknum,pbnum))
      ALLOCATE(Rhodum(ptnum),Vpdum(ptnum),Vsdum(ptnum))

      DO db = 1, tdmax

      CALL H5open_f (error)

      IF(db == 1) THEN
         IF(eosmod == 1) CALL H5Fopen_f("../DATABASES/MMA-EoS/dunite.h5", H5F_ACC_RDONLY_F, file_id, error)
         IF(eosmod == 2) CALL H5Fopen_f("../DATABASES/MMA-EoS/hartzburgite.h5", H5F_ACC_RDONLY_F, file_id, error)
         IF(eosmod == 3) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyrolite.h5", H5F_ACC_RDONLY_F, file_id, error)
         IF(eosmod == 4) CALL H5Fopen_f("../DATABASES/MMA-EoS/morb.h5", H5F_ACC_RDONLY_F, file_id, error)
         IF(eosmod == 5) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyroxenite.h5", H5F_ACC_RDONLY_F, file_id, error)
      END IF
      !Add below other databases...
      !IF(db == 2) THEN
      !END IF
      !...

      CALL H5Gopen_f(file_id, "/rock", group_id, error)

      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum,'Rho',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum,'Vp',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum,'Vs',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)

      DO j = 1 , tknum ; DO i = 1 , pbnum
 
         gi = j + (i-1)*tknum
         td_rho(db,j,i) = Rhodum(gi) !kg/m^3
         td_vp(db,j,i)  = Vpdum(gi)  !km/s
         td_vs(db,j,i)  = Vsdum(gi)  !km/s

      END DO ; END DO

      END DO

      DEALLOCATE(Rhodum,Vpdum,Vsdum)

   END IF

!!! Set initial rocktype

   !Load initial P-T fields
   IF(ptmod > 1) CALL load(1,Tend,input_dir)
   DO m = 1, marknum

      CALL rocktypecheck(m)

   END DO

   END IF
   !end of IF(fsemod == 0) THEN

   !Delete markers of YANG grid overlapping with YIN grid
   IF(yy == 2) THEN

     DO m=1,marknum
        IF(mYY(m) == 2) THEN
           CALL yin2yang(mx1(m),mx3(m),lr,cr)
           IF(lr>-pi/2.0 .AND. lr<pi/2.0 .AND. cr>0.25*pi .AND. cr<0.75*pi) rocktype(m)= rocktype(m) + 100
           !IF(lr>X1(3) .AND. lr<X1(nx1-2).AND. cr>X3(3) .AND. cr<X3(nx3-2)) rocktype(m)= rocktype(m) + 100
        END IF
     END DO   

   END IF

   !Save initial setup
   DO m=1,marknum
       mx123(m,1,1) = REAL(mx1(m),4)
       mx123(m,1,2) = REAL(mx2(m),4)
       mx123(m,1,3) = REAL(mx3(m),4)
       mrtYY(m,1,1) = INT(rocktype(m),1)
       mrtYY(m,1,2) = INT(mYY(m),1)
   END DO

   if ( rankMPI .eq. 1 ) then
      write(*,*)
      write(*,'(a,1f6.2,a,1f6.2,a)') ' Total memory needed for Eulerian grid is ',memnodes*8*sizeMPI/1d9,' GB, of which ',memnodes*8/1d9,' per CPU node'
  
      write(*,*)
      write(*,'(a,1f6.2,a,1f6.2,a)') ' Total memory needed for crystal aggregates is ',((memmark-1)*8+1)*sizeMPI/1d9,' GB, of which ',((memmark-1)*8+1)/1d9,' per CPU node'
      ! Add 5% to account for other variables declared
      write(*,*)
      write(*,'(a,1f6.2,a,1f6.2,a)') ' Total memory needed is about ',((memnodes+memmark)*8*sizeMPI/1d9)*1.05d0, ' GB, of which ',((memnodes+memmark)*8/1d9)*1.05d0,' per CPU node'
      write(*,*)
   endif

   if ( rankMPI .eq. 1 ) then
      write(*,"(a)") ' INITIALIZE EULERIAN GRID, AGGREGATES, ETC. OK!'
      write(*,*)
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
   endif

   RETURN

   END SUBROUTINE init0
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE set_fossil_fabric 

   USE comvar
   USE omp_lib
   USE hdf5    

   IMPLICIT NONE

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   INTEGER :: m,nx(3),texmod0,size0
   DOUBLE PRECISION :: lr,rr,cr
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf0,odf_ens0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs00,acs_ens00
   DOUBLE PRECISION, DIMENSION(3,3) :: fse0

   ! matrixes of pre-existing crystal volume fraction, orientation and FSE
   INTEGER(HID_T)  :: file_id ! Handles
!  INTEGER(HID_T)  :: dcpl, group_id, dataspace_id, dataset_id, attr_id, memtype ! unused
   INTEGER     ::   error  ! Error flag

   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if ( rankMPI .eq. 1 ) then
      print *,'Impose preexisting fabric'
      print *
   endif

   ALLOCATE(acs00(size,3,3),acs_ens00(size,3,3),odf0(size),odf_ens0(size))

   !!! Read infos from hdf5 file format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   !  Open file using default properties.

   CALL H5Fopen_f('../D-REX_S/fossilfabric.h5', H5F_ACC_RDONLY_F, file_id, error)

   CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,texmod0,'texmod',0)

   IF(texmod == 0 .AND. texmod0 == 1) THEN
      print *,'Fossil fabric for DREX aggregates computed with SBFTEX in D-REX_S'
      STOP
   END IF
 
   IF(texmod == 0) THEN

      CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,size0,'size',0)
      IF(size0 /= size) THEN
         print *,'Fossil fabric has different number of crystals'
         STOP
      END IF

      nx(1)=size
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs00,'acs',0)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf0,'odf',0)
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens00,'acs_ens',0)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens0,'odf_ens',0)

   END IF

   nx(1)=3
   nx(2)=3
   CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,fse0,'Fij',0)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   DO m = 1 , marknum

      lr = mx1(m) ; rr = mx2(m) ; cr = mx3(m)
      IF(mYY(m) == 2) CALL yin2yang(mx1(m),mx3(m),lr,cr)

      IF(lr >= mx1minfab .AND. lr <= mx1maxfab) THEN
      IF(rr >= mx2minfab .AND. rr <= mx2maxfab) THEN
      IF(cr >= mx3minfab .AND. cr <= mx3maxfab) THEN

         IF(fsemod == 0) THEN

         odf(m,:) = odf0
         acs(:,:,:,m) = acs00
         IF(rocktype(m) == 1) THEN
            acs_ens(:,:,:,m) = acs_ens00
            odf_ens(m,:) = odf_ens0
         ELSE
            acs_ens(:,:,:,m) = acs0
            odf_ens(m,:) = 1d0/REAL(size3**3)
         END IF

         END IF

         Fij(:,:,m) = fse0

         IF(cartspher == 2) CALL rottensor(m)

      END IF; END IF; END IF

   END DO

   DEALLOCATE(acs00,acs_ens00,odf0,odf_ens0)

   END SUBROUTINE set_fossil_fabric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rottensor(m)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m,i,j,k,p,q 

   DOUBLE PRECISION :: phi1,theta,phi2,angle0,angle01,lr,cr
   DOUBLE PRECISION, DIMENSION(3) :: a,b
   DOUBLE PRECISION, DIMENSION(3,3) :: acs1,acs2,acs_ens2,fse1
   ! eulerian angles

   !In D-REX_S, the vertical direction is axis 2, which coincide with the
   !Y axis in cartesian coordinates, so remove 90° from colatitude 
   phi1 = mx1(m) - pi*3.0/2.0 ; theta = mx3(m) - pi/2.0 ; phi2 = 0d0
   IF(mYY(m) == 2) THEN
      CALL yin2yang(mx1(m),mx3(m),lr,cr)
      phi1 = lr - pi*3.0/2.0 ; theta = cr - pi/2.0
   END IF

   !Rotate fse with a single rotation 
   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acs1(1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs1(2,1)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
   acs1(3,1)=SIN(phi2)*SIN(theta)
   
   acs1(1,2)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
   acs1(2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
   acs1(3,2)=COS(phi2)*SIN(theta)
   
   acs1(1,3)=SIN(theta)*SIN(phi1)
   acs1(2,3)=-SIN(theta)*COS(phi1)
   acs1(3,3)=COS(theta)

   fse1(:,:) = Fij(:,:,m)
   Fij(:,:,m) = 0d0

   DO j=1,3; DO k=1,3; DO p=1,3; DO q=1,3
      Fij(j,k,m) = Fij(j,k,m) + acs1(j,p)*acs1(k,q)*fse1(p,q)
   END DO; END DO; END DO; END DO

   !Rotate cosine direction matrix with 2 separate rotations
   IF(fsemod == 0) THEN

   !the minus sign is because we need counter-clockwise rotation
   angle0=-phi1/deg2rad
   angle01=-theta/deg2rad
   !First axis of rotation: z-axis
   a=0d0
   a(3)=1;
   !Second axis of rotation: the new x-axis
   acs2=0d0 
   acs2(1,1)=1
   acs2(2,2)=1
   acs2(3,3)=1
   CALL rot3D(acs2,acs1,a,-angle0)
   b(:)=acs1(:,1)
  
   DO i=1,size

      ! First rotate around the z-axis by long
      acs2(:,:) = acs(i,:,:,m)
      acs_ens2(:,:) = acs_ens(i,:,:,m)
      CALL rot3D(acs2,acs(i,:,:,m),a,angle0)
      CALL rot3D(acs_ens2,acs_ens(i,:,:,m),a,angle0)

      ! Second rotate around the new x-axis by colat
      acs2(:,:) = acs(i,:,:,m)
      acs_ens2(:,:) = acs_ens(i,:,:,m)
      CALL rot3D(acs2,acs(i,:,:,m),b,angle01)
      CALL rot3D(acs_ens2,acs_ens(i,:,:,m),b,angle01)
   
   END DO

   END IF

   END SUBROUTINE rottensor

