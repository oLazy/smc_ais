!==============================================================================
!
!                Sequential Reversible Jump MCMC Sampling 
!                   with AIS bridging distributions
!           for relfection coefficient inversion along track
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, November 07 2011
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: September 18 2010
!
!  Sequential particle filter based on Gilks and Berzuini 2001
!  RJMCMC based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
!  Agostinetti Malinverno 2010
!
!==============================================================================

MODULE DATA_TYPE
   IMPLICIT NONE
   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.0D0), SP=KIND(0.0)
   REAL(KIND=RP),   PARAMETER :: PI  = 3.141592653589793238462643383279502884197_RP
END MODULE DATA_TYPE

!==============================================================================
MODULE RJMCMC_COM
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE

!
! General switches
!
   INTEGER(KIND=IB), PARAMETER :: ISETSEED   = 0     ! Fix the random seed 
   INTEGER(KIND=IB), PARAMETER :: ISINPRIOR  = 0     ! Use sine prior with depth
   INTEGER(KIND=IB), PARAMETER :: ISPHER     = 0     ! Spherical refl. coeff. 
   INTEGER(KIND=IB), PARAMETER :: ILINROT    = 0     ! 1 = turn on linear rotation proposal
   INTEGER(KIND=IB), PARAMETER :: ICOOL      = 0     ! Turn cooling in RJMCMC on
   INTEGER(KIND=IB), PARAMETER :: ISLICE     = 0     ! Run slice univariate sampling as MCMC in fixed dimension
   INTEGER(KIND=IB), PARAMETER :: IMAP       = 0     ! WRITE REPLICA AND EXIT
   INTEGER(KIND=IB), PARAMETER :: ISIM       = 0     ! WRITE REPLICA WITH NOISE AND EXIT
   INTEGER(KIND=IB), PARAMETER :: ICOUPLE_CR = 1     ! Constrain discontinuities of rho and v to the same sign
   INTEGER(KIND=IB), PARAMETER :: subsmp     = 2_IB  ! Subsample Sommerfeld integral plane wave part
   INTEGER(KIND=IB), PARAMETER :: NPAVE      = 1_IB  ! # pings to skip

!!
!!  AUV DATA RJMCMC trial ping
!!
   INTEGER(KIND=IB), PARAMETER :: ICOV    = 1         ! 1 = Sample over sigma
                                                      ! 2 = Use sigma from ping ave
                                                      ! 3 = Use sigma from ping ave but scale by factor (free parameter)
   INTEGER(KIND=IB), PARAMETER :: IdB     = 0         ! 1=Carry out computation in dB
   INTEGER(KIND=IB), PARAMETER :: NANG    = 32
   INTEGER(KIND=IB), PARAMETER :: NLMX    = 20
   INTEGER(KIND=IB), PARAMETER :: NPL     = 4         ! Number parameters per layer
   INTEGER(KIND=IB), PARAMETER :: NAVEF   = 3         ! # freq per band
   REAL(KIND=RP),    PARAMETER :: frbw    = 1._RP/15._RP! Frac. bandwidth for freq ave.
   INTEGER(KIND=IB), PARAMETER :: ipingst = 1         ! Ping # to start at
   INTEGER(KIND=IB), PARAMETER :: NPING   = 100       ! Number of pings
   INTEGER(KIND=IB), PARAMETER :: NPART   = 20      ! Number of particles
   INTEGER(KIND=IB), PARAMETER :: NPARTAIS= 10      ! Number of particles to perform AIS on

!   CHARACTER(len=64)                 :: filebasegl         = 'sim_SRJMH'
!   CHARACTER(len=64)                 :: particle_init_file = 'p001_initial_particles.txt'
   CHARACTER(len=64)                 :: filebasegl         = 'AUV_SRJMH'
   CHARACTER(len=64)                 :: particle_init_file = 'p001_initial_particles.txt'
!   CHARACTER(len=64)                 :: particle_init_file = 'p0003_pave_10_1000_2400_particles.txt'
!   CHARACTER(len=64)                 :: particle_init_file = 'p0070_initial_particles.txt'
   CHARACTER(len=64),DIMENSION(NPING):: pingfilebase
   CHARACTER(len=64),DIMENSION(NPING):: particlefile
   INTEGER(KIND=IB)                  :: filebaselen
   INTEGER(KIND=IB)                  :: filebaselengl = 9
   INTEGER(KIND=IB),PARAMETER  :: IGA     = 0    ! Type of averaging (0=intensity, 1=Gaussian, 2=1/3 octave intensity)
   INTEGER(KIND=IB), PARAMETER :: NBAND   = 3
   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1250._RP, 2400._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 4
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP, 2400._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 6
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP, 2400._RP, &
!                                                 2800._RP, 3200._RP/)
!   REAL(KIND=RP),DIMENSION(NBAND),PARAMETER:: sdtrgt = (/0.015_RP, 0.015_RP, 0.015_RP/)
   REAL(KIND=RP),DIMENSION(NBAND),PARAMETER:: sdtrgt = (/0.017_RP, 0.017_RP, 0.017_RP/)
   REAL(KIND=RP) :: FBW = 25.0_RP

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),DIMENSION(NPL):: minlim   = 0._RP
   REAL(KIND=RP),DIMENSION(NPL):: maxlim   = 0._RP
   INTEGER(KIND=IB)            :: kmin     = 0       ! Min number of layers
   INTEGER(KIND=IB)            :: kmax     = 0       ! Max number of layers
   REAL(KIND=RP),PARAMETER     :: hmin     = 0.01_RP ! Min allowed layer thickness
   REAL(KIND=RP),PARAMETER     :: fact     = 3.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),DIMENSION(NPL):: maxpert  = 0._RP
   REAL(KIND=RP),DIMENSION(NPL):: pertsd   = 0._RP 
   REAL(KIND=RP),DIMENSION(NPL):: pertsdsc = (/ 20._RP,50._RP,40._RP,20._RP /)
   REAL(KIND=RP),DIMENSION(NPL):: w_sl     = 0._RP   ! initial slice width
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE:: fr,fstep ! Total frequency array
   REAL(KIND=RP)                 :: vs01,vs02,vsh1,vsh2,numin,numax
   REAL(KIND=RP)                 :: z_t
   REAL(KIND=RP)                 :: cw
   REAL(KIND=RP)                 :: rw
   REAL(KIND=RP)                 :: hmx
   REAL(KIND=RP)                 :: logLmin
   REAL(KIND=RP)                 :: logLmax
   REAL(KIND=RP)                 :: logLtrgt,logLtrgtsum
   INTEGER(KIND=IB)              :: kmintmp
   INTEGER(KIND=IB)              :: kmaxtmp
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(len=64) :: repfile
!!
!!  Standard deviation prior variables:
!!
   REAL(KIND=RP), DIMENSION(NBAND):: minlimsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NBAND):: maxlimsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NBAND):: maxpertsd = 0._RP
   REAL(KIND=RP), DIMENSION(NBAND):: pertsdsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NBAND):: pertsdsdsc= 18._RP
!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX = (NLMX * NPL) + (NPL-1)
   INTEGER(KIND=IB)           :: NVV        = ((NLMX*NPL)+NPL-1)*((NLMX*NPL)+NPL-1)*NLMX
   INTEGER(KIND=IB)           :: NSDEVM     = ((NLMX*NPL)+NPL-1)*NLMX
   INTEGER(KIND=IB)           :: len_snd,len_snd2
   INTEGER(KIND=IB)           :: ioutside = 0
   INTEGER(KIND=IB)           :: ireject = 0, iaccept = 0
   INTEGER(KIND=IB)           :: ireject_bd = 0, iaccept_bd = 0
   INTEGER(KIND=IB)           :: irejectlf = 0, iacceptlf = 0
   INTEGER(KIND=IB)           :: i_bd         ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_zpert      ! Z perturb track (0=nothing, 1=birth-death, 2=perturb 1 z)
   INTEGER(KIND=IB)           :: i_sdpert = 0 ! if sigma is perturbed, don't compute forward model
   REAL(KIND=RP)              :: acc_ratio

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iskipais = 0       ! Skip AIS if initial MCMC and resampling successful

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: MCMC_INI   = 2e1_IB ! Small # jrMCMC steps to try before AIS
   INTEGER(KIND=IB),PARAMETER    :: MCMC_BI    = 2e1_IB ! # balance steps in traget region burn-in (exits if target reached)
   INTEGER(KIND=IB),PARAMETER    :: MCMC_STEPS = 5e2_IB ! # balance steps in traget region (always does these)
   INTEGER(KIND=IB)              :: NreAIS     = 2_IB   ! # resampling steps along trajectory
   INTEGER(KIND=IB),PARAMETER    :: npercnt    = 10_IB  ! % particles that must be in logL target region
   INTEGER(KIND=IB),PARAMETER    :: NCOOLTRAJ  = 2_IB   ! # cooling trajectories to try before giving up
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 4_IB   ! Misc parameters in sample (for bookeeping)
   INTEGER(KIND=IB),PARAMETER    :: NDM        = 60     ! No. steps in lin rot est
!!
!! Annealing burn-in parameters (sets beta schedule)
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: beta                 ! Inverse Temprerature
   INTEGER(KIND=IB),PARAMETER            :: NTEMP    = 2e2       ! Number of beta values
   INTEGER(KIND=IB),PARAMETER            :: NBAL     = 3e0       ! Number of balanching steps at each T
   REAL(KIND=RP)                         :: beta1gl  = 0.01_RP   ! Global inverse T to define annealing schedule (start)
   REAL(KIND=RP)                         :: beta4gl  = 1.00_RP   ! Global inverse T to define annealing schedule (end for 3 exp legs)

!!
!!  Rotation matrices 
!!
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,(NLMX*NPL)+NPL-1,NLMX):: Cov0   ! Linear parameter cov mat
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,(NLMX*NPL)+NPL-1,NLMX):: VV     ! Linear rotation
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,NLMX)                 :: sdevm  ! Std dev for perturbations
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,NLMX)                 :: dmbest ! Best step-size for derivative est
   INTEGER(KIND=IB),DIMENSION(NLMX)                               :: idmbest! is best step-size already computed in this dimension?

!!
!!  Structures for objects and data 
!!
   TYPE :: objstruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1):: par     ! Forward parameters
      REAL(KIND=RP),DIMENSION(NLMX)         :: z          ! Depth partition
      REAL(KIND=RP),DIMENSION(NLMX)         :: h          ! Depth partition
      REAL(KIND=RP),DIMENSION(NBAND)        :: sdpar      ! Std dev parameters
      REAL(KIND=RP),DIMENSION(NPL-1)        :: g          ! Acoustic parameters birth-death layer
      REAL(KIND=RP),DIMENSION(NPL-1)        :: gp         ! Acoustic parameters birth-death layer perturbed
      REAL(KIND=RP),DIMENSION(NPL-1,NPL-1)  :: Chat,Chati ! Covariance matrix for perturbing one BD layer
      REAL(KIND=RP)                         :: detChat
      INTEGER(KIND=IB)                      :: k          ! Layer dimension
      INTEGER(KIND=IB)                      :: NFP        ! Number forward parameters
      INTEGER(KIND=IB)                      :: nfail      ! Counter to keep track of how often particle failed
      REAL(KIND=RP)                         :: logL       ! log likelihood
      REAL(KIND=RP)                         :: logP       ! log likelihood
      REAL(KIND=RP)                         :: logwt      ! log wt (likelihood ratio when moving from one ping to next)
      REAL(KIND=RP)                         :: logPr      ! log Prior probability ratio
      REAL(KIND=RP)                         :: lognorm    ! Data covariance matrices
      REAL(KIND=RP),DIMENSION(NBAND,NANG)   :: res = 0._RP! Data residuals for current object
   END TYPE objstruc

   TYPE :: datastruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Robs   = 0._RP ! Observed data
      REAL(KIND=RP),DIMENSION(NANG)           :: angobs = 0._RP ! Observed angles
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rrep   = 0._RP ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: res    = 0._RP ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: sigma  = 0._RP ! Index for bad points/data gaps
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cdi    = 0._RP ! Inverse data covariance matrices
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cd     = 0._RP ! Data covariance matrices
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: CdCi   = 0._RP ! Inverse data covariance matrices
      REAL(KIND=RP),DIMENSION(NBAND)          :: lognorm=0._RP  ! Data lognorm
      REAL(KIND=RP),DIMENSION(NBAND)          :: logdet =0._RP  ! Data log-determinant
      INTEGER(KIND=IB), DIMENSION(NBAND)      :: NDPF           ! # data per freq
      INTEGER(KIND=IB)                        :: NANG   = NANG  ! # data/angles (struc copy)
      INTEGER(KIND=IB)                        :: NBAND  = NBAND ! # freq bands (struc copy)
      INTEGER(KIND=IB),DIMENSION(NBAND)       :: NDAT   = 0     ! Observed angles
   END TYPE datastruc

   REAL(KIND=RP),DIMENSION(NBAND):: sdev = 0._RP                ! Standard devs
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: icount

!!
!!  Global variables for spherical reflection coeff computation
!!
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE    :: N,Nr
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE       :: drTh
   COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE    :: diTh
   COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE  :: RPs
   TYPE :: btRstruc
      SEQUENCE
      COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: btR2
      COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE   :: rTh2
   END TYPE btRstruc
   TYPE(btRstruc),DIMENSION(NAVEF*NBAND)       :: Sarg
!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ncount,ncount2,ierr
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,tag,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   INTEGER(KIND=IB)            :: isize1,isize2,isize3

   INTERFACE
      FUNCTION RANDPERM(num)
         USE data_type, ONLY : IB
         IMPLICIT NONE
         INTEGER(KIND=IB), INTENT(IN) :: num
         INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
      END FUNCTION RANDPERM
   END INTERFACE

END MODULE RJMCMC_COM
  
!=======================================================================

PROGRAM  SRJMCMC_PLANE

!=======================================================================
USE MPI
USE RJMCMC_COM
USE NR
USE M_VALMED
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: iping,iiping,ipart,ipartsnd,ipartrcv,iloop,ithread,idest,imcmc
INTEGER(KIND=IB)  :: idone,icountmcmc,ithin,ifail,ik,iipart,itrgtcnt
INTEGER(KIND=IB)  :: isource,NBAL1,nfailtmp,irealloc,iais,ireAIS,ifr,ipar
INTEGER(KIND=IB)  :: NTEMP1,NTEMPtmp
INTEGER(KIND=IB),DIMENSION(NPART)     :: idxfail,idxpart
INTEGER(KIND=IB),DIMENSION(NLMX)      :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX):: kidx

TYPE(objstruc)                             :: obj           ! Objects in RJMH chain
TYPE(objstruc)                             :: objmax        ! Object with max logL (memory for cooling)
TYPE(objstruc)                             :: objnew        ! Objects in RJMH chain
TYPE(objstruc),DIMENSION(NPART)            :: particles     ! All particles (current ping)
TYPE(objstruc),DIMENSION(NPART)            :: particles_old ! All particles (current ping)
TYPE(objstruc),DIMENSION(NPART)            :: particles_new ! All particles (current ping)
TYPE(datastruc),DIMENSION(NPING)           :: dat           ! Data

REAL(KIND=RP)                              :: ran_uni    ! Uniform random number
REAL(KIND=RP)                              :: beta1,beta4! Inverse Ts to pass to MAKE_BETA
REAL(KIND=RP)                              :: logLG      ! Global best log likelihood
REAL(KIND=RP)                              :: LOGPLUS    ! Function: Addition carried out in log space

REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sample_p01 ! PPD for ping 01 (to get initial particles
REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: betatmp    ! PPD for ping 01 (to get initial particles
REAL(KIND=RP),DIMENSION(NPART,NBAND)       :: sdtmp
REAL(KIND=RP),DIMENSION(NBAND)             :: sdmead,logLtrgttmp

!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

REAL(KIND=RP)               :: tstart, tend              ! Overall time 
REAL(KIND=RP)               :: tiping1, tiping2          ! Overall time 
REAL(KIND=RP)               :: tstart2, tend2            ! Time for one forward model computation
REAL(KIND=RP)               :: tstartsnd, tendsnd        ! Communication time
REAL(KIND=RP)               :: tstartcmp, tendcmp, tcmp  ! Forward computation time

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )

!!
!! All threads READ ALL DATA for all pings
!!
iiping = ipingst
DO iping = 1,NPING
   WRITE(pingfilebase(iping),'(a,i3.3,a)')'p',iiping,'_1000_2400'
   filebaselen = 14
   PRINT*, pingfilebase(iping),particlefile(iping)
!   WRITE(pingfilebase(iping),'(a,i4.4,a)')'p',iiping,'_pave_10_1000_2400'
!   filebaselen = 23
   !! Read data for new ping
   CALL READ_DATA(dat(iping),pingfilebase(iping),particlefile(iping))
   iiping = iiping+NPAVE
ENDDO

!!
!!  Prior bounds
!!
!! AUV SIMULATION:
minlim = (/ hmin, 1450._RP, 1.20_RP, 0.001_RP /)
maxlim = (/ hmx,  1700._RP, 2.10_RP, 1.000_RP /)
vs01 = 1475._RP
vs02 = 1560._RP
vsh1 = 1510._RP
vsh2 = 1610._RP
numin= 0.4_RP
numax= 0.8_RP

kmin = 1
kmax = NLMX
kmintmp = kmin
kmaxtmp = kmax
maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc
w_sl = maxpert/10._RP

IF(ICOV == 1)THEN
   !! Set prior and proposal scaling for AR model:
   minlimsd   = 0.00_RP
   maxlimsd   = 0.1_RP
   pertsdsdsc = 30._RP
   maxpertsd  = maxlimsd-minlimsd
   pertsdsd   = maxpertsd/pertsdsdsc
ENDIF

!!
!! Print some stats to the screen
!!
IF(rank == src)THEN
   WRITE(6,211) 'ICOV                 :   ',ICOV
   IF(ICOV == 1)THEN
      WRITE(6,216) '...sampling over standard deviations.                   '
   ENDIF
   IF(ISLICE == 1)THEN
      WRITE(6,216) '...using slice sampling as MCMC step in fixed dimension.'
   ENDIF
  WRITE(6,211) 'ILINROT              : ',ILINROT
  WRITE(6,211) 'NDM                  : ',NDM
  WRITE(6,211) 'ISPHER               : ',ISPHER
  WRITE(6,211) 'Number of angles     : ',NANG
  WRITE(6,211) 'Number of frequencies: ',NBAND
  WRITE(6,211) 'Number of f ave      : ',NAVEF
  WRITE(6,211) 'IGA                  : ',IGA
  WRITE(6,211) 'NFPMX                : ',NFPMX
  WRITE(6,214) 'ang(min)             : ',dat(1)%angobs(1)
  WRITE(6,214) 'ang(max)             : ',dat(1)%angobs(NANG)
  WRITE(6,211) 'kmin                 : ',kmin
  WRITE(6,211) 'kmax                 : ',kmax
  WRITE(6,211) 'subsmp               : ',subsmp
  WRITE(6,214) 'hmin                 : ',hmin
  CALL FLUSH(6)
  WRITE(6,*) ''
  WRITE(6,*) 'minlim:  '
  WRITE(6,215) minlim
  WRITE(6,*) 'maxlim:  '
  WRITE(6,215) maxlim
  WRITE(6,*) 'pertsdsc:'
  WRITE(6,215) pertsdsc
  WRITE(6,*) 'maxlim sigma:  '
  WRITE(6,215) maxlimsd
  WRITE(6,*) ''
  WRITE(6,*) 'Done reading data.'
  WRITE(6,*) ''
  CALL FLUSH(6)
  IF(IGA == 0)THEN
    WRITE(*,*)' !!!   NARROW BAND INTENSITY FREQ AVERAGING   !!!'
    WRITE(*,*)' FBW = ',FBW
  ELSEIF(IGA == 1)THEN
    WRITE(*,*)' !!!   GAUSSIAN FREQ AVERAGING   !!!'
  ELSEIF(IGA == 2)THEN
    WRITE(*,*)' !!!   1/3 octave INTENSITY FREQ AVERAGING   !!!'
  ELSEIF(IGA == 3)THEN
    WRITE(*,*)' !!!   FRACTIONAL BANDWIDTH INTENSITY FREQ AVERAGING   !!!'
    WRITE(*,*)' frbw = ',frbw
  ENDIF
  WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''
  WRITE(6,208),'Done reading data from files ',pingfilebase(1),' to ',pingfilebase(NPING),'.'
ENDIF

!!
!! Read PPD of first ping
!!
CALL READ_PPD_P01(particles,particle_init_file)
obj = particles(1)

!! 1000 & 2400
!obj%par   = 0._RP
!obj%k   = 4
!obj%NFP = (obj%k * NPL) + (NPL-1)
!obj%par(1:obj%NFP) = (/ 0.5_RP, 1480._RP, 1.25_RP, 0.010_RP,&
!                        1.2_RP, 1530._RP, 1.45_RP, 0.085_RP,&
!                        2.1_RP, 1580._RP, 1.65_RP, 0.300_RP,&
!                        2.5_RP, 1580._RP, 1.90_RP, 0.550_RP,&
!                                1580._RP, 1.75_RP, 0.800_RP /)
!IF(ICOV == 1) obj%sdpar = 0.03_RP
CALL CALCH(obj)

CALL WRITE_LOG()      !! Write sampling parameters to log
CALL PARALLEL_SEED()  !! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)
CALL INIT_BD(obj)     !! Init birth-death covariance matrix

ALLOCATE( icount(NTHREAD) )
icount = 0

tstart = MPI_WTIME()

CALL INIT_FREQ(dat(1))  !! Compute freq for freq-averaging, Initialize Sommerfeld integral for Spherical R

IF(IMAP == 1)THEN
   IF(rank == src)WRITE(6,*) 'IMAP activated, exiting after computing replica for MAP.'
   CALL COMPUTE_MAP(obj,dat(1))
   STOP
ELSE
   IF(rank == src)WRITE(6,*) 'logL = ',obj%logL
   IF(rank == src)WRITE(6,*) 'Starting model:'
   IF(rank == src)CALL PRINTPAR(obj)
   tstart2 = MPI_WTIME()
   CALL LOGLHOOD(obj,dat(1))
   icount(rank+1) = icount(rank+1) + 1
   tend2 = MPI_WTIME()
   IF(rank == src)WRITE(6,*) 'logL = ',obj%logL
   IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
ENDIF
CALL FLUSH(6)

!!
!! Initialize rotation matrices
!!
VV    = 0._RP
DO ik = kmin,kmax
DO ipar = 1,NFPMX
   VV(ipar,ipar,ik) = 1._RP
ENDDO
ENDDO

!------------------------------------------------------------------------
!
!       ************ Sequential RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
IF(rank == src)THEN
   WRITE(6,*) 'Starting Sequential RJMCMC sampling...'
   OPEN(UNIT=44,FILE=logfile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
!   OPEN(UNIT=77,FILE='resampling.txt',FORM='formatted',STATUS='REPLACE', &
!   ACTION='WRITE',RECL=1024)
ENDIF
logLtrgt = 0._RP
ifail = 0

iiping = ipingst+NPAVE
DO iping = 2,NPING
   !! Make cooling schedule
   NTEMP1 = NTEMP
   NBAL1  = NBAL
   ALLOCATE( beta(NTEMP1),betatmp(NTEMP1) )
   beta1 = beta1gl
   beta4 = beta4gl
   CALL MAKE_BETA(beta1,beta4,NTEMP1)
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   irealloc = 0
   betatmp = beta

IF(rank==src)THEN
   tiping1 = MPI_WTIME()
   iais = 1
   WRITE(6,*)   ''
   WRITE(6,*)   ''
   WRITE(6,207) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(6,207) '~~                                                                                                                 ~~'
   WRITE(6,206) '~~                                       Working on ping',iiping,'                                                          ~~'
   WRITE(6,207) '~~                                                                                                                 ~~'
   WRITE(6,207) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

   !!
   !! Randomize order of particles since NPARTAIS << NPART
   !!
   idxpart = RANDPERM(NPART)
   DO ipart = 1,NPART
      particles_old(ipart) = particles(idxpart(ipart))
   ENDDO
   particles = particles_old

   !!
   !! Estimate target logL region
   !!
   logLmin = HUGE(1._RP)
   logLmax = -HUGE(1._RP)
   DO ipart = 1,NPART
      CALL LOGLHOOD(particles(ipart),dat(iping))
      IF(particles_old(ipart)%logL > logLmax)logLmax = particles_old(ipart)%logL
      IF(particles_old(ipart)%logL < logLmin)logLmin = particles_old(ipart)%logL
      DO ifr = 1,NBAND
         sdtmp(ipart,ifr) = particles_old(ipart)%sdpar(ifr)
      ENDDO
   ENDDO
   DO ifr = 1,NBAND
      sdmead(ifr) = VALMED(sdtmp(:,ifr))
!      logLtrgttmp(ifr) = -(REAL(NANG,RP)/2._RP)*LOG(2._RP*PI) -(REAL(NANG,RP)*LOG(sdmead(ifr))+REAL(NANG,RP)/2._RP)
      logLtrgttmp(ifr) = -(REAL(NANG,RP)/2._RP)*LOG(2._RP*PI) -(REAL(NANG,RP)*LOG(sdtrgt(ifr))+REAL(NANG,RP)/2._RP)
   ENDDO
   logLtrgt = SUM(logLtrgttmp)
   WRITE(*,213) 'median sigma(ifr) = ',sdmead
   WRITE(*,213) 'logLtrgt(ifr)     = ',logLtrgttmp
   WRITE(*,213) 'logLtrgt          = ',logLtrgt
ENDIF !! MPI IF

!!
!! Compute linear rotation for each dimension (using best fit particle), then BCAST
!!
IF(ILINROT == 1)THEN
  IF(rank == src)CALL COMPUTE_VV(particles,dat(iping))
  PRINT*,rank,'START SND_VV 1'
  CALL SND_VV()
  PRINT*,rank,'DONE SND_VV 1'
ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!                                                                           !!
!! DO SOME MCMC WITH NEW DATA, THEN RESAMPLE                                 !!
!!                                                                           !!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
   WRITE(6,*) 'rjMCMC:'
   WRITE(6,203) '   iping,  iloop,particle,idone,        logL,  k,iacc_bd, acc_ratio, time(prtcle),source, time(comp),ifail,nfail'
   WRITE(6,203) '----------------------------------------------------------------------------------------------------------------'
   ipartsnd = 0
   idxfail = 0
   DO iloop = 1,MIN(NTHREAD-1,NPARTAIS)
      ipartsnd = ipartsnd + 1
      idest = ipartsnd
      idone = 0
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
   ENDDO ! PARTICLE LOOP
   iloop = 0
   DO 
      iloop = iloop + 1
      tstartsnd = MPI_WTIME()
      CALL RCV_PARTICLE(objnew,isource,ipartrcv,ifail,tcmp)
      tendsnd = MPI_WTIME()
      particles(ipartrcv) = objnew
      idest = isource

      IF(MOD(iloop,100)==0)THEN
         WRITE(6,201) iiping,iloop,ipartrcv,ipartsnd,particles(ipartrcv)%logL,particles(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
         WRITE(44,201) iiping,iloop,ipartrcv,ipartsnd,particles(ipartrcv)%logL,particles(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
      ENDIF
      IF(ipartsnd < NPARTAIS)THEN
         tendsnd = MPI_WTIME()
         ipartsnd = ipartsnd + 1
         idone = 0
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ELSE
         idone = 1
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ENDIF

      IF(iloop == NPARTAIS)EXIT
   ENDDO ! PARTICLE LOOP

   kmintmp = NLMX
   kmaxtmp = 0
   DO ipart = 1,NPART
      !! Compute likelihoods for new ping for all particles
!      CALL LOGLHOOD(particles(ipart),dat(iping))
      !! Compute particles log-weights (ratio of likelihoods or difference in logL)
!      particles(ipart)%logwt = particles(ipart)%logL-particles_old(ipart)%logL
      !! THIS IS WHAT YARDIM DOES I THINK; SO ONLY LIKELIHOOD...
      particles(ipart)%logwt = particles(ipart)%logL
      kmintmp = MIN(kmintmp,particles(ipart)%k)
      kmaxtmp = MAX(kmaxtmp,particles(ipart)%k)
   ENDDO ! PARTICLE LOOP
   IF(kmintmp > 1)    kmintmp = kmintmp - 1
   IF(kmaxtmp < NLMX) kmaxtmp = kmaxtmp + 1

   WRITE(6,*) 'kmin for ping',iiping,'= ',kmintmp
   WRITE(6,*) 'kmax for ping',iiping,'= ',kmaxtmp
   !!
   !! Resample new set of particles according to weight
   !! resampled ones are particles_new
   !!
   CALL RESAMPLE(particles,particles_new,kcount,kidx)
   particles = particles_new

   !!
   !! Check if at least some particles are in target region. If not, repeat with AIS three times as many interp. dists.
   !!
   itrgtcnt = 0
   logLmin = HUGE(1._RP)
   logLmax = -HUGE(1._RP)
   DO ipart = 1,NPART
      IF(particles(ipart)%logL > logLtrgt)itrgtcnt = itrgtcnt + 1
      IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
      IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
   ENDDO
   irealloc = 0
   IF((itrgtcnt < npercnt*FLOOR(REAL(NPART,RP)/100._RP)).AND.(iais < NCOOLTRAJ))THEN
     iskipais = 0
     WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
     WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
     WRITE(6,*) itrgtcnt,'particles are in target region. Proceeding with AIS.'
   ELSE
     iskipais = 1
     WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
     WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
     WRITE(6,*) itrgtcnt,'particles are in target region. Skipping AIS.'
   ENDIF
ELSE
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!    INITIAL rjMCMC SLAVE PART
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   DO
      ifail = 0 
      idone = 0
      !! Reset acceptance counters:
      iaccept_bd     = 0
      ireject_bd     = 0
      i_bd           = 0
      i_zpert        = 0
      tstartcmp = MPI_WTIME()
      !! Get new particle from master
      CALL SND_PARTICLE(obj,ipart,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)

      IF(idone == 1) EXIT

      DO imcmc = 1,MCMC_INI
         CALL EXPLORE_MH(obj,dat(iping),1._RP)
      ENDDO
      tendcmp = MPI_WTIME()
   
      tcmp = tendcmp-tstartcmp
      !! Get updated particle back to master
      CALL RCV_PARTICLE(obj,isource,ipart,ifail,tcmp)
   ENDDO
ENDIF !! MPI IF
CALL MPI_BCAST( iskipais,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )

IF(iskipais == 0)THEN
!!
!! Compute linear rotation for each dimension (using best fit particle), then BCAST
!!
IF(ILINROT == 1)THEN
  IF(rank == src)CALL COMPUTE_VV2(particles,dat(iping))
  PRINT*,rank,'START SND_VV 2'
  CALL SND_VV()
  PRINT*,rank,'DONE SND_VV 2'
ENDIF
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Compute AIS weights
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
999 CONTINUE
   WRITE(6,*) 'AIS:'
   WRITE(6,210) '   iping,  iloop,particle,    logL new,    logL old,   logwt new,   logwt old,source, time(comp)'
   WRITE(6,210) '------------------------------------------------------------------------------------------------'

   DO ireAIS = 1,NreAIS

      NTEMPtmp = (NTEMP1/NreAIS)
      irealloc = 1
      DEALLOCATE(betatmp)
      ALLOCATE(betatmp(NTEMPtmp))
      betatmp  = beta((ireAIS-1)*NTEMPtmp+1:ireAIS*NTEMPtmp)
      ipartsnd = 0
      idxfail = 0
      DO ithread = 1,MIN(NTHREAD-1,NPARTAIS)
         ipartsnd = ipartsnd + 1
         idest = ipartsnd
         idone = 0
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMPtmp,NBAL1,irealloc)
      ENDDO ! PARTICLE LOOP
      iloop = 0
      ifail = 0
      DO 
         iloop = iloop + 1
         CALL RCV_PARTICLE(objnew,isource,ipartrcv,ifail,tcmp)
         particles(ipartrcv) = objnew
         idest = isource
         IF(MOD(iloop,100)==0)THEN
            WRITE(6,209) iiping,iloop,ipartrcv,particles(ipartrcv)%logL,particles_old(ipartrcv)%logL,&
                         particles(ipartrcv)%logwt,particles_old(ipartrcv)%logwt,isource,tcmp
            CALL FLUSH(6)
         ENDIF
         IF(ipartsnd < NPARTAIS)THEN
            tendsnd = MPI_WTIME()
            ipartsnd = ipartsnd + 1
            idone = 0
            CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMPtmp,NBAL1,irealloc)
         ENDIF

         IF(iloop == NPARTAIS)EXIT
      ENDDO ! PARTICLE LOOP
      !!
      !! Resample step after ireAIS temperature block
      !!
      CALL RESAMPLE_AIS(particles,particles_new,particles_old,kcount,kidx)

      particles = particles_new
      logLmin = HUGE(1._RP)
      logLmax = -HUGE(1._RP)
      DO ipart = 1,NPART
         IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
         IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
      ENDDO
      WRITE(*,212)'Resampled',ireAIS,' times',logLmin,logLmax

   ENDDO ! AIS resample loop
   !!
   !! Leapfrog a bit
   !!
!   CALL EXPLORE_LEAPFROG(particles,dat)
!   PRINT*,'leapfrog accepted:',iacceptlf

   !!
   !! Check if at least some particles are in target region. If not, repeat with AIS three times as many dists.
   !!
   itrgtcnt = 0
   logLmin = HUGE(1._RP)
   logLmax = -HUGE(1._RP)
   DO ipart = 1,NPART
      IF(particles(ipart)%logL > logLtrgt)itrgtcnt = itrgtcnt + 1
      IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
      IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
   ENDDO
   irealloc = 0
   IF((itrgtcnt < npercnt*FLOOR(REAL(NPART,RP)/100._RP)).AND.(iais < NCOOLTRAJ))THEN
      iais = iais + 1_IB
      DEALLOCATE(beta)
      NTEMP1 = 4_IB*NTEMP1
!      NBAL1 = FLOOR(1.25_RP*REAL(NBAL1,RP))
!      NBAL1 = NBAL1+1_IB
      ALLOCATE(beta(NTEMP1))
      beta1 = beta1/4._RP
      beta4 = beta4gl
      CALL MAKE_BETA(beta1,beta4,NTEMP1)
      irealloc = 1
      WRITE(6,*) 'Repeating AIS for ping',iiping
      WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
      WRITE(6,*) 'NTEMP1 = ',NTEMP1
      WRITE(6,*) 'NBAL1 = ',NBAL1
      WRITE(6,*) 'beta1 = ',beta1
      WRITE(6,*) 'beta4 = ',beta4
      WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
      GOTO 999
   ELSE
      WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
      WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
      WRITE(6,*) itrgtcnt,'particles are in target region. Proceeding to rjMCMC.'
   ENDIF
   CALL FLUSH(6)
   !!
   !! Send all slaveѕ the done if AIS was success.
   !!
   DO ithread = 1,MIN(NTHREAD-1,NPART)
      idone = 1
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,ithread,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
   ENDDO
   !!
   !! Resample step
   !!
!   CALL RESAMPLE_AIS(particles,particles_new,particles_old,kcount,kidx)

ELSE
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!    AIS SLAVE PART
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   DO
      ifail = 0
      idone = 0
      tstartcmp = MPI_WTIME()
      CALL SND_PARTICLE(obj,ipart,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      IF(idone == 1) EXIT
      CALL AIS(obj,dat(iping),NTEMP1,NBAL1)
      tendcmp = MPI_WTIME()
      tcmp = tendcmp-tstartcmp
      !! Get updated particle back to master
      CALL RCV_PARTICLE(obj,isource,ipart,ifail,tcmp)
   ENDDO

ENDIF !! MPI IF
ENDIF !! iskipais IF

!!
!! Compute linear rotation for each dimension (using best fit particle), then BCAST
!!
IF(ILINROT == 1)THEN
  IF(rank == src)CALL COMPUTE_VV2(particles,dat(iping))
  PRINT*,rank,'START SND_VV 3'
  CALL SND_VV()
  PRINT*,rank,'DONE SND_VV 3'
ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!                                                                           !!
!! AGAIN DO SOME MCMC after AIS, THEN RESAMPLE                               !!
!!                                                                           !!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
   WRITE(6,*) 'rjMCMC:'
   WRITE(6,203) '   iping,  iloop,particle,idone,        logL,  k,iacc_bd, acc_ratio, time(prtcle),source, time(comp),ifail,nfail'
   WRITE(6,203) '----------------------------------------------------------------------------------------------------------------'
   ipartsnd = 0
   idxfail = 0
   particles = particles_new  !! Save the resampled particles so we can get them back in case model fails.
   DO iloop = 1,MIN(NTHREAD-1,NPARTAIS)
      ipartsnd = ipartsnd + 1
      idest = ipartsnd
      idone = 0
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
   ENDDO ! PARTICLE LOOP
   iloop = 0
   DO 
      iloop = iloop + 1
      tstartsnd = MPI_WTIME()
      CALL RCV_PARTICLE(objnew,isource,ipartrcv,ifail,tcmp)
      tendsnd = MPI_WTIME()
      particles(ipartrcv) = objnew
      idest = isource

      IF(MOD(iloop,100)==0)THEN
         WRITE(6,201) iiping,iloop,ipartrcv,ipartsnd,particles(ipartrcv)%logL,particles(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
         WRITE(44,201) iiping,iloop,ipartrcv,ipartsnd,particles(ipartrcv)%logL,particles(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
      ENDIF
      IF(ipartsnd < NPARTAIS)THEN
         tendsnd = MPI_WTIME()
         ipartsnd = ipartsnd + 1
         idone = 0
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ELSE
         idone = 1
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ENDIF

      IF(iloop == NPARTAIS)EXIT
   ENDDO ! PARTICLE LOOP

   kmintmp = NLMX
   kmaxtmp = 0
   DO ipart = 1,NPART
      !! Compute likelihoods for new ping for all particles
!      CALL LOGLHOOD(particles(ipart),dat(iping))
      !! Compute particles log-weights (ratio of likelihoods or difference in logL)
!      particles(ipart)%logwt = particles(ipart)%logL-particles_old(ipart)%logL
      !! THIS IS WHAT YARDIM DOES I THINK; SO ONLY LIKELIHOOD...
      particles(ipart)%logwt = particles(ipart)%logL
      kmintmp = MIN(kmintmp,particles(ipart)%k)
      kmaxtmp = MAX(kmaxtmp,particles(ipart)%k)
   ENDDO ! PARTICLE LOOP
   IF(kmintmp > 1)    kmintmp = kmintmp - 1
   IF(kmaxtmp < NLMX) kmaxtmp = kmaxtmp + 1

   WRITE(6,*) 'kmin for ping',iiping,'= ',kmintmp
   WRITE(6,*) 'kmax for ping',iiping,'= ',kmaxtmp
   !!
   !! Resample new set of particles according to weight
   !! resampled ones are particles_new
   !!
   CALL RESAMPLE(particles,particles_new,kcount,kidx)
   particles = particles_new

   !!
   !! Check if at least some particles are in target region. If not, repeat with AIS three times as many dists.
   !!
   itrgtcnt = 0
   logLmin = HUGE(1._RP)
   logLmax = -HUGE(1._RP)
   DO ipart = 1,NPART
      IF(particles(ipart)%logL > logLtrgt)itrgtcnt = itrgtcnt + 1
      IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
      IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
   ENDDO
   irealloc = 0
   IF((itrgtcnt < npercnt*FLOOR(REAL(NPART,RP)/100._RP)).AND.(iais < NCOOLTRAJ))THEN
      iskipais = 0
      WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
      WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
      WRITE(6,*) itrgtcnt,'particles are in target region. Proceeding with AIS.'
   ELSE
      iskipais = 1
      WRITE(6,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
      WRITE(6,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
      WRITE(6,*) itrgtcnt,'particles are in target region. Skipping AIS.'
   ENDIF
ELSE
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!    INITIAL rjMCMC SLAVE PART
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   DO
      ifail = 0 
      idone = 0
      !! Reset acceptance counters:
      iaccept_bd     = 0
      ireject_bd     = 0
      i_bd           = 0
      i_zpert        = 0
      tstartcmp = MPI_WTIME()
      !! Get new particle from master
      CALL SND_PARTICLE(obj,ipart,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)

      IF(idone == 1) EXIT

      DO imcmc = 1,MCMC_INI
         CALL EXPLORE_MH(obj,dat(iping),1._RP)
      ENDDO
      tendcmp = MPI_WTIME()
   
      tcmp = tendcmp-tstartcmp
      !! Get updated particle back to master
      CALL RCV_PARTICLE(obj,isource,ipart,ifail,tcmp)
   ENDDO
ENDIF !! MPI IF

CALL MPI_BCAST( logLtrgt,1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( kmintmp,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( kmaxtmp,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

IF(ICOOL == 1)THEN
   !! Make cooling schedule
   NTEMP1 = NTEMP
   NBAL1  = NBAL
   ALLOCATE( beta(NTEMP1) )
   beta1 = beta1gl
   beta4 = beta4gl
   CALL MAKE_BETA(beta1,beta4,NTEMP1)
ENDIF

!!
!! Compute linear rotation for each dimension (using best fit particle), then BCAST
!!
IF(ILINROT == 1)THEN
  IF(rank == src)CALL COMPUTE_VV(particles,dat(iping))
  PRINT*,rank,'START SND_VV 4'
  CALL SND_VV()
  PRINT*,rank,'DONE SND_VV 4'
ENDIF

IF(rank==src)THEN
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!     MASTER PART main rjMCMC balancing
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   WRITE(6,*) 'Target logL = ',logLtrgt
   WRITE(6,*)   ''

   DO ipart = 1,NPART
      particles_new(ipart)%nfail = 0
   ENDDO ! PARTICLE LOOP
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!                                                                                !!
   !! Send new particles to slaves, then receive answers (auto-load balancing):      !!
   !!                                                                                !!
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   
   particles = particles_new  !! Save the resampled particles so we can get them back in case model fails.
   WRITE(6,*) 'rjMCMC:'
   WRITE(6,203) '   iping,  iloop,particle,idone,        logL,  k,iacc_bd, acc_ratio, time(prtcle),source, time(comp),ifail,nfail'
   WRITE(6,203) '----------------------------------------------------------------------------------------------------------------'
   ipartsnd = 0
   idxfail = 0
   DO iloop = 1,MIN(NTHREAD-1,NPART)
      ipartsnd = ipartsnd + 1
      idest = ipartsnd
      idone = 0
      CALL SND_PARTICLE(particles_new(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
   ENDDO ! PARTICLE LOOP
   iloop = 0
   DO 
      iloop = iloop + 1
!      IF(MOD(iloop,500+1)==0)THEN
!         WRITE(6,*) 'rjMCMC:'
!         WRITE(6,203) '   iping,  iloop,particle,idone,        logL,  k,iacc_bd, acc_ratio, time(prtcle),source, time(comp),ifail,nfail'
!         WRITE(6,203) '----------------------------------------------------------------------------------------------------------------'
!      ENDIF
      tstartsnd = MPI_WTIME()
      CALL RCV_PARTICLE(objnew,isource,ipartrcv,ifail,tcmp)
      tendsnd = MPI_WTIME()
      particles_new(ipartrcv) = objnew
      idest = isource

      IF(MOD(iloop,500)==0)THEN
         WRITE(6,201) iiping,iloop,ipartrcv,ipartsnd,particles_new(ipartrcv)%logL,particles_new(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
         WRITE(44,201) iiping,iloop,ipartrcv,ipartsnd,particles_new(ipartrcv)%logL,particles_new(ipartrcv)%k,&
                      iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
      ENDIF
      IF(ipartsnd < NPART)THEN
         tendsnd = MPI_WTIME()
         ipartsnd = ipartsnd + 1
         idone = 0
         CALL SND_PARTICLE(particles_new(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ELSE
         idone = 1
         CALL SND_PARTICLE(particles_new(ipartsnd),ipartsnd,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
      ENDIF

      IF(iloop == NPART)EXIT
   ENDDO ! PARTICLE LOOP
   particles = particles_new
   !!
   !! Leapfrog a bit
   !!
!   logLmin = HUGE(1._RP)
!   logLmax = -HUGE(1._RP)
!   DO ipart = 1,NPART
!      IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
!      IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
!   ENDDO
!   WRITE(*,'(a,F12.4,F12.4)')'MCMC: Before LF logL',logLmin,logLmax
!   CALL EXPLORE_LEAPFROG(particles,dat)
!   logLmin = HUGE(1._RP)
!   logLmax = -HUGE(1._RP)
!   DO ipart = 1,NPART
!      IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
!      IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
!   ENDDO
!   WRITE(*,'(a,F12.4,F12.4)')'MCMC: Before LF logL',logLmin,logLmax

   !!
   !!  File to save posterior particle sample
   !!
   OPEN(UNIT=40,FILE=particlefile(iping),FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
   WRITE(6,*) 'Writing to file ',particlefile(iping)
   DO ipart = 1,NPART
      WRITE(40,202) particles_new(ipart)%logL, particles_new(ipart)%logPr, particles_new(ipart)%lognorm, &
                    REAL(particles_new(ipart)%k,RP),particles_new(ipart)%par,particles_new(ipart)%sdpar
   ENDDO
   CALL FLUSH(40)
   CLOSE(40)
   tiping2 = MPI_WTIME()
   WRITE(6,*) '# particles failed at this ping:',SUM(idxfail)
   WRITE(6,*) 'Total time for ping:            ',tiping2-tiping1

ELSE
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!    MAIN rjMCMC RE-BALANCING SLAVE PART
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   DO               !! Keep working until master tells you you're done
      ifail = 0 
      idone = 0
      !! Reset acceptance counters:
      iaccept_bd     = 0
      ireject_bd     = 0
      i_bd           = 0
      i_zpert        = 0
      tstartcmp = MPI_WTIME()
      !! Get new particle from master
      CALL SND_PARTICLE(obj,ipart,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)

      IF(idone == 1) EXIT

      IF(ICOOL == 1)THEN
         objmax = obj
         DO imcmc = 1,NTEMP1
            DO ithin = 1,NBAL
               CALL EXPLORE_MH(obj,dat(iping),beta(imcmc))
               IF(obj%logL > objmax%logL) objmax=obj
            ENDDO
            IF(objmax%logL > logLtrgt)EXIT
         ENDDO
         obj = objmax
      ENDIF
      imcmc = 0
      DO 
         CALL EXPLORE_MH(obj,dat(iping),1._RP)
         imcmc = imcmc + 1
         IF(obj%logL > logLtrgt)EXIT
         IF(imcmc >= MCMC_BI)EXIT
      ENDDO
      DO imcmc = 1,MCMC_STEPS
         CALL EXPLORE_MH(obj,dat(iping),1._RP)
      ENDDO
      tendcmp = MPI_WTIME()
   
      tcmp = tendcmp-tstartcmp
      !! Get updated particle back to master
      CALL RCV_PARTICLE(obj,isource,ipart,ifail,tcmp)
   ENDDO

ENDIF !! MPI ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
iiping = iiping+NPAVE
ENDDO !! PING LOOP
CLOSE(44)  !! Closing log file

IF(rank == src) WRITE(6,*) 'SRJMH run finished'

201 FORMAT(I8,I8,I9,I6,F13.4,I4,I8,F11.4,F14.4,I7,F12.4,I6,I5)
202 FORMAT(500ES18.8)
203 FORMAT(A112)
204 FORMAT(I6,A4,F12.6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
205 FORMAT(A11,I3,A11,I3)
206 FORMAT(A52,I5,A60)
207 FORMAT(A117)
208 FORMAT(A29,A15,A5,A15,A1)
209 FORMAT(I8,I8,I9,4F13.4,I7,F12.4)
210 FORMAT(A96)
211 FORMAT(A,I4)
212 FORMAT(A,I3,A,2F8.2)
213 FORMAT(A,3F10.4)
214 FORMAT(A,F8.2)
215 FORMAT(6F12.6)
216 FORMAT(A56)

CALL MPI_FINALIZE( ierr ) 
END PROGRAM SRJMCMC_PLANE

!=======================================================================

SUBROUTINE LOGLHOOD(obj,dat)
!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB) :: ifreq,ifr,iang,ncra,ipar,iidx
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP),DIMENSION(obj%NFP)            :: m_inr
REAL(KIND=RP),DIMENSION(NAVEF,obj%k+1)      :: vp,alf1dB
REAL(KIND=RP),DIMENSION(NAVEF,NANG)         :: Rpltry
REAL(KIND=RP), DIMENSION(NFPMX)             :: mtmp
REAL(KIND=RP),DIMENSION(NBAND)              :: Etmp
INTEGER(KIND=IB),DIMENSION(obj%k)           :: idxh
REAL(KIND=RP)                               :: fref,cref,alfrdB,y

!!
!!  Compute plane wave refl. coeff. (band averaged)
!!
!IF(i_sdpert == 0)THEN

CALL CALCH(obj)
mtmp = 0._RP
mtmp(1:obj%NFP) = obj%par(1:obj%NFP)
idxh = 0
idxh = (/1:obj%k*NPL:NPL/)
mtmp(idxh) = obj%h(1:obj%k)

   DO ifreq = 1,NBAND
      IF(ISPHER == 0)THEN
         IF(NAVEF > 1)THEN
            CALL REF_NLAY3(dat%angobs,mtmp(1:obj%NFP),fr((ifreq-1)*NAVEF+1:ifreq*NAVEF),&
                           Rpltry,cw,rw,NAVEF,NANG,obj%NFP)
         ELSE
            CALL REF_NLAY3(dat%angobs,mtmp(1:obj%NFP),fr(ifreq),&
                           Rpltry,cw,rw,NAVEF,NANG,obj%NFP)
         ENDIF
      ELSE
         IF(NAVEF > 1)THEN
            CALL SPH_REF_NLAY(dat%angobs,mtmp(1:obj%NFP),fr((ifreq-1)*NAVEF+1:ifreq*NAVEF),&
                              Rpltry,cw,rw,z_t,NAVEF,NANG,obj%NFP,ifreq)
         ELSE
            CALL SPH_REF_NLAY(dat%angobs,mtmp(1:obj%NFP),fr(ifreq),&
                              Rpltry,cw,rw,z_t,NAVEF,NANG,obj%NFP,ifreq)
         ENDIF
      ENDIF

      IF(NAVEF > 1)THEN
         IF(IGA == 1)THEN
            DO iang = 1,NANG
               !!
               !! GAUSSIAN AVERAGE
               !!
               dat%Rrep(ifreq,iang) = SUM(Rpltry(:,iang)*EXP(-(fr((ifreq-1)*NAVEF+1:ifreq*NAVEF)-bands(ifreq))**2._RP/ &
                                      (frbw*bands(ifreq))**2._RP)*fstep(ifreq),1)/ &
                                      SUM(EXP(-(fr-bands(ifreq))**2._RP/ &
                                     (frbw*bands(ifreq))**2._RP)*fstep(ifreq),1);
            ENDDO 
         ELSE
            !!
            !! INTENSITY AVERAGE
            !!
            dat%Rrep(ifreq,:) = SQRT(SUM(Rpltry*Rpltry,1)/REAL(NAVEF,KIND=RP))
         ENDIF
      ELSE
         dat%Rrep(ifreq,:) = ABS(Rpltry(1,:))
      ENDIF
   ENDDO
   IF(IdB == 1)THEN
      dat%Rrep = -20._RP*LOG10(ABS(dat%Rrep))
   ENDIF

   dat%res = (dat%Robs-dat%Rrep)!*dat%Rex
   obj%res = dat%res
!ELSE
!   dat%res = obj%res
!ENDIF

!!
!!  Compute log likelihood
!!
IF(ICOV == 1)THEN
   !!
   !! Sample over sigma (one for all freqs)
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = -(REAL(NANG,RP)/2._RP)*LOG((2._RP*PI)) &
          -(SUM(dat%res(ifreq,:)**2._RP)/(2._RP*obj%sdpar(ifreq)**2._RP)+REAL(NANG,RP)*LOG(obj%sdpar(ifreq)))
   ENDDO
ELSEIF(ICOV == 2)THEN
   !!
   !! Keep sigma fixed at ping avergaed estimate
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = dat%lognorm(ifreq) &
          -(SUM(dat%res(ifreq,:)**2._RP/(2._RP*dat%sigma(ifreq,:)**2._RP)))
   ENDDO
ELSEIF(ICOV == 3)THEN
   !!
   !! Sample over scaling (one for all freqs) for ping averaged sigma estimate
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = -(REAL(NANG,RP)/2._RP)*LOG((2._RP*PI)) &
          -(SUM(dat%res(ifreq,:)**2._RP/(2._RP*(obj%sdpar(1)*dat%sigma(ifreq,:))**2._RP))&
          +REAL(NANG,RP)*LOG(obj%sdpar(ifreq))+SUM(LOG(dat%sigma(ifreq,:))))
   ENDDO
ENDIF
obj%logL = SUM(Etmp)
obj%logP = LOG(PRODUCT(1._RP/maxpert))

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD
!==============================================================================

SUBROUTINE LOGLHOOD2(obj,dat)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat

obj%logL = 1._RP

END SUBROUTINE LOGLHOOD2
!==============================================================================

SUBROUTINE EXPLORE_MH(obj,dat,beta_mh)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,NFPnew,idel,ipar,ilay,idxz,infloop
INTEGER(KIND=IB)                            :: ncra
TYPE(objstruc)                              :: obj,objnew
TYPE(datastruc)                             :: dat
REAL(KIND=RP)                               :: logPLratio,logPQ,logy
REAL(KIND=RP)                               :: ran_uni,ran_uni_sd,ran_uni_BD, ran_unik
REAL(KIND=RP)                               :: znew,beta_mh
INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp

objnew = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)

!! Do normal MCMC with 0.5 probability
IF(ran_uni_BD <= 0.5_RP)THEN
!! Do nothing to k and z
   i_zpert = 0
   i_bd = 0
ENDIF

!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(ran_uni_BD > 0.5_RP)THEN 
   i_zpert = 1
   !! Perturbing k:
   CALL RANDOM_NUMBER(ran_unik)
   i_bd = 0
   IF(obj%k == kmaxtmp)THEN  !! If k == kmax, no birth allowed
                             !! kmaxtmp is determined on the fly from first particle cloud
      IF(ran_unik>=0.5_RP) i_bd = 2
   ELSEIF(obj%k == kmintmp)THEN  !! If k == kmin, no death allowed
                                 !! kmintmp is determined on the fly from first particle cloud
      IF(ran_unik>=0.5_RP) i_bd = 1
   ELSE
      IF((ran_unik <= 1._RP/3._RP)) i_bd = 1
      IF((ran_unik > 1._RP/3._RP) .AND. (ran_unik <= 2._RP/3._RP)) i_bd = 2
   ENDIF
   IF(i_bd == 1)THEN
      CALL BIRTH(obj,objnew)
   ELSEIF(i_bd == 2)THEN
      CALL DEATH(obj,objnew)
   ENDIF
ENDIF

!!
!! Do Metropolis-Hastings
!!
IF(obj%k /= objnew%k)THEN
   IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)
   CALL CHECKBOUNDS(objnew)
   IF(ioutside == 0)THEN
      CALL LOGLHOOD(objnew,dat)
      IF(i_bd == 1)THEN
         !! BIRTH ACCEPTANCE RATIO
         logPQ = -objnew%logPr + objnew%lognorm + &
                 (0.5_RP*DOT_PRODUCT(MATMUL((objnew%gp-objnew%g),objnew%Chati),(objnew%gp-objnew%g)))
!         logPLratio =  logPQ + (objnew%logL - obj%logL)*beta_mh
         !!
         !! Raising the whole acceptance to the power of beta should balance properly over all model choices at all temperatures.
         !!
         logPLratio =  (logPQ + objnew%logL - obj%logL)*beta_mh
      ELSEIF(i_bd == 2)THEN
         !! DEATH ACCEPTANCE RATIO
         logPQ = objnew%logPr - objnew%lognorm - &
                 (0.5_RP*DOT_PRODUCT(MATMUL((objnew%gp-objnew%g),objnew%Chati),(objnew%gp-objnew%g)))
         !!
         !! Raising the whole acceptance to the power of beta should balance properly over all model choices at all temperatures.
         !!
         logPLratio = (logPQ + objnew%logL - obj%logL)*beta_mh
      ENDIF
      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
         objnew = obj
         ireject_bd = ireject_bd + 1
      ELSE
         obj = objnew
         iaccept_bd = iaccept_bd + 1
         IF(idmbest(obj%k) == 0)THEN
           IF(ILINROT == 1)THEN
             CALL LINROT(obj,dat)
!             WRITE(*,*) rank,'updating Linrot in BD, k=',obj%k
           ENDIF
         ENDIF
      ENDIF 
   ELSE
      objnew = obj
      ireject_bd = ireject_bd + 1
      ioutside = 0
   ENDIF
ELSE  ! k-change if
   idxrand = 0
   idxrand(1:obj%NFP) = RANDPERM(obj%NFP)
!   IF(ISLICE == 0)THEN
      !!
      !! Do Metropolis-Hastings update on c, rho, alpha
      !!
      DO ipar = 1,obj%NFP
         iwhich = idxrand(ipar)
         CALL PROPOSAL(obj,objnew,iwhich)
         IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)
         CALL CHECKBOUNDS(objnew)
         IF(ioutside == 0)THEN
            CALL LOGLHOOD(objnew,dat)
            logPLratio = (objnew%logL - obj%logL)*beta_mh
            CALL RANDOM_NUMBER(ran_uni)
            IF(ran_uni >= EXP(logPLratio))THEN
               objnew = obj
               ireject = ireject + 1
            ELSE
               obj = objnew
               iaccept = iaccept + 1
            ENDIF
         ELSE
            objnew = obj
            ireject = ireject + 1
            ioutside = 0
         ENDIF 
      ENDDO
!   ELSE ! SLICE if
!      !!
!      !! Do slice sampling update on c, rho, alpha
!      !!
!      DO ipar = 1,ncra
!         iwhich = idxcra(idxrand(ipar))
!         CALL RANDOM_NUMBER(ran_uni)
!         logy = obj%logL + LOG(ran_uni)
!         CALL EXPLORE_SLICE(obj,dat,logy,iwhich)
!      ENDDO
!      obj = objnew
!   ENDIF
   !!
   !! Do Metropolis-Hastings on data-error standard deviations
   !!
   IF(ICOV == 1)THEN
      !! Perturb std devs with .25 probability
      CALL RANDOM_NUMBER(ran_uni_sd)
      IF(ran_uni_sd>=0.25_RP)THEN
         DO ipar = 1,NBAND
            CALL PROPOSAL_SD(obj,objnew,ipar)
!            i_sdpert = 1
            IF(ioutside == 0)THEN
               CALL LOGLHOOD(objnew,dat)
               logPLratio = (objnew%logL - obj%logL)*beta_mh
               CALL RANDOM_NUMBER(ran_uni)
               IF(ran_uni >= EXP(logPLratio))THEN
                  objnew = obj
                  ireject = ireject + 1
               ELSE
                  obj = objnew
                  iaccept = iaccept + 1
               ENDIF
            ELSE
               objnew = obj
               ireject = ireject + 1
               ioutside = 0
            ENDIF
            i_sdpert = 0
         ENDDO
      ENDIF
   ENDIF ! ICOV if
ENDIF ! k-change if

END SUBROUTINE EXPLORE_MH
!!==============================================================================

SUBROUTINE EXPLORE_LEAPFROG(particles,dat)
!!==============================================================================
!!
!! Leapfrog MCMC according to Skilling/MacKay 2003
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: iipart,ipart,ik,ilf
TYPE(objstruc),DIMENSION(NPART)             :: particles     ! All particles (current ping)
TYPE(objstruc),DIMENSION(NPART)             :: particlesnew  ! All particles (current ping)
TYPE(objstruc)                              :: objs,objt,objsp
TYPE(datastruc)                             :: dat
REAL(KIND=RP)                               :: logPLratio,ran_uni
INTEGER(KIND=IB),DIMENSION(NLMX)            :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX)      :: kidx
INTEGER(KIND=IB),DIMENSION(NPART)           :: idxrand

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Find no. particles per k
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
kcount = 0
kidx   = 0
DO ik = 1,NLMX
   iipart = 0
   DO ipart = 1,NPART
      IF(particles(ipart)%k == ik)THEN
         iipart = iipart + 1
         kcount(ik) = kcount(ik)+1
         kidx(iipart,ik) = ipart
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Leapfrog models in each dimension
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
DO ik = 1,NLMX
   IF(kcount(ik) >= 100)THEN

      DO ilf = 1,FLOOR(REAL(kcount(ik),RP)/10._RP)
         idxrand = 0
         idxrand = RANDPERM(kcount(ik))
         objs = particles(kidx(idxrand(1),ik))
         objt = particles(kidx(idxrand(2),ik))
         !!
         !! Propose leapfrog step:
         !!
         objsp = objs
         objsp%par = 2._RP*objt%par-objs%par
         CALL CALCH(objsp)
         IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objsp)
         CALL CHECKBOUNDS(objsp)
         IF(MAXVAL(objsp%par-objs%par) == 0._RP) ioutside = 1
         !!
         !! Evaluate with Metropolis-Hastings:
         !!
         IF(ioutside == 0)THEN
            CALL LOGLHOOD(objsp,dat)
            logPLratio = (objsp%logL - objs%logL)
            CALL RANDOM_NUMBER(ran_uni)
            IF(ran_uni >= EXP(logPLratio))THEN
               objsp = objs
               irejectlf = irejectlf + 1
            ELSE
               objs = objsp
               iacceptlf = iacceptlf + 1
            ENDIF
         ELSE
            objsp = objs
            irejectlf = irejectlf + 1
            ioutside = 0
         ENDIF 
         particles(kidx(idxrand(1),ik)) = objs
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE EXPLORE_LEAPFROG
!=======================================================================

SUBROUTINE EXPLORE_SLICE(obj,dat,logy,iwhich)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=IB)  :: j,k,jcount,iwhich,ilay,ipar
TYPE (objstruc) :: obj,try
TYPE (datastruc):: dat
REAL(KIND=RP)  :: logy
REAL(KIND=RP)  :: ran_uni,u
REAL(KIND=RP)  :: L,R,x0,x1,gx0

ilay = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
ipar = iwhich-(ilay-1)*NPL
IF(ilay > obj%k) ipar = ipar + 1

jcount = 0
try = obj
gx0 = obj%logL
x0  = obj%par(iwhich)

!!
!!  Slice level = logy
!!  Find initial sampling interval
!!
CALL RANDOM_NUMBER(ran_uni)
u = ran_uni * w_sl(ipar)
L = x0 - u
R = x0 + (w_sl(ipar)-u)

!!
!!  Expand interval until outside slice
!!
DO
   IF(L <= minlim(ipar))EXIT
   try%par(iwhich) = L
   CALL LOGLHOOD(try,dat)
   icount(rank+1) = icount(rank+1) + 1
   jcount = jcount + 1
   IF(try%logL <= logy)EXIT
   L = L - w_sl(ipar)
   IF(jcount > 5000)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop1'
ENDDO
jcount = 0
DO
   IF(R >= maxlim(ipar))EXIT
   try%par(iwhich) = R
   CALL LOGLHOOD(try,dat)
   icount(rank+1) = icount(rank+1) + 1
   jcount = jcount + 1
   IF(try%logL <= logy)EXIT
   R = R + w_sl(ipar)
   IF(jcount > 5000)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop2'
ENDDO

!!
!!  Shrink interval to lower and upper bounds
!!
IF(L < minlim(ipar)) L = minlim(ipar)
IF(R > maxlim(ipar)) R = maxlim(ipar)

!!
!!  Sample from interval, shrinking it on each rejection
!!
jcount = 0

DO
   CALL RANDOM_NUMBER(ran_uni)
   x1 = L+(R-L)*ran_uni
   try%par(iwhich) = x1
   CALL LOGLHOOD(try,dat)
   icount(rank+1) = icount(rank+1) + 1
   jcount = jcount + 1
   IF(try%logL >= logy)EXIT
   IF(x1 > x0)THEN
      R = x1
   ELSE
      L = x1
   ENDIF
   IF(jcount > 100)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop3'
ENDDO
obj = try

RETURN
END SUBROUTINE EXPLORE_SLICE
!!==============================================================================

SUBROUTINE COMPUTE_VV(particles,dat)
!!==============================================================================
!!
!! Compute the 
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: iipart,ipart,ik,ilf,ipar,imaxpart
TYPE(objstruc),DIMENSION(NPART)             :: particles     ! All particles (current ping)
TYPE(objstruc)                              :: obj
TYPE(datastruc)                             :: dat
INTEGER(KIND=IB),DIMENSION(NLMX)            :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX)      :: kidx

!!
!! Find no. particles per k
!!
kcount = 0
kidx   = 0
DO ik = 1,NLMX
   iipart = 0
   DO ipart = 1,NPART
      IF(particles(ipart)%k == ik)THEN
         iipart = iipart + 1
         kcount(ik) = kcount(ik)+1
         kidx(iipart,ik) = ipart
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO
!!
!! Compute linear rotation for each k with sufficient models
!!
VV    = 0._RP
DO ik = kmin,kmax
DO ipar = 1,NFPMX
   VV(ipar,ipar,ik) = 1._RP
ENDDO
ENDDO

DO ik = kmin,kmax
  idmbest = 0
  dmbest  = 0._RP
  IF(kcount(ik)>0_IB)THEN
    logLmax = -HUGE(1._RP)
    DO ipart = 1,kcount(ik)
      IF(particles(kidx(ipart,ik))%logL > logLmax)THEN
        imaxpart = kidx(ipart,ik)
        logLmax = particles(kidx(ipart,ik))%logL
      ENDIF
    ENDDO
!
!    obj = particles(imaxpart)
!    CALL LOGLHOOD(obj,dat)
!    PRINT*,'Check logLmax:'
!    PRINT*,obj%logL,particles(imaxpart)%logL
!
    CALL LINROT(particles(imaxpart),dat)
  ENDIF
ENDDO

END SUBROUTINE COMPUTE_VV
!!==============================================================================

SUBROUTINE COMPUTE_VV2(particles,dat)
!!==============================================================================
!!
!! Compute the 
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: iipart,ipart,ik,ilf,ipar,imaxpart
TYPE(objstruc),DIMENSION(NPART)             :: particles     ! All particles (current ping)
TYPE(objstruc)                              :: obj
TYPE(datastruc)                             :: dat
INTEGER(KIND=IB),DIMENSION(NLMX)            :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX)      :: kidx

!!
!! Find no. particles per k
!!
kcount = 0
kidx   = 0
DO ik = 1,NLMX
   iipart = 0
   DO ipart = 1,NPART
      IF(particles(ipart)%k == ik)THEN
         iipart = iipart + 1
         kcount(ik) = kcount(ik)+1
         kidx(iipart,ik) = ipart
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO
!!
!! Compute linear rotation for each k with sufficient models
!!
VV    = 0._RP
DO ik = kmin,kmax
DO ipar = 1,NFPMX
   VV(ipar,ipar,ik) = 1._RP
ENDDO
ENDDO

DO ik = kmin,kmax
  idmbest = 0
  dmbest  = 0._RP
  IF(kcount(ik)>0_IB)THEN
    logLmax = -HUGE(1._RP)
    DO ipart = 1,kcount(ik)
      IF(particles(kidx(ipart,ik))%logL > logLmax)THEN
        imaxpart = kidx(ipart,ik)
        logLmax = particles(kidx(ipart,ik))%logL
      ENDIF
    ENDDO
    CALL LINROT2(particles(imaxpart),dat)
  ENDIF
ENDDO

END SUBROUTINE COMPUTE_VV2
!=======================================================================

SUBROUTINE CHECKBOUNDS(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: ip,iwhich,ipar,ilay,ncra
INTEGER(KIND=IB)                            :: ih
TYPE(objstruc)                              :: obj
REAL(KIND=RP)                               :: vspmin,vspmax,zhere
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp

IF(obj%par((obj%k-1)*4+1) > maxlim(1))ioutside = 1
IF(obj%par(1) < minlim(1))ioutside = 1

DO ilay = 1,obj%k
   IF(hmin > obj%h(ilay))ioutside = 1
   IF(maxlim(1) < obj%h(ilay))ioutside = 1
ENDDO

CALL GETIDXCRA(obj,idxcra,ncra)
DO ip = 1,ncra
   iwhich = idxcra(ip)
   ilay = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
   ipar = iwhich-(ilay-1)*NPL
   IF(ilay > obj%k) ipar = ipar + 1
   IF(ISINPRIOR == 0)THEN
      IF(((obj%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - obj%par(iwhich)) < 0._RP))THEN
         ioutside = 1
      ENDIF
   ELSE
      IF(ipar == 2)THEN
         zhere = obj%z(ilay)
         IF(ilay > obj%k) zhere = obj%z(obj%k)
         vspmin = vs01 + (vsh1-vs01)*SIN(zhere/maxlim(1)*PI/2._RP)**numin
         vspmax = vs02 + (vsh2-vs02)*SIN(zhere/maxlim(1)*PI/2._RP)**numax
         IF(((obj%par(iwhich) - vspmin) < 0._RP).OR.((vspmax - obj%par(iwhich)) < 0._RP))THEN
            ioutside = 1
         ENDIF
      ELSE
         IF(((obj%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - obj%par(iwhich)) < 0._RP))THEN
            ioutside = 1
         ENDIF
      ENDIF
   ENDIF
ENDDO
RETURN
END SUBROUTINE CHECKBOUNDS
!!=======================================================================
!
!SUBROUTINE PROPOSALZ(obj,objnew,idxz,iwhich)
!!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB)              :: i,iwhich,idxz,iloop
!TYPE(objstruc)                :: obj,objnew
!REAL(KIND=RP)                 :: ran_uni, ran_nor
!REAL(KIND=RP),DIMENSION(obj%k):: ztmp
!
!!! Perturb only one z with Gaussian proposal
!!CALL GASDEVJ(ran_nor)
!!objnew%z(idxz) = obj%z(idxz) + pertsd(1)*ran_nor
!!! Perturb only one z with Cauchy proposal
!CALL RANDOM_NUMBER(ran_uni)
!objnew%z(idxz) = obj%z(idxz) + pertsd(1)*TAN(PI*(ran_uni-0.5_RP))
!
!CALL CALCPAR(objnew)
!
!!! Check that z is not perturbed into another layer; reject if it happened
!IF(idxz ==1)THEN
!   IF(objnew%z(idxz) > obj%z(idxz+1))ioutside = 1
!ELSEIF(idxz == objnew%k)THEN
!   IF(objnew%z(idxz) < objnew%z(idxz-1))ioutside = 1
!ELSE
!   IF((objnew%z(idxz) < objnew%z(idxz-1)).OR.(objnew%z(idxz) > obj%z(idxz+1)))ioutside = 1
!ENDIF
!
!RETURN
!END SUBROUTINE PROPOSALZ
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,iwhich,ipar,ipar2,iloop,ilay,ilay2
TYPE(objstruc) :: obj,objnew,objtmp
REAL(KIND=RP)  :: ran_uni, ran_nor
REAL(KIND=RP),DIMENSION(obj%NFP) :: m2p,m1
REAL(KIND=RP),DIMENSION(obj%NFP,obj%NFP) :: VVtmp

VVtmp = 0._RP
VVtmp = VV(1:obj%NFP,1:obj%NFP,obj%k)
objnew = obj
ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
ipar  = iwhich-(ilay-1)*NPL
IF(ilay > obj%k) ipar = ipar + 1
!!
!! Gaussian proposal
!!
!CALL GASDEVJ(ran_nor)
!objnew%par(iwhich) = obj%par(iwhich) + pertsd(ipar)*ran_nor
!!
!! CAUCHY proposal
!!
CALL RANDOM_NUMBER(ran_uni)

IF(ILINROT == 0)THEN
   objnew%par(iwhich) = obj%par(iwhich) + pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
ELSE
   IF(idmbest(objnew%k) == 0)THEN
     objnew%par(iwhich) = obj%par(iwhich) + pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
   ELSE
     !!
     !! Propose in rotated space:
     !!
     DO i=1,obj%NFP
       ilay2  = CEILING(REAL(i,RP)/REAL(NPL,RP))
       ipar2  = i-(ilay2-1)*NPL
       IF(ilay2 > obj%k) ipar2 = ipar2 + 1
       m1(i) = (obj%par(i) - minlim(ipar2))/(maxlim(ipar2)-minlim(ipar2))
     ENDDO
     m2p = MATMUL(TRANSPOSE(VVtmp),m1)
     m2p(iwhich) = m2p(iwhich)+fact*sdevm(iwhich,obj%k)*TAN(PI*(ran_uni-0.5_RP))
     m1 = MATMUL(VVtmp,m2p)
     DO i=1,obj%NFP
       ilay2  = CEILING(REAL(i,RP)/REAL(NPL,RP))
       ipar2  = i-(ilay2-1)*NPL
       IF(ilay2 > obj%k) ipar2 = ipar2 + 1
       objnew%par(i) = m1(i)*(maxlim(ipar2)-minlim(ipar2))+minlim(ipar2)
     ENDDO
   ENDIF
ENDIF
CALL CALCH(objnew)

RETURN
END SUBROUTINE PROPOSAL
!=======================================================================

SUBROUTINE PROPOSAL_SD(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdpar(iwhich) = obj%sdpar(iwhich) + pertsdsd(iwhich)*ran_nor
IF(((objnew%sdpar(iwhich) - minlimsd(iwhich)) < 0._RP).OR.((maxlimsd(iwhich) - objnew%sdpar(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SD
!=======================================================================

SUBROUTINE GETIDXCRA(obj,idxcra,ncra)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,ilay,ipl,isd,ncra
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)) :: idxcra

!! Layers:
idxcra = 0
ipar = 1
DO ilay = 1,obj%k
   DO ipl = 2,NPL
      idxcra(ipar) = (ilay-1)*NPL+ipl
      ipar = ipar + 1
   ENDDO
ENDDO
!! Half-space:
DO ipl = 1,NPL-1
   idxcra(ipar) = obj%k*NPL+ipl
   ipar = ipar + 1
ENDDO
ncra = (obj%k+1)*(NPL-1)

END SUBROUTINE GETIDXCRA
!!=======================================================================
!
!SUBROUTINE CALCPAR(obj)
!!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: i
!TYPE(objstruc) :: obj
!INTEGER(KIND=IB),DIMENSION(obj%k) :: idxh
!REAL(KIND=IB),DIMENSION(obj%k) :: h
!
!idxh = 0
!idxh = (/1:obj%k*4:4/)
!h(1) = obj%z(1)
!DO i = 2,obj%k
!   h(i) = obj%z(i)-obj%z(i-1)
!ENDDO
!obj%par(idxh) = h
!
!END SUBROUTINE CALCPAR
!=======================================================================

SUBROUTINE CALCH(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxz
REAL(KIND=IB),DIMENSION(obj%k) :: h

idxz = 0
idxz = (/1:obj%k*NPL:NPL/)

obj%z = 0._RP
obj%z(1:obj%k) = obj%par(idxz)

obj%h = 0._RP
obj%h(1) = obj%par(1)
DO i = 2,obj%k
   obj%h(i) = obj%par(idxz(i))-obj%par(idxz(i-1))
ENDDO

END SUBROUTINE CALCH
!=======================================================================

SUBROUTINE CALCM(obj,mtmp)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxh
REAL(KIND=IB),DIMENSION(NFPMX) :: mtmp

mtmp = 0._RP
mtmp = obj%par
idxh = 0
idxh = (/1:obj%k*NPL:NPL/)
mtmp(idxh) = obj%h(1:obj%k)

END SUBROUTINE CALCM
!=======================================================================

SUBROUTINE DEATH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=RP),DIMENSION(NPL-1)   :: cra_ave

objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)
!! Pick random layer:
idxdeath = 0
idxdeath(1:obj%k) = RANDPERM(obj%k)
idel = idxdeath(1)

!! Delete interface (setting new layer properties to one of old layers properties at random ):
CALL RANDOM_NUMBER(ran_uni)
IF(ran_uni>=0.5_RP)idel = idel+1

objnew%par = 0._RP  ! Take care here that there are no old "left overs" from other dimentsions
cra_ave = 0._RP
IF(idel == 1)THEN
   !! Killed interface is first interface:
   cra_ave = (obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL) + obj%par(idel*NPL+2:idel*NPL+NPL))/2._RP
   objnew%par(1:objnew%NFP) = obj%par(NPL+1:obj%NFP)
   objnew%par(2:NPL) = cra_ave
!   objnew%par(1) = obj%par(1) + obj%par(NPL+1)
   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par(2:NPL)
   objnew%gp = objnew%par(2:NPL)
ELSEIF(idel >= obj%k)THEN
   !! Killed interface is last interface:
   cra_ave = (obj%par(obj%NFP-((NPL-1)*2-1):obj%NFP-(NPL-1)) + obj%par(obj%NFP-(NPL-2):obj%NFP))/2._RP
   objnew%par(1:objnew%NFP-(NPL-1)) = (/ obj%par(1:obj%NFP-(NPL*2-1)) /)
   objnew%par(objnew%NFP-(NPL-2):objnew%NFP) = cra_ave
   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par(obj%NFP-(NPL-2):obj%NFP)
   objnew%gp = objnew%par(objnew%NFP-(NPL-2):objnew%NFP)
ELSE
   !! Killed interface is in normal layer stack:
   cra_ave = (obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL) + obj%par(idel*NPL+2:idel*NPL+NPL))/2._RP
   objnew%par(1:objnew%NFP) = (/ obj%par(1:(idel-1)*NPL), obj%par(idel*NPL+1:obj%NFP) /)
   objnew%par((idel-1)*NPL+2:(idel-1)*NPL+NPL) = cra_ave
!   objnew%par((idel-1)*NPL+1) = obj%par((idel-1)*NPL+1) + obj%par((idel-1)*NPL+(NPL+1))
   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
   objnew%gp = objnew%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
ENDIF
CALL CALCH(objnew)

END SUBROUTINE DEATH
!=======================================================================

SUBROUTINE BIRTH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)              :: i, iznew, ipert, iloop, iwhich, ipar
TYPE(objstruc)                :: obj,objnew,objnew2
REAL(KIND=RP)                 :: znew,hnew1,hnew2,ran_uni,ran_nor
REAL(KIND=RP),DIMENSION(obj%k):: ztmp

objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)
!!
!! Draw new z until new layer > hmin
!!
CALL RANDOM_NUMBER(ran_uni)
znew = maxpert(1)*ran_uni
ztmp = obj%z(1:obj%k) - znew
IF(hmin > MINVAL(ABS(ztmp))) ioutside = 1
IF(hmin > znew) ioutside = 1
!!
!! Insert new interface
!!
iznew = 0
DO i = 1,obj%k
   IF(ztmp(i) > 0._RP) EXIT
   iznew = i
ENDDO
iznew = iznew + 1
objnew%z   = 0._RP
objnew%par = 0._RP
IF(iznew == 1)THEN
   !! New layer is first layer:
!   hnew1 = znew
!   hnew2 = obj%z(iznew)-znew
   objnew%par(1:objnew%NFP) = (/ obj%par(1:NPL),obj%par /)
   objnew%par(1)            = znew
!   objnew%par(NPL+1)        = hnew2
ELSEIF(iznew > obj%k)THEN
   !! New layer is created in half-space:
!   hnew1 = znew - obj%z(iznew-1)
   objnew%par(1:objnew%NFP) = (/ obj%par(1:obj%NFP-(NPL-1)),hnew1, &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP), &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP) /)
   objnew%par(objnew%NFP-((NPL*2)-2)) = znew
ELSE
   !! New layer is within layer stack:
!   hnew1 = znew - obj%z(iznew-1)
!   hnew2 = obj%z(iznew)-znew
   objnew%par(1:objnew%NFP) = (/ obj%par(1:(iznew)*NPL), &
                                 obj%par(((iznew-1)*NPL)+1:iznew*NPL), &
                                 obj%par(iznew*NPL+1:obj%NFP) /)
   objnew%par(((iznew-1)*NPL)+1) = znew
!   objnew%par((iznew*NPL)+1)     = hnew2
ENDIF
CALL CALCH(objnew)
!!
!! Pick one of the two new layers at random and perturb:
!!
CALL RANDOM_NUMBER(ran_uni)
IF(ran_uni <= 0.5_RP)THEN
   ipert = iznew
ELSE
   ipert = iznew+1
ENDIF

!! Perturb only one parameter (determined in EXPLORE_MH call)
objnew2 = objnew
DO ipar = 2,NPL
   iloop = 0
   iwhich = (ipert-1)*NPL+ipar
   IF(ipert > objnew%k) iwhich = iwhich - 1
   !!
   !! Gaussian proposal
   !!
   CALL GASDEVJ(ran_nor)
   objnew2%par(iwhich) = objnew%par(iwhich) + pertsd(ipar)*ran_nor

   !! This records the perturbation for the bd acceptance ratio
   objnew2%g(ipar-1)  = objnew%par(iwhich)
   objnew2%gp(ipar-1) = objnew2%par(iwhich)

   IF(((objnew2%par(iwhich) - minlim(ipar)) < 0._RP).OR. &
      ((maxlim(ipar) - objnew%par(iwhich)) < 0._RP)) ioutside = 1
ENDDO
objnew = objnew2
END SUBROUTINE BIRTH
!=======================================================================

SUBROUTINE READ_PPD_P01(obj,samplefiletmp)
!=======================================================================
!! Only master should read the sample. Particles are sent to slaves.
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: ipart
REAL(KIND=RP),DIMENSION(NAP+NFPMX+NBAND):: sample_p01
TYPE(objstruc),DIMENSION(NPART):: obj
CHARACTER(len=64):: samplefiletmp

OPEN(UNIT=33,FILE=samplefiletmp,FORM='formatted',ACTION='READ')

IF(rank == src)WRITE(6,*)'Reading particles from file ',samplefiletmp

DO ipart = 1,NPART
   READ(33,*) sample_p01
   obj(ipart)%k   = INT(sample_p01(4))
   obj(ipart)%logL= sample_p01(1)
   obj(ipart)%NFP = (obj(ipart)%k * NPL) + (NPL-1)
   obj(ipart)%par = 0._RP
   obj(ipart)%par(1:obj(ipart)%NFP) = sample_p01(5:4+obj(ipart)%NFP)
   obj(ipart)%sdpar = 0._RP
   IF(ICOV == 1)THEN
      obj(ipart)%sdpar(1:NBAND) = sample_p01(4+NFPMX+1:4+NFPMX+NBAND)
   ENDIF
   CALL CALCH(obj(ipart))   !! Compute layer thickness for z-partition for object obj
   CALL INIT_BD(obj(ipart)) !! Init birth-death covariance matrix for each particle
ENDDO
CLOSE(33)
RETURN
END SUBROUTINE READ_PPD_P01
!=======================================================================

SUBROUTINE READ_DATA(dat,filebasenow,particlefiletmp)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifreq,iang
TYPE (datastruc) :: dat
CHARACTER(len=64):: filebasenow,infile,sdfile,particlefiletmp

!OPEN(UNIT=20,FILE=filebasenow,FORM='formatted',STATUS='OLD',ACTION='READ')
!READ(20,*) filebaselen
!READ(20,*) filebase
!CLOSE(20)

infile          = filebasenow(1:filebaselen) // '.txt'
sdfile          = filebasenow(1:filebaselen) // '_sigma.txt'
particlefiletmp = filebasenow(1:filebaselen) // '_particles.txt'
logfile         = filebasegl(1:filebaselengl) // '_SRJMH.log'
seedfile        = filebasegl(1:filebaselengl) // '_seeds.log'
mapfile         = filebasegl(1:filebaselengl) // '_map.dat'
repfile         = filebasegl(1:filebaselengl) // '_rep.dat'

!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
!IF(rank == src)WRITE(6,209) 'Loading data from file...',infile
OPEN(UNIT=20,FILE=infile,FORM='formatted',STATUS='OLD',ACTION='READ')
!IF(ISPHER == 1)THEN
READ(20,*) z_t
!ENDIF
READ(20,*) cw
READ(20,*) rw
READ(20,*) hmx
DO ifreq = 1,NBAND
    READ(20,*) (dat%Robs(ifreq,iang),iang=1,NANG)
ENDDO
READ(20,*) (dat%angobs(iang),iang=1,NANG)
!DO ifreq = 1,NBAND
!    READ(20,*) (dat%Rex(ifreq,iang),iang=1,NANG)
!ENDDO
CLOSE(20)
IF(ICOV >= 2)THEN
   OPEN(UNIT=20,FILE=sdfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   DO ifreq = 1,NBAND
      READ(20,*) (dat%sigma(ifreq,iang),iang=1,NANG)
      dat%logdet(ifreq) = SUM(2._RP*LOG(dat%sigma(ifreq,:)))
      dat%lognorm(ifreq) = -0.5_RP*REAL(NANG,KIND=RP)*LOG(2._RP*PI) - 0.5_RP*dat%logdet(ifreq)
   ENDDO
   CLOSE(20)
ENDIF

!!
!! Length of particle sending array
!!
! 6 -> logL, logwt etc
! NFP, NLMX -> par, z
! 7 -> g, gp, detChat
! 18 -> Chat, Chati
! 6 -> acceptance pars, rank...
! 1 -> nfail counter...
len_snd  = 6+NFPMX+NBAND+NLMX+7+18+6+1
len_snd2 = 6+NFPMX+NBAND+NLMX+7+18+1

209 FORMAT(A26,A40)
RETURN
END SUBROUTINE READ_DATA
!=======================================================================

SUBROUTINE INIT_BD(obj)
!=======================================================================
!! INIT Birth-Death covariance matrix
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: ipar
TYPE(objstruc) :: obj
!!
!! Some factors for the acceptance ratio in birth/death case:
!!
obj%Chat    = 0._RP
obj%Chati   = 0._RP
obj%detChat = 1._RP
DO ipar = 1,NPL-1
   obj%Chat(ipar,ipar)  = pertsd(ipar+1)**2._RP
   obj%Chati(ipar,ipar) = 1._RP/pertsd(ipar+1)**2._RP
   obj%detChat = obj%detChat * obj%Chat(ipar,ipar)
ENDDO
obj%logPr = SUM(LOG(maxlim(2:NPL)-minlim(2:NPL))) ! prior ratio for change in dimension of 1
obj%lognorm = LOG(SQRT((2._RP*PI)**REAL(NPL-1,RP)*obj%detChat))

RETURN
END SUBROUTINE INIT_BD
!=======================================================================

SUBROUTINE INIT_FREQ(dat)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: ifreq
TYPE(datastruc) :: dat
REAL(KIND=RP)   :: flo,fhi

ALLOCATE( fr(NBAND*NAVEF),fstep(NBAND) )
fr = 0._RP
IF(NAVEF == 1)THEN
   fr = bands
ELSE
   DO ifreq = 1,NBAND
      IF(IGA == 1)THEN
         flo = bands(ifreq) - bands(ifreq)*frbw
         fhi = bands(ifreq) + bands(ifreq)*frbw
      ELSEIF(IGA == 2)THEN
         flo = bands(ifreq) / (2._RP**(1._RP/6._RP)) ! Use 1/3 octave bandwidth
         fhi = bands(ifreq) * (2._RP**(1._RP/6._RP)) !
      ELSE
         flo = bands(ifreq) - FBW
         fhi = bands(ifreq) + FBW
      ENDIF
      fstep(ifreq) = (fhi-flo)/(NAVEF-1)
      fr((ifreq-1)*NAVEF+1:ifreq*NAVEF) = flo + ((/ 1:NAVEF /)-1) * fstep(ifreq)
   ENDDO
ENDIF
IF(ISPHER == 1)THEN
   CALL SOMMER_INTEGRAND(dat%angobs,fr,cw,rw,z_t,NBAND*NAVEF,NANG)
ENDIF

RETURN
END SUBROUTINE INIT_FREQ
!=======================================================================

SUBROUTINE COMPUTE_MAP(obj,dat)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j
TYPE(datastruc):: dat
TYPE(objstruc) :: obj
REAL(KIND=RP):: tstart2,tend2,ran_nor

tstart2 = MPI_WTIME()
CALL LOGLHOOD(obj,dat)
tend2 = MPI_WTIME()
IF(ISIM == 1)THEN
   !! Make simulated data:
   DO i = 1,NBAND
      DO j = 1,NANG
         CALL GASDEVJ(ran_nor)
         dat%Rrep(i,j) = dat%Rrep(i,j)+ran_nor*obj%sdpar(i)
      ENDDO
   ENDDO
ENDIF
IF(rank == src)CALL SAVEREPLICA(obj,dat,repfile)
IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
END SUBROUTINE COMPUTE_MAP
!=======================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i
TYPE(objstruc) :: obj

DO i=1,obj%k
   WRITE(6,201) obj%par((i-1)*NPL+1:i*NPL)
ENDDO
WRITE(6,202) '            ',obj%par(obj%k*NPL+1:obj%k*NPL+(NPL-1))
IF(ICOV == 1)THEN
   WRITE(6,*) 'SD parameters:'
   WRITE(6,203) obj%sdpar
ENDIF

201 FORMAT(4F12.4)
202 FORMAT(A12,4F12.4)
203 FORMAT(20F12.4)
END SUBROUTINE PRINTPAR
!=======================================================================

SUBROUTINE PRIOR(obj,dat)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=IB) :: i,ipar,ilay
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP)    :: ran_uni

DO i=1,obj%NFP
   ilay = CEILING(REAL(i,RP)/REAL(NPL,RP))
   ipar = i-(ilay-1)*NPL
   IF(ilay > obj%k) ipar = ipar + 1
   CALL RANDOM_NUMBER(ran_uni)
   obj%par(i) = minlim(ipar) + maxpert(ipar) * ran_uni
ENDDO
IF(ICOV == 1)THEN
   DO i=1,NBAND
      CALL RANDOM_NUMBER(ran_uni)
      obj%sdpar(i) = minlimsd(i) + maxpertsd(i) * ran_uni
   ENDDO
ENDIF
!IF(IAR == 1)THEN
!   DO i=1,NARFP*NBAND
!      CALL RANDOM_NUMBER(ran_uni)
!      obj%arpar(i) = minlimar(i) + maxpertar(i) * ran_uni
!   ENDDO
!ENDIF
CALL CALCH(obj)
CALL LOGLHOOD(obj,dat)

RETURN
END SUBROUTINE PRIOR
!=======================================================================

SUBROUTINE LINROT(objst,dat)
!!
!!  Estimates linearized model covariance matrix 
!!
!=======================================================================
USE MPI
USE NR
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ifq,ip,ibest,ntot
INTEGER(KIND=IB):: ilay, ipar, NFP

TYPE (objstruc),DIMENSION(2):: obj          ! Objects 
TYPE (objstruc)             :: objst        ! initial object 
TYPE (datastruc)            :: dat          ! Data

REAL(KIND=RP)                                :: best,test

REAL(KIND=RP), DIMENSION(NDM)                :: dm
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: dm_start
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: p1,p2
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: dat1,dat2
REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: dRdm

REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Jac,JactJac,Mpriorinv,Temp
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: bufJac,WW,Jacscale
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: VVtmp,Winv
REAL(KIND=RP)                                :: Wmax,Tstar,sdevtmp
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Ctmp

!IF(rank == src)WRITE(*,*)'Computing linear rotation estimate...'
Tstar = 1._RP
ntot = NBAND*NANG
obj(1) = objst
obj(2) = objst

NFP = objst%NFP
ALLOCATE( dat1(ntot),dat2(ntot),dRdm(NFP,ntot,NDM),       &
          dm_start(NFP),p1(NBAND,NANG),p2(NBAND,NANG),    &
          Jac(ntot,NFP),JactJac(NFP,NFP),bufJac(ntot),    &
          Jacscale(NFP),Mpriorinv(NFP,NFP),Temp(ntot,NFP),&
          Ctmp(NFP,NFP),VVtmp(NFP,NFP),Winv(NFP,NFP),WW(NFP) )
!dm_start  = 1._RP
dm_start = objst%par(1:NFP)/10._RP

DO ip=1,NFP
   dm = 0._RP
   dm(1) = dm_start(ip) ! Estimate deriv for range of dm values
   DO i=1,NDM
      obj(1)%par = objst%par
      obj(1)%par(ip) = objst%par(ip)+dm(i)
      CALL CALCH(obj(1))
      CALL LOGLHOOD(obj(1),dat)
      p1 = dat%Rrep
      DO ifq=1,NBAND
        IF(ICOV /= 4)THEN
          dat1((ifq-1)*NANG+1:(ifq)*NANG) = p1(ifq,:)/objst%sdpar(ifq)
        ELSE
          dat1((ifq-1)*NANG+1:(ifq)*NANG) = MATMUL(dat%CdCi(:,:,ifq),p1(ifq,:))
        ENDIF
      ENDDO

      obj(2)%par = objst%par
      obj(2)%par(ip) = objst%par(ip)-dm(i)
      CALL CALCH(obj(2))
      CALL LOGLHOOD(obj(2),dat)
      p2 = dat%Rrep
      DO ifq=1,NBAND
        IF(ICOV /= 4)THEN
          dat2((ifq-1)*NANG+1:(ifq)*NANG) = p2(ifq,:)/objst%sdpar(ifq)
        ELSE
          dat2((ifq-1)*NANG+1:(ifq)*NANG) = MATMUL(dat%CdCi(:,:,ifq),p2(ifq,:))
        ENDIF
      ENDDO
      DO j=1,ntot
         IF(ABS((dat1(j)-dat2(j))/(dat1(j)+dat2(j))) > 1.0e-7_RP)THEN
            dRdm(ip,j,i) = (dat1(j)-dat2(j))/(2._RP*dm(i))
         ELSE
            dRdm(ip,j,i) = 0._RP
         ENDIF
      ENDDO
      IF(i < NDM) dm(i+1) = dm(i)/1.2_RP
   ENDDO

   DO j=1,ntot   ! For each datum, choose best derivative estimate
      best = 1.e10_RP
      ibest = 1
      DO i=1,NDM-2
         IF(ABS(dRdm(ip,j,i+0)) < 1.0e-7_RP .OR.  &
            ABS(dRdm(ip,j,i+1)) < 1.0e-7_RP .OR.  &
            ABS(dRdm(ip,j,i+2)) < 1.0e-7_RP)THEN
            test = 1.e20_RP
         ELSE
            test = ABS((dRdm(ip,j,i+0)/dRdm(ip,j,i+1) + &
                        dRdm(ip,j,i+1)/dRdm(ip,j,i+2))/2._RP-1.0_RP)
         ENDIF
         IF((test < best) .AND. (test > 1.e-7_RP))THEN
            best  = test
            ibest = i+1
         ENDIf
      ENDDO
      Jac(j,ip) = dRdm(ip,j,ibest)  ! Best deriv into Jacobian
      IF (best == 1.e20_RP)THEN
        Jac(j,ip) = 0._RP
      ENDIF
      dmbest(ip,objst%k) = dm(ibest)
   ENDDO
ENDDO
do i=1,NFP   ! Scale columns of Jacobian for stability
   ilay = CEILING(REAL(i,RP)/REAL(NPL,RP))
   ipar = i-(ilay-1)*NPL
   IF(ilay > objst%k) ipar = ipar + 1
   IF(ipar /= 1)THEN
      do j=1,ntot
         Jac(j,i) = Jac(j,i)*(maxlim(ipar)-minlim(ipar))
      enddo
   ENDIF
enddo

!OPEN(UNIT=20,FILE='jacfile.txt',FORM='formatted',STATUS='REPLACE', &
!     ACTION='WRITE',RECL=1024)
!DO i = 1,ntot
!    WRITE(20+objst%k,205) Jac(i,:)
!ENDDO
!CLOSE(20+objst%k)

JactJac = MATMUL(TRANSPOSE(Jac),Jac)

Mpriorinv = 0._RP
do i=1,NFP
   Mpriorinv(i,i) = (1.0_RP/(1.0_RP/4.0_RP))**2._RP
enddo

sdevtmp = 4.6_RP
!Ctmp = JactJac
Ctmp = JactJac+ Mpriorinv
!Ctmp = JactJac/sdevtmp**2/Tstar + Mpriorinv

CALL SVDCMP(Ctmp,WW,VVtmp)

Wmax = MAXVAL(WW,1)
Winv = 0._RP
DO i=1,NFP
   IF (WW(i)/Wmax > 1.e-9_RP) THEN
      Winv(i,i) = 1._RP/WW(i)
   ELSE
      WRITE(6,*) 'Warning: Small singular value'
      WRITE(6,*) WW/Wmax
   ENDIF
ENDDO

Cov0(1:NFP,1:NFP,objst%k) = MATMUL(MATMUL(VVtmp,Winv),TRANSPOSE(VVtmp))
Ctmp = Cov0(1:NFP,1:NFP,objst%k)
CALL SVDCMP(Ctmp,WW,VVtmp)
sdevm(:,objst%k) = 0._RP
sdevm(1:NFP,objst%k) = SQRT(WW)
VV(:,:,objst%k) = 0._RP
VV(1:NFP,1:NFP,objst%k) = VVtmp

idmbest(objst%k) = 1

204 FORMAT(i7,6(f16.5))
205 FORMAT(500ES40.22)

DEALLOCATE( dat1,dat2,dRdm,p1,p2,Jac,JactJac,bufJac,Mpriorinv, &
            Ctmp,VVtmp,Winv,WW )
RETURN
END SUBROUTINE LINROT
!=======================================================================

SUBROUTINE LINROT2(objst,dat)
!!
!!  Estimates linearized model covariance matrix when dmbest already known
!!
!=======================================================================
USE MPI
USE NR
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ifq,ip,ntot
INTEGER(KIND=IB):: ilay,ipar,NFP

TYPE (objstruc),DIMENSION(2):: obj          ! Objects 
TYPE (objstruc)             :: objst        ! initial object 
TYPE (datastruc)            :: dat          ! Data

REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: p1,p2
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: dat1,dat2
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: dRdm

REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Jac,JactJac,Mpriorinv,Temp
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: bufJac,WW,Jacscale
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: VVtmp,Winv
REAL(KIND=RP)                                :: Wmax,Tstar,sdevtmp
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Ctmp

Tstar = 1._RP
ntot = NBAND*NANG
obj(1) = objst
obj(2) = objst

NFP = objst%NFP
ALLOCATE( dat1(ntot),dat2(ntot),dRdm(NFP,ntot),       &
          p1(NBAND,NANG),p2(NBAND,NANG),    &
          Jac(ntot,NFP),JactJac(NFP,NFP),bufJac(ntot),    &
          Jacscale(NFP),Mpriorinv(NFP,NFP),Temp(ntot,NFP),&
          Ctmp(NFP,NFP),VVtmp(NFP,NFP),Winv(NFP,NFP),WW(NFP) )

DO ip=1,NFP
   obj(1)%par = objst%par
   obj(1)%par(ip) = objst%par(ip)+dmbest(ip,objst%k)
   CALL CALCH(obj(1))
   CALL LOGLHOOD(obj(1),dat)
   p1 = dat%Rrep
   DO ifq=1,NBAND
     dat1((ifq-1)*NANG+1:(ifq)*NANG) = p1(ifq,:)/objst%sdpar(ifq)
   ENDDO
   obj(2)%par = objst%par
   obj(2)%par(ip) = objst%par(ip)-dmbest(ip,objst%k)
   CALL CALCH(obj(2))
   CALL LOGLHOOD(obj(2),dat)
   p2 = dat%Rrep
   DO ifq=1,NBAND
     dat2((ifq-1)*NANG+1:(ifq)*NANG) = p2(ifq,:)/objst%sdpar(ifq)
   ENDDO
   DO j=1,ntot
     IF(ABS((dat1(j)-dat2(j))/(dat1(j)+dat2(j))) > 1.0e-7_RP)THEN
       dRdm(ip,j) = (dat1(j)-dat2(j))/(2._RP*dmbest(ip,objst%k))
     ELSE
       dRdm(ip,j) = 0._RP
     ENDIF
     Jac(j,ip) = dRdm(ip,j)  ! Best deriv into Jacobian
     IF(ABS(dRdm(ip,j)) < 1.0e-7_RP)THEN
       Jac(j,ip) = 0._RP
     ENDIF
   ENDDO
ENDDO
DO i=1,NFP   ! Scale columns of Jacobian for stability
  ilay = CEILING(REAL(i,RP)/REAL(NPL,RP))
  ipar = i-(ilay-1)*NPL
  IF(ilay > objst%k) ipar = ipar + 1
  IF(ipar /= 1)THEN
    DO j=1,ntot
      Jac(j,i) = Jac(j,i)*(maxlim(ipar)-minlim(ipar))
    ENDDO
  ENDIF
ENDDO
!OPEN(UNIT=20,FILE='jacfile2.txt',FORM='formatted',STATUS='REPLACE', &
!     ACTION='WRITE',RECL=1024)
!DO j = 1,ntot
!    WRITE(20,205) Jac(j,:)
!ENDDO
!CLOSE(20)
JactJac = MATMUL(TRANSPOSE(Jac),Jac)

Mpriorinv = 0._RP
do j=1,NFP
   Mpriorinv(j,j) = (1.0_RP/(1.0_RP/4.0_RP))**2._RP
enddo

sdevtmp = 4.6_RP
!Ctmp = JactJac
Ctmp = JactJac+ Mpriorinv
!Ctmp = JactJac/sdevtmp**2/Tstar + Mpriorinv

CALL SVDCMP(Ctmp,WW,VVtmp)

Wmax = MAXVAL(WW,1)
Winv = 0._RP
DO j=1,NFP
   IF (WW(j)/Wmax > 1.e-9_RP) THEN
      Winv(j,j) = 1._RP/WW(j)
   ELSE
      WRITE(6,*) 'Warning: Small singular value'
      WRITE(6,*) WW/Wmax
   ENDIF
ENDDO

Cov0(1:NFP,1:NFP,objst%k) = MATMUL(MATMUL(VVtmp,Winv),TRANSPOSE(VVtmp))
Ctmp = Cov0(1:NFP,1:NFP,objst%k)
CALL SVDCMP(Ctmp,WW,VVtmp)
sdevm(:,objst%k) = 0._RP
sdevm(1:NFP,objst%k) = SQRT(WW)
VV(:,:,objst%k) = 0._RP
VV(1:NFP,1:NFP,objst%k) = VVtmp

idmbest(objst%k) = 1

204 FORMAT(i7,6(f16.5))
205 FORMAT(500ES40.22)

DEALLOCATE( dat1,dat2,dRdm,p1,p2,Jac,JactJac,bufJac,Mpriorinv, &
            Ctmp,VVtmp,Winv,WW )
RETURN
END SUBROUTINE LINROT2
!=======================================================================

SUBROUTINE PRINTVV()
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i,j,l

WRITE(6,*) 'Rotation matrices:'
DO l=3,7
WRITE(6,*) 'k=',l
DO i=1,l*NPL+3
   WRITE(6,201) VV(i,:,l)
ENDDO
WRITE(6,*) 'sdevm:'
WRITE(6,201) sdevm(:,l)
WRITE(6,*) ''
ENDDO
WRITE(6,*) 'Cov0 matrices:'
DO l=3,7
WRITE(6,*) 'k=',l
DO i=1,l*NPL+3
   WRITE(6,201) Cov0(i,:,l)
ENDDO
ENDDO

201 FORMAT(100ES16.8)
END SUBROUTINE PRINTVV
!=======================================================================

SUBROUTINE SOMMER_INTEGRAND(thd,freq,cw,rw,z_t,nfrq,NANG)
!
!     Compute spherical wave reflection coeff
!     CANNOT HANDLE HALFSPACE AT THIS POINT
!     Based on code by Charles W. Holland and John Camin ARL PSU
!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: kf,i,j,imx
INTEGER(KIND=IB) :: nfrq, NANG, Nmx, Nrmx, Nr_fr
REAL(KIND=RP) :: cw,rw,z_t,w,exp_arg,umax
REAL(KIND=RP)    :: dbesj0
REAL(KIND=RP)    :: ASINH
REAL(KIND=RP),Dimension(:),ALLOCATABLE :: m,mi,k
COMPLEX(KIND=RP),Dimension(:),ALLOCATABLE :: tmp,tmpi
COMPLEX(KIND=RP),Dimension(:),ALLOCATABLE :: exp_termR
COMPLEX(KIND=RP),Dimension(:,:),ALLOCATABLE :: tmp2
COMPLEX(KIND=RP),Dimension(:,:),ALLOCATABLE :: kr_sinRTh
REAL(KIND=RP), DIMENSION(nfrq) :: freq
REAL(KIND=RP), DIMENSION(NANG) :: thd,x,r,R1
REAL(KIND=RP), DIMENSION(2,NANG) :: rr
COMPLEX(KIND=RP),DIMENSION(NANG) :: rRP,iRP
COMPLEX(KIND=RP) :: CACOS,CTAN

!Spatial Variables
x = z_t/TAN(thd*PI/180._RP)        !offset
r = x                              !xy plane, 2D range
R1 = SQRT(x*x + z_t*z_t) !Radial distance of Reflected Path

rr(1,:)= r
rr(2,:)= 0._RP

ALLOCATE( N(nfrq),Nr(nfrq),k(nfrq),drTh(nfrq),diTh(nfrq),RPs(NANG,nfrq) )
N = 0_IB
Nr = 0_IB
drTh = 0._RP
diTh = 0._RP
RPs = 0._RP
DO kf=1,nfrq

   w = 2._RP*PI*freq(kf)
   k(kf) = w/cw
   RPs(:,kf) = EXP(CMPLX(0._RP,1._RP,RP)*k(kf)*R1)/R1
   RPs(:,kf) = RPs(:,kf)/k(kf)

   ! Set limit for integration along imag axis where exponential goes to 
   ! 1e-4 of its value at pi/2;
   exp_arg=9.2_RP/(k(kf)*z_t)
   umax=1.1_RP*ASINH(exp_arg)

   !Fix N from analytic considerations, 5 samples per num oscillations 0->pi/2
   Nr_fr=FLOOR(k(kf)*MAXVAL(r+z_t)/(2._RP*PI))
   Nr(kf) = MAXVAL((/Nr_fr*5_IB,80_IB/))
   IF(MOD(REAL(Nr(kf),RP),2._RP) == 0._RP)THEN
      Nr(kf) = Nr(kf)+1_IB
   ENDIF

   ! Theta Spacing in Complex Plane...
   ! int_end= pi/2-.001i;
   ! int_end= pi/2;rTh = linspace(0,int_end,N);
   ! iTh = linspace(int_end,pi/2-i*umax,N);
   ! drTh = mean(diff(rTh));diTh = mean(diff(iTh));

   drTh(kf) = PI/2._RP/REAL(Nr(kf)-1,RP)
   diTh(kf) = -CMPLX(0._RP,1._RP,RP)*drTh(kf)
   N(kf) = Nr(kf)+FLOOR(umax/drTh(kf))
   ! force the total to also be odd then the int along imag will also 
   ! be odd N-Nr+1
   IF(MOD(N(kf),2)==0)THEN
      N(kf) = N(kf)-1_IB
   ENDIF

ENDDO

!!
!! Find matrix dimensions:
!!
Nmx = MAXVAL(N)
!IF(rank==src)PRINT*,N

DO kf=1,nfrq
   ALLOCATE(Sarg(kf)%btR2(N(kf),NANG),Sarg(kf)%rTh2(N(kf)))
   Sarg(kf)%btR2 = CMPLX(0._RP,0._RP,RP)
   Sarg(kf)%rTh2 = CMPLX(0._RP,0._RP,RP)
ENDDO

DO kf=1,nfrq
   ALLOCATE( tmp2(N(kf),2),exp_termR(N(kf)),kr_sinRTh(N(kf),NANG) )
   tmp2 = 0._RP
   exp_termR = 0._RP
   kr_sinRTh = 0._RP
   Sarg(kf)%rTh2(1:Nr(kf)) = (/ ((i*drTh(kf))-drTh(kf) ,i=1,Nr(kf)) /)
   Sarg(kf)%rTh2(Nr(kf)+1:N(kf)) = PI/2._RP-CMPLX(0._RP,1._RP,RP)*(/ (i*drTh(kf),i=1,N(kf)-Nr(kf)) /)
   Sarg(kf)%rTh2(1:N(kf)) = Sarg(kf)%rTh2(1:N(kf)) + TINY(1._RP)

   !Numerical Integration masking scheme
   ALLOCATE( m(Nr(kf)),mi(N(kf)-Nr(kf)+1),tmp(Nr(kf)),tmpi(N(kf)-Nr(kf)+1) )
   tmp = 0._RP
   tmpi = 0._RP
   m           = 1._RP      !Masking array for Simpsons Rule
   m(2:Nr(kf)-1)   = 4._RP      !(1,4,2,4,2,4,1)
   m(3:Nr(kf)-1:2) = 2._RP

   mi = 1._RP               !Masking array for Simpsons Rule
   mi(2:N(kf)-Nr(kf)-1+1)   = 4._RP !(1,4,2,4,2,4,1)
   mi(3:N(kf)-Nr(kf)-1+1:2) = 2._RP

   tmp2(:,1) = SIN(Sarg(kf)%rTh2(1:N(kf)))
   tmp2(:,2) = 0._RP            !! Had to add second line here to allow for use of MATMUL...
   kr_sinRTh=k(kf)*MATMUL(tmp2,rr)

   ! Bessel function first kind, zeroth order (far field approx. and then exact near field)
   Sarg(kf)%btR2(1:N(kf),:) = SQRT(2._RP/(PI*kr_sinRTh))*COS(kr_sinRTh-0.25_RP*PI) !large arg approx >4;
   !!
   !!  ??????BESSEL FUNCTION FOR COMPLEX VARIABLE??????
   !!
   DO i = 1,NANG
   DO j = 1,N(kf)
      IF(ABS(kr_sinRTh(j,i)) < 4._RP)THEN 
         Sarg(kf)%btR2(j,i)=dbesj0(REAL(kr_sinRTh(j,i),RP))
      ENDIF
   ENDDO
   ENDDO

   exp_termR=EXP(CMPLX(0._RP,1._RP,RP)*k(kf)*(z_t)*COS(Sarg(kf)%rTh2(1:N(kf))))*SIN(Sarg(kf)%rTh2(1:N(kf)))
   tmp = m*exp_termR(1:Nr(kf))
   tmpi = mi*exp_termR(Nr(kf):N(kf))

   DO i = 1,NANG
      Sarg(kf)%btR2(1:Nr(kf),i) = Sarg(kf)%btR2(1:Nr(kf),i)*tmp
      Sarg(kf)%btR2(Nr(kf):N(kf),i) = Sarg(kf)%btR2(Nr(kf):N(kf),i)*tmpi
   ENDDO

   DEALLOCATE( tmp2,m,mi,tmp,tmpi,kr_sinRTh,exp_termR )
ENDDO

207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SOMMER_INTEGRAND
!=======================================================================

SUBROUTINE SPH_REF_NLAY(thd,m_rg,freq,ref_sph,cw,rw,z_t,nfrq,NANG,NFP,idx)
!
!     Compute spherical wave reflection coeff
!     CANNOT HANDLE HALFSPACE AT THIS POINT
!     Based on code by Charles W. Holland and John Camin ARL PSU
!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: kf,lf,i,j,l,idx,Ni
INTEGER(KIND=IB) :: nfrq, NANG,NFP,N_tmp
REAL(KIND=RP) :: cw,rw,z_t
REAL(KIND=RP)    :: dbesj0
REAL(KIND=RP)    :: ASINH
REAL(KIND=RP), DIMENSION(nfrq) :: freq
REAL(KIND=RP), DIMENSION(NFP)  :: m_rg
REAL(KIND=RP), DIMENSION(NANG) :: thd
COMPLEX(KIND=RP),DIMENSION(NANG) :: rRP,iRP
REAL(KIND=RP),DIMENSION(nfrq,NANG) :: ref_sph
COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE :: ref,ref2,rTh_tmp
COMPLEX(KIND=RP) :: CACOS,CTAN

DO kf=1,nfrq
   lf = (idx-1)*nfrq+kf
!   PRINT*,'lf=',lf
   Ni = N(lf)-Nr(lf) + 4
   !!
   !! Only subsample real part!! Leave imaginary angles alone.
   !!
   IF(MOD(N(lf)-Ni,subsmp) == 0._RP)THEN
      N_tmp = CEILING(REAL(N(lf)-Ni-1,RP)/REAL(subsmp,RP))
   ELSE
      N_tmp = CEILING(REAL(N(lf)-Ni-1,RP)/REAL(subsmp,RP))+1
   ENDIF
   N_tmp = N_tmp + Ni
   ALLOCATE( ref(N_tmp),rTh_tmp(N_tmp),ref2(N(lf)) )
   IF(subsmp > 1)THEN
      ref = 0._RP
      ref2= 0._RP
      rTh_tmp = 0._RP
      j = 1
      DO i=1,N(lf)-Ni,subsmp
         rTh_tmp(j) = Sarg(lf)%rTh2(i)
         j = j + 1
      ENDDO
      rTh_tmp(N_tmp-Ni:N_tmp) = Sarg(lf)%rTh2(N(lf)-Ni:N(lf))

      CALL REF_NLAY4(90._RP-rTh_tmp*180._RP/PI,m_rg,freq(kf),ref,cw,rw,1,N_tmp,NFP)

      l = 1
      DO i=1,N_tmp-1-Ni
         DO j=1,subsmp
            ref2(l) = ref(i)+(ref(i+1)-ref(i))*REAL(j-1,RP)/REAL(subsmp,RP)
            l = l + 1
         ENDDO
      ENDDO
      ref2(N(lf)-Ni:N(lf)) = ref(N_tmp-Ni:N_tmp)
   ELSE
      ref2 = 0._RP
      CALL REF_NLAY4(90._RP-Sarg(lf)%rTh2(1:N(lf))*180._RP/PI,m_rg,freq(kf),ref2,cw,rw,1,N(lf),NFP)
   ENDIF

   rRP = drTh(lf)/3._RP*MATMUL(ref2(1:Nr(lf)),Sarg(lf)%btR2(1:Nr(lf),:))
   iRP = diTh(lf)/3._RP*MATMUL(ref2(Nr(lf):N(lf)),Sarg(lf)%btR2(Nr(lf):N(lf),:))

   ref_sph(kf,:) = ABS(CMPLX(0._RP,1._RP,RP)*(rRP+iRP)/RPs(:,lf))
   DEALLOCATE( ref,ref2,rTh_tmp )
ENDDO
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SPH_REF_NLAY
!=======================================================================

SUBROUTINE REF_NLAY4(thd,m_rg,freq,ref,cw,rw,nfrq,NANG,NFP)
!
!     Compute plane-wave refl coeff
!     CANNOT HANDLE HALFSPACE AT THIS POINT
!     Based on code by Charles W. Holland ARL PSU
!
!=======================================================================
USE DATA_TYPE
!USE MPI
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr, iang, il
INTEGER(KIND=IB) :: nfrq, NLAY, NLAY2, NANG, NFP
REAL(KIND=RP) :: cw,rw,dB2nep,pi180
REAL(KIND=RP),DIMENSION(nfrq)            :: freq
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: freq2
REAL(KIND=RP),DIMENSION(NFP)  :: m_rg
COMPLEX(KIND=RP), DIMENSION(NANG) :: thd,th1,sinth1,costh1
COMPLEX(KIND=RP), DIMENSION(nfrq,NANG) :: ref
COMPLEX(KIND=RP), DIMENSION(NANG,nfrq) :: reftmp,z1
REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: c,alf,r,rc,d
COMPLEX(KIND=RP), DIMENSION(:), ALLOCATABLE :: v,rv,vc1
COMPLEX(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: zz,th,sinth,k
COMPLEX(KIND=RP), DIMENSION(NANG,nfrq) :: znlay,znlay1,zin
COMPLEX(KIND=RP), DIMENSION(NANG,nfrq) :: tankd
COMPLEX(KIND=RP) :: CACOS,CTAN
LOGICAL :: ISNAN

NLAY = ((NFP-3)/4)+1
NLAY2 = NLAY+1

ALLOCATE(c(NLAY2),alf(NLAY2),r(NLAY2),rc(NLAY2),d(NLAY-1))
ALLOCATE(v(NLAY2),vc1(NLAY2),rv(NLAY2),zz(NANG,NLAY2),th(NANG,NLAY2),sinth(NANG,NLAY2))
ALLOCATE(k(nfrq,NLAY2))

c = 0.0_RP; alf = 0.0_RP; r = 0.0_RP; d = 0._RP

c   = (/cw, m_rg((/2:NFP-3:4/)), m_rg(NFP-2)/)
alf = (/0._RP, m_rg((/4:NFP-3:4/)), m_rg(NFP) /)
r   = (/rw, m_rg((/3:NFP-3:4/)), m_rg(NFP-1)/)
d   = m_rg((/1:NFP-3:4/))

ref = 0._RP
reftmp = 0._RP
dB2nep=2._RP*PI*20._RP*1000._RP/LOG(10._RP)

!
! DOUBLE CHECK THIS
!
v=1._RP/CMPLX(1._RP/c,alf/dB2nep,RP) ! force radiation condition to be satisfied

th = 0._RP
sinth = 0._RP
zz = 0._RP

pi180 = PI/180._RP
rc = r*c
rv = r*v
vc1 = v/c(1)
th1 = thd*pi180
sinth1 = SIN(th1)
costh1 = COS(th1)

DO il = 1,NLAY2
   DO iang = 1,NANG
      IF(il == 1)THEN
         zz(iang,il) = rc(il)/sinth1(iang) !since incident angles are real so must z1
      ELSE
         th(iang,il)    = CACOS(costh1(iang)*vc1(il))
         sinth(iang,il) = SIN(th(iang,il))
         zz(iang,il)    = rv(il)/sinth(iang,il)
      ENDIF
   ENDDO
ENDDO

!  COMPUTE input impedance
DO ifr = 1,nfrq
   k(ifr,:) = 2._RP*PI*freq(ifr)/v(:)
   DO iang = 1,NANG
      IF(NLAY == 1)THEN ! HALFSPACE
         reftmp(iang,ifr)=(zz(iang,2) - zz(iang,1) )/(zz(iang,2) + zz(iang,1));
      ELSE
         tankd(iang,ifr) = CTAN( k(ifr,NLAY)*d(NLAY-1)*sinth(iang,NLAY) )
         znlay(iang,ifr)  = zz(iang,NLAY)
         znlay1(iang,ifr) = zz(iang,NLAY+1)

         zin(iang,ifr) = znlay(iang,ifr)*(znlay1(iang,ifr)-CMPLX(0._RP,1._RP,RP)*znlay(iang,ifr)*tankd(iang,ifr))/ &
                        (znlay(iang,ifr) - CMPLX(0._RP,1._RP,RP)*znlay1(iang,ifr)*tankd(iang,ifr))
         DO il = NLAY,3,-1
            tankd(iang,ifr) = CTAN( k(ifr,il-1)*d(il-2)*sinth(iang,il-1) )
            zin(iang,ifr) = zz(iang,il-1)*(zin(iang,ifr) -CMPLX(0._RP,1._RP,RP)*zz(iang,il-1)*tankd(iang,ifr))/ &
                           (zz(iang,il-1) - CMPLX(0._RP,1._RP,RP)*zin(iang,ifr)*tankd(iang,ifr));
         ENDDO
         z1(iang,ifr)=zz(iang,1)
         reftmp(iang,ifr) = (zin(iang,ifr)-z1(iang,ifr))/(zin(iang,ifr)+z1(iang,ifr))
         ! At 0 degrees, the reflection coefficeint is -1
         IF(ISNAN(ABS(reftmp(iang,ifr))) == .TRUE.) reftmp(iang,ifr)=-1._RP
      ENDIF
   ENDDO
ENDDO

ref = TRANSPOSE(reftmp)

DEALLOCATE(c,alf,r,rc,d,v,vc1,rv,zz,th,sinth,k)
RETURN
END SUBROUTINE REF_NLAY4

!=======================================================================

SUBROUTINE REF_NLAY3(thd,m_rg,freq,ref,cw,rw,nfrq,NANG,NFP)
!
!     Compute plane-wave refl coeff
!     CANNOT HANDLE HALFSPACE AT THIS POINT
!     Based on code by Charles W. Holland ARL PSU
!
!=======================================================================
USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB) :: i, j, l
INTEGER(KIND=IB) :: nfrq, NLAY, NLAY2, NANG, NFP
REAL(KIND=RP) :: cw,rw,dB2nep
REAL(KIND=RP), DIMENSION(nfrq) :: freq
REAL(KIND=RP), DIMENSION(NFP)  :: m_rg
REAL(KIND=RP), DIMENSION(NANG) :: thd
REAL(KIND=RP), DIMENSION(nfrq,NANG) :: ref
REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: c,alf,r,d,th1
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: freq2,z1
COMPLEX(KIND=RP), DIMENSION(:), ALLOCATABLE :: v,th1tmp
COMPLEX(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: zz,th,v2,k,reftmp
COMPLEX(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: znlay,znlay1,zin,tankd,tankdtmp
COMPLEX(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: thtmp,ktmp,zm1
COMPLEX(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: z
COMPLEX(KIND=RP) :: CACOS,CTAN

NLAY = ((SIZE(m_rg,1)-3)/4)+1
NLAY2 = NLAY+1

ALLOCATE(c(NLAY2),alf(NLAY2),r(NLAY2),d(NLAY-1),th1(NANG))
ALLOCATE(th1tmp(NANG))
ALLOCATE(v(NLAY2),v2(NLAY2,nfrq),zz(NLAY2,NANG),th(NLAY2,NANG))
ALLOCATE(z(NLAY2,NANG,nfrq),freq2(NLAY2,nfrq),k(nfrq,NLAY2))
ALLOCATE(reftmp(NANG,nfrq),znlay(NANG,nfrq),znlay1(NANG,nfrq))
ALLOCATE(zin(NANG,nfrq),zm1(NANG,nfrq),z1(NANG,nfrq))
ALLOCATE(tankdtmp(nfrq,NANG),tankd(nfrq,NANG),thtmp(nfrq,NANG),ktmp(nfrq,NANG))

c = 0.0_RP; alf = 0.0_RP; r = 0.0_RP; d = 0._RP

c   = (/cw, m_rg((/2:NFP-3:4/)), m_rg(NFP-2)/)
alf = (/0._RP, m_rg((/4:NFP-3:4/)), m_rg(NFP) /)
r   = (/rw, m_rg((/3:NFP-3:4/)), m_rg(NFP-1)/)
d   = m_rg((/1:NFP-3:4/))

ref = 0.0_RP
reftmp = 0.0_RP
dB2nep=2._RP*PI*20._RP*1000._RP/LOG(10._RP)

!
! DOUBLE CHECK THIS
!
v=1._RP/CMPLX(1._RP/c,alf/dB2nep,RP) ! force radiation condition to be satisfied
th1=thd*PI/180._RP
th = 0._RP
zz = 0._RP
zz(1,:) = r(1)*c(1)/SIN(th1) !since incident angles are real so must z1

DO j = 2,SIZE(c,1)
   DO i = 1,NANG
      th1tmp(i) = CACOS(COS(th1(i))*v(j)/c(1))
   ENDDO
   th(j,:)=th1tmp
   zz(j,:)=r(j)*v(j)/SIN(th(j,:))
ENDDO

DO i = 1,nfrq
   z(:,:,i) = zz
   v2(:,i) = v
ENDDO

DO i = 1,NLAY2
   freq2(i,:) = freq
ENDDO

k = TRANSPOSE(2._RP*PI*freq2/v2)

IF(NLAY == 1)THEN ! HALFSPACE

    reftmp=(z(2,:,:) - z(1,:,:) )/(z(2,:,:) + z(1,:,:));

ELSE

   !  COMPUTE input impedance
   DO i = 1,nfrq
      thtmp(i,:) = SIN(th(NLAY,:))
   ENDDO
   DO i = 1,NANG
      ktmp(:,i) = k(:,NLAY)
   ENDDO
!   tankd = CTAN( ktmp*d(NLAY-1)*thtmp )
   tankdtmp = ktmp*d(NLAY-1)*thtmp
   DO i = 1,nfrq
      DO j = 1,NANG
         tankd(i,j) = CTAN( tankdtmp(i,j) )
      ENDDO
   ENDDO

   znlay  = z(NLAY,:,:)
   znlay1 = z(NLAY+1,:,:)

   zin = znlay*(znlay1-CMPLX(0._RP,1._RP,RP)*znlay*TRANSPOSE(tankd))/ &
         (znlay - CMPLX(0._RP,1._RP,RP)*znlay1*TRANSPOSE(tankd))

   DO i = NLAY,3,-1
      DO j = 1,nfrq
         thtmp(j,:) = SIN(th(i-1,:))
      ENDDO
      DO j = 1,NANG
         ktmp(:,j) = k(:,i-1)
      ENDDO
!      tankd = CTAN( ktmp*d(i-2)*thtmp )
      tankdtmp = ktmp*d(i-2)*thtmp
      DO l = 1,nfrq
         DO j = 1,NANG
            tankd(l,j) = CTAN( tankdtmp(l,j) )
         ENDDO
      ENDDO
      zm1 = z(i-1,:,:)

      zin = zm1*(zin -CMPLX(0._RP,1._RP,RP)*zm1*TRANSPOSE(tankd))/ &
           (zm1 - CMPLX(0._RP,1._RP,RP)*zin*TRANSPOSE(tankd));
   ENDDO
   z1=z(1,:,:)

   reftmp = (zin-z1)/(zin+z1)
ENDIF

ref = ABS(TRANSPOSE(reftmp))

DEALLOCATE(c,alf,r,d,th1,v,v2,zz,th,z,freq2,k,reftmp,znlay, &
           znlay1,zin,zm1,z1,tankd,thtmp,ktmp)

RETURN
END SUBROUTINE REF_NLAY3
!=======================================================================

SUBROUTINE KK_WATERS_dBmkHz(y,fr,alfrdB,cr,f1,N,alf1dB,vp)
!=======================================================================
!function [alf1dB vp]=KK_Waters_dBmkHz(n,fr,alfrdB,cr,f1)
!! vp=KK_Waters(n,fr,alfr,cr,f1)
!!
!! %Kramers Kronig dispersion from Waters et al. (JASA 108, pg 556)
!!
!! INPUTS
!! n = frequency exponent of attenuation assuming power law
!! fr = reference frequency for attenuation
!! alfr = attenuation in dB/m/kHz at reference frequency fr
!! cr = sound speed at reference frequency, fr
!! f1 = frequency of interest (can be a vector)
!!
!! OUTPUTS
!! alf1dB = attenuation in (dB/m/kHz) for frequency vector f1
!! vp = compressional velocity in (m/s)
USE DATA_TYPE
REAL(KIND=RP) :: y,fr,alfrdB,cr,wo
REAL(KIND=RP),DIMENSION(N) :: f1,alf1dB,vp,alfr,w1,Q2

alfr = alfrdB*(fr/1000._RP)**y/(20._RP*LOG10(EXP(1._RP))) !reference attenuation in nepers/m
wo = 2._RP*PI*fr
w1 = 2._RP*PI*f1

IF(y == 1._RP)THEN
 Q2 = -2._RP/PI*LOG(f1/fr)
ELSE
  p = y-1
 Q2 = (w1**p-wo**p)*TAN(y*PI/2._RP)
ENDIF

IF(y > 3)THEN
   PRINT*,'y must be less than 3!'
!   STOP
ELSE
    vp=1._RP/(1._RP/cr + Q2*alfr/(2._RP*PI*fr)**y)
    alf1dB= alfrdB*(f1/fr)**(y-1._RP)
ENDIF

RETURN
END SUBROUTINE KK_WATERS_dBmkHz
!=======================================================================

SUBROUTINE SND_PARTICLE(obj,ipart,idest,idone,logLtrgt,betatmp,NTEMP1,NBAL1,irealloc)
!=======================================================================
!!
!! Exchanging particles
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB):: ipart,idest,idone,NTEMP1,NBAL1,irealloc
REAL(KIND=RP):: t1,t2
TYPE(objstruc):: obj
REAL(KIND=RP),DIMENSION(len_snd2):: obj_tmp
REAL(KIND=RP),DIMENSION(NTEMP1)  :: betatmp
REAL(KIND=RP)                    :: logLtrgt

!!
!!  Sending particle to slave
!!
IF(rank == src)THEN
   CALL MPI_SSEND( idone, 1, MPI_INTEGER, idest, &
                   ipart, MPI_COMM_WORLD, ierr )
   IF(idone == 0)THEN
      CALL MPI_SSEND( irealloc, 1, MPI_INTEGER, idest, &
                      ipart, MPI_COMM_WORLD, ierr )
      IF(irealloc == 1)THEN
         CALL MPI_SSEND( NTEMP1, 1, MPI_INTEGER, idest, &
                         ipart, MPI_COMM_WORLD, ierr )
         CALL MPI_SSEND( NBAL1, 1, MPI_INTEGER, idest, &
                         ipart, MPI_COMM_WORLD, ierr )
         CALL MPI_SSEND( betatmp, NTEMP1, MPI_DOUBLE_PRECISION, idest, &
                         ipart, MPI_COMM_WORLD, ierr )
      ENDIF
      CALL MPI_SSEND( logLtrgt, 1, MPI_DOUBLE_PRECISION, idest, &
                      ipart, MPI_COMM_WORLD, ierr )
      obj_tmp = 0._RP
      obj_tmp =  (/ obj%logL, obj%logwt, obj%logPr, obj%lognorm, REAL(obj%k,RP), REAL(obj%NFP,RP), &
                    obj%par, obj%sdpar, obj%z, obj%g, obj%gp, obj%detChat, &
                    obj%Chat(1,:), obj%Chat(2,:), obj%Chat(3,:), &
                    obj%Chati(1,:),obj%Chati(2,:),obj%Chati(3,:),REAL(obj%nfail,RP) /)
      CALL MPI_SSEND( obj_tmp, len_snd2, MPI_DOUBLE_PRECISION, idest, &
                      ipart, MPI_COMM_WORLD, ierr )
   ENDIF
ELSE
   CALL MPI_RECV(idone, 1, MPI_INTEGER, src,&
                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
   IF(idone == 0)THEN
      obj_tmp = 0._RP
      CALL MPI_RECV(irealloc, 1, MPI_INTEGER, src,&
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
      IF(irealloc == 1)THEN
         CALL MPI_RECV(NTEMP1, 1, MPI_INTEGER, src,&
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
         CALL MPI_RECV(NBAL1, 1, MPI_INTEGER, src,&
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
         DEALLOCATE(beta)
         ALLOCATE(beta(NTEMP1))
         CALL MPI_RECV(beta, NTEMP1, MPI_DOUBLE_PRECISION, src,&
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
      ENDIF
      CALL MPI_RECV(logLtrgt, 1, MPI_DOUBLE_PRECISION, src,&
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV(obj_tmp, len_snd2, MPI_DOUBLE_PRECISION, src,&
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
      ipart = status(MPI_TAG)  !! Particla number is passed along as tag
      obj%logL       = obj_tmp(1)
      obj%logwt      = obj_tmp(2)
      obj%logPr      = obj_tmp(3)
      obj%lognorm    = obj_tmp(4)
      obj%k          = INT(obj_tmp(5))
      obj%NFP        = INT(obj_tmp(6))
      obj%par        = 0._RP
      obj%par        = obj_tmp(7:7+NFPMX-1)
      obj%sdpar      = 0._RP
      obj%sdpar      = obj_tmp(7+NFPMX:7+NFPMX+NBAND-1)
      obj%z          = 0._RP
      obj%z          = obj_tmp(7+NFPMX+NBAND:7+NFPMX+NBAND+NLMX-1)
      obj%g          = obj_tmp(7+NFPMX+NBAND+NLMX:7+NFPMX+NBAND+NLMX+3-1)
      obj%gp         = obj_tmp(7+NFPMX+NBAND+NLMX+3:7+NFPMX+NBAND+NLMX+6-1)
      obj%detChat    = obj_tmp(7+NFPMX+NBAND+NLMX+6)
      obj%Chat(1,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1:7+NFPMX+NBAND+NLMX+6+1+3-1)
      obj%Chat(2,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+3:7+NFPMX+NBAND+NLMX+6+1+6-1)
      obj%Chat(3,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+6:7+NFPMX+NBAND+NLMX+6+1+9-1)
      obj%Chati(1,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+9:7+NFPMX+NBAND+NLMX+6+1+12-1)
      obj%Chati(2,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+12:7+NFPMX+NBAND+NLMX+6+1+15-1)
      obj%Chati(3,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+15:7+NFPMX+NBAND+NLMX+6+1+18-1)
      obj%nfail      = INT(obj_tmp(7+NFPMX+NBAND+NLMX+6+1+18))
   ENDIF
ENDIF !  MPI
CALL FLUSH(6)
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SND_PARTICLE
!=======================================================================

SUBROUTINE RCV_PARTICLE(obj,isource,ipart,ifail,tcmp)
!=======================================================================
!!
!! Exchanging particles
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB):: isource,ipart,ifail
REAL(KIND=RP):: logLG,tcmp
REAL(KIND=RP):: LOGPLUS
REAL(KIND=RP):: t1,t2
TYPE(objstruc):: obj
REAL(KIND=RP),DIMENSION(len_snd):: obj_tmp

!!
!!  Sending samples to master
!!
IF(rank == src)THEN
   obj_tmp = 0._RP
   CALL MPI_RECV(obj_tmp, len_snd, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
   isource = status(MPI_SOURCE)
   ipart   = status(MPI_TAG)
   CALL MPI_RECV(tcmp, 1, MPI_DOUBLE_PRECISION, isource,&
                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
   CALL MPI_RECV(ifail, 1, MPI_INTEGER, isource,&
                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
   obj%logL       = obj_tmp(1)
   obj%logwt      = obj_tmp(2)
   obj%logPr      = obj_tmp(3)
   obj%lognorm    = obj_tmp(4)
   obj%k          = INT(obj_tmp(5))
   obj%NFP        = INT(obj_tmp(6))
   obj%par        = 0._RP
   obj%par        = obj_tmp(7:7+NFPMX-1)
   obj%sdpar      = 0._RP
   obj%sdpar      = obj_tmp(7+NFPMX:7+NFPMX+NBAND-1)
   obj%z          = 0._RP
   obj%z          = obj_tmp(7+NFPMX+NBAND:7+NFPMX+NBAND+NLMX-1)
   obj%g          = obj_tmp(7+NFPMX+NBAND+NLMX:7+NFPMX+NBAND+NLMX+3-1)
   obj%gp         = obj_tmp(7+NFPMX+NBAND+NLMX+3:7+NFPMX+NBAND+NLMX+6-1)
   obj%detChat    = obj_tmp(7+NFPMX+NBAND+NLMX+6)
   obj%Chat(1,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1:7+NFPMX+NBAND+NLMX+6+1+3-1)
   obj%Chat(2,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+3:7+NFPMX+NBAND+NLMX+6+1+6-1)
   obj%Chat(3,:)  = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+6:7+NFPMX+NBAND+NLMX+6+1+9-1)
   obj%Chati(1,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+9:7+NFPMX+NBAND+NLMX+6+1+12-1)
   obj%Chati(2,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+12:7+NFPMX+NBAND+NLMX+6+1+15-1)
   obj%Chati(3,:) = obj_tmp(7+NFPMX+NBAND+NLMX+6+1+15:7+NFPMX+NBAND+NLMX+6+1+18-1)
   acc_ratio      = obj_tmp(7+NFPMX+NBAND+NLMX+25)
   iaccept_bd     = INT(obj_tmp(7+NFPMX+NBAND+NLMX+26))
   ireject_bd     = INT(obj_tmp(7+NFPMX+NBAND+NLMX+27))
   i_bd           = INT(obj_tmp(7+NFPMX+NBAND+NLMX+28))
   i_zpert        = INT(obj_tmp(7+NFPMX+NBAND+NLMX+29))
   obj%nfail      = INT(obj_tmp(7+NFPMX+NBAND+NLMX+30))
ELSE
   obj_tmp = 0._RP
   obj_tmp =  (/ obj%logL, obj%logwt, obj%logPr, obj%lognorm, REAL(obj%k,RP), REAL(obj%NFP,RP), &
                 obj%par , obj%sdpar , obj%z, obj%g, obj%gp, obj%detChat,  &
                 obj%Chat(1,:), obj%Chat(2,:), obj%Chat(3,:),  &
                 obj%Chati(1,:),obj%Chati(2,:),obj%Chati(3,:), &
                 REAL(iaccept,RP)/REAL(ireject,RP),REAL(iaccept_bd,RP),&
                 REAL(ireject_bd,RP),REAL(i_bd,RP),REAL(i_zpert,RP),REAL(rank,RP),&
                 REAL(obj%nfail,RP) /)
   CALL MPI_SSEND( obj_tmp, len_snd, MPI_DOUBLE_PRECISION, src, &
                   ipart, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( tcmp, 1, MPI_DOUBLE_PRECISION, src, &
                   ipart, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( ifail, 1, MPI_INTEGER, src, &
                   ipart, MPI_COMM_WORLD, ierr )
ENDIF !  MPI
CALL FLUSH(6)
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE RCV_PARTICLE
!=======================================================================

SUBROUTINE SND_VV()
!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
REAL(KIND=RP),DIMENSION(NVV)   :: buf_snd_VV      ! Buffers for MPI sending rotation matrix VV
REAL(KIND=RP),DIMENSION(NSDEVM):: buf_snd_sdevm   ! Buffers for MPI sending step size

!!
!!  BCAST from master to all slaves
!!
buf_snd_VV = 0._RP
buf_snd_VV = PACK(VV(:,:,:),.true.)
CALL MPI_BCAST( buf_snd_VV,NVV, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierr )
buf_snd_sdevm = 0._RP
buf_snd_sdevm = PACK(sdevm(:,:),.true.)
CALL MPI_BCAST( buf_snd_sdevm,NSDEVM, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierr )

CALL MPI_BCAST( idmbest,NLMX, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )

IF(rank /= src)THEN
   VV = 0._RP
   VV = RESHAPE(buf_snd_VV,(/ (NLMX*NPL)+NPL-1,(NLMX*NPL)+NPL-1,NLMX /))
   sdevm = 0._RP
   sdevm = RESHAPE(buf_snd_sdevm,(/ (NLMX*NPL)+NPL-1,NLMX /))
ENDIF
RETURN
END SUBROUTINE SND_VV
!=======================================================================

SUBROUTINE COUPLE_CR(obj)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: i,j,ncra
TYPE (objstruc)                             :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+1)          :: idxc,idxr
REAL(KIND=RP),DIMENSION(kmax)               :: dc,dr
idxcra= 0
idxc  = 0
idxr  = 0
CALL GETIDXCRA(obj,idxcra,ncra)
idxc = idxcra(1:ncra:(NPL-1))
idxr = idxcra(2:ncra:(NPL-1))
dc = 0
dr = 0
dc(1:obj%k) = obj%par(idxc(2:obj%k+1)) - obj%par(idxc(1:obj%k))
dr(1:obj%k) = obj%par(idxr(2:obj%k+1)) - obj%par(idxr(1:obj%k))
!!
!! ioutside must not be set to 0 here!! Would destroy the 
!! perturbation setting from earlier (in PROPOSAL)
!!
DO i=1,obj%k
   IF(dc(i) >= 0._RP)THEN
      IF(dr(i) <  0._RP) ioutside = 1
   ELSE
      IF(dr(i) >= 0._RP) ioutside = 1
   ENDIF
ENDDO
128 FORMAT(20F16.6)
END SUBROUTINE COUPLE_CR
!=======================================================================

SUBROUTINE SAVEREPLICA(obj,dat,filename)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
TYPE (objstruc)  :: obj      ! Best object
TYPE (datastruc) :: dat      ! Best object
CHARACTER(len=64):: filename ! replica file name

WRITE(6,*) 'Global best model:'
CALL PRINTPAR(obj)
WRITE(6,*) 'Global best logL = ',obj%logL
OPEN(UNIT=50,FILE=filename,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
DO i = 1,NBAND
   WRITE(50,208) dat%Rrep(i,:)
ENDDO
WRITE(50,208) dat%angobs
!DO i = 1,NBAND
!   WRITE(50,208) dat%Rex(i,:)
!ENDDO
CLOSE(50)

208 FORMAT(50ES16.8)
RETURN
END SUBROUTINE SAVEREPLICA
!=======================================================================

SUBROUTINE WRITE_LOG()
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
!!------------------------------------------------------------------------
!!
!!  Print sampling parameters to screen for logging
!!
IF(rank == src)THEN
   WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(6,*) '~~~                                                        ~~~  '
   WRITE(6,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
   WRITE(6,*) '~~~                                                        ~~~  '
   WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(6,*) '...running on ',NTHREAD,' cores'
   WRITE(6,*) ''
   WRITE(6,210) 'ICOV                 :   ',ICOV
   WRITE(6,210) 'Number of angles     :   ',NANG
   WRITE(6,210) 'Number of frequencies:   ',NBAND
   CALL FLUSH(6)
   WRITE(6,209) 'Particle files:           ',particlefile
   WRITE(6,*) ''
   WRITE(6,*) 'NPART    = ',NPART
   WRITE(6,*) 'MCMC_STEP= ',MCMC_STEPS
   WRITE(6,*) 'MCMC_BI  = ',MCMC_BI
   WRITE(6,*) 'NTEMP    = ',NTEMP
   WRITE(6,*) 'NBAL     = ',NBAL 
   WRITE(6,*) 'NFPMX    = ',NFPMX
   WRITE(6,*) 'kmin     = ',kmin
   WRITE(6,*) 'kmax     = ',kmax
   WRITE(6,*) 'minlim:  '
   WRITE(6,201) minlim
   WRITE(6,*) 'maxlim:  '
   WRITE(6,201) maxlim
   WRITE(6,*) 'pertsdsc:'
   WRITE(6,201) pertsdsc
   WRITE(6,*) 'subsmp = ',subsmp
   WRITE(6,*) ''
   WRITE(6,*) 'Done reading data.'
   WRITE(6,*) ''
   CALL FLUSH(6)
   IF(IGA == 0)THEN
      WRITE(*,*)' !!!   NARROW BAND INTENSITY FREQ AVERAGING   !!!'
      WRITE(*,*)' FBW = ',FBW
   ELSEIF(IGA == 1)THEN
      WRITE(*,*)' !!!   GAUSSIAN FREQ AVERAGING   !!!'
   ELSEIF(IGA == 2)THEN
      WRITE(*,*)' !!!   1/3 octave INTENSITY FREQ AVERAGING   !!!'
   ENDIF
   WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(6,*) ''
ENDIF
RETURN
201 FORMAT(200F12.4)
209 FORMAT(A26,A64)
210 FORMAT(A26,I4)
END SUBROUTINE WRITE_LOG
!=======================================================================

SUBROUTINE MAKE_BETA(beta1,beta4,NTEMP1)
!!
!! Make cooling schedule (three geometrically spaced legs)
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: itemp,NTEMP1
REAL(KIND=RP)                  :: beta1,beta2,beta3,beta4
REAL(KIND=RP),DIMENSION(NTEMP1):: logbeta

beta2 = beta1+0.2_RP*ABS(beta1-beta4)
beta3 = beta1+0.75_RP*ABS(beta1-beta4)
!IF(rank == src)THEN
!   WRITE(6,*)'Computing cooling schedule (inverse temperature,geometrical spacing):'
!   WRITE(6,*)'beta 1 =',beta1,'1 to ',NTEMP/20
!   WRITE(6,*)'beta 2 =',beta2,NTEMP/20+1,' to ',NTEMP/4
!   WRITE(6,*)'beta 3 =',beta3,NTEMP/4+1,' to ',NTEMP
!ENDIF
logbeta = 0._RP
beta = 0._RP
logbeta(1) = LOG(beta1)
DO itemp=2,NTEMP1/20
   logbeta(itemp) = logbeta(itemp-1)+(LOG(beta2)-LOG(beta1))/(REAL(NTEMP1/20,RP)-1)
ENDDO
DO itemp=NTEMP1/20+1,NTEMP1/3
   logbeta(itemp)=logbeta(itemp-1)+(LOG(beta3)-LOG(beta2))/(REAL(NTEMP1/3,RP)-REAL(NTEMP1/20,RP))
ENDDO
DO itemp=NTEMP1/3+1,NTEMP1
   logbeta(itemp)=logbeta(itemp-1)+(LOG(beta4)-LOG(beta3))/(REAL(NTEMP1,RP)-REAL(NTEMP1/3,RP))
ENDDO
beta = EXP(logbeta)
!IF(rank == src)THEN
!   OPEN(UNIT=60,FILE='beta.txt',FORM='formatted',STATUS='UNKNOWN', &
!   ACTION='WRITE',POSITION='REWIND',RECL=1024)
!   DO itemp=1,NTEMP1
!     WRITE(60,*) beta(itemp)
!   ENDDO 
!   CALL FLUSH(60)
!   CLOSE(60)
!ENDIF
!STOP

209 FORMAT(50ES16.8)
RETURN
END SUBROUTINE MAKE_BETA
!==============================================================================

SUBROUTINE AIS(obj,dat,NTEMP1,NBAL1)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,itemp,imcmc,NTEMP1,NBAL1
TYPE(objstruc) :: obj
TYPE(datastruc):: dat
REAL(KIND=RP)  :: logPLratio,logwt
INTEGER(KIND=IB),DIMENSION(NFPMX) :: idxrand
idxrand = 0
logwt = -obj%logP
!!
!! Cool chain according to beta
!!
DO itemp = 1,NTEMP1
   !!
   !! Run MCMC chain for NBAL steps
   !! (AIS does not require fully balanced chain,
   !!  so NBAL = 1 is OK)
   !!
   DO imcmc = 1,NBAL1
      CALL EXPLORE_MHAIS(obj,logPLratio,beta(itemp),dat)
      logwt = logwt - logPLratio
   ENDDO

ENDDO
!! Balance a bit more at last temp
!DO imcmc = 1,NTEMP1
!   CALL EXPLORE_MHAIS(obj,logPLratio,beta(NTEMP1),dat)
!   logwt = logwt - logPLratio
!ENDDO

obj%logwt = logwt + obj%logP + obj%logL
CALL CALCH(obj)   !! Compute z-partition for AIS object obj
CALL INIT_BD(obj) !! Init birth-death covariance matrix for each AIS particle
RETURN
END SUBROUTINE AIS
!==============================================================================

SUBROUTINE EXPLORE_MHAIS(obj,logPLratio,beta_mh,dat)
!==============================================================================
USE MPI
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,NFPnew,idel,ipar,ilay,idxz,infloop
INTEGER(KIND=IB)                            :: ncra
TYPE(objstruc)                              :: obj,objnew
TYPE(datastruc)                             :: dat
REAL(KIND=RP)                               :: logPLratio,ran_uni,ran_uni_Z,ran_uni_sd
REAL(KIND=RP)                               :: znew,beta_mh
INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP)   ,DIMENSION(obj%k)           :: ztmp

objnew = obj
!!
!! Do Metropolis-Hastings
!!
!! Perturb only c, rho, alpha, and sigma
CALL GETIDXCRA(obj,idxcra,ncra)
idxrand = 0
idxrand(1:ncra) = RANDPERM(ncra)
DO ipar = 1,ncra
   IF(i_zpert == 2)THEN
      CALL CHECKBOUNDS(objnew)
      IF(ioutside == 0)THEN
         CALL LOGLHOOD(objnew,dat)
!         logPLratio = (objnew%logL - obj%logL)*beta_mh
         logPLratio = (objnew%logP - obj%logP) + (objnew%logL - obj%logL)*beta_mh
         CALL RANDOM_NUMBER(ran_uni)
         IF(ran_uni >= EXP(logPLratio))THEN
            objnew = obj
            logPLratio = 0._RP
            ireject = ireject + 1
         ELSE
            obj = objnew
            iaccept = iaccept + 1
         ENDIF
      ELSE
         logPLratio = 0._RP
         objnew = obj
         ireject = ireject + 1
         ioutside = 0
      ENDIF 
   ELSE
      iwhich = idxcra(idxrand(ipar))
      CALL PROPOSAL_AIS(obj,objnew,dat,beta_mh,iwhich)
      IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)
      CALL CHECKBOUNDS(objnew)
      IF(ioutside == 0)THEN
         CALL LOGLHOOD(objnew,dat)
!         logPLratio = (objnew%logL - obj%logL)*beta_mh
         logPLratio = (objnew%logP - obj%logP) + (objnew%logL - obj%logL)*beta_mh
         CALL RANDOM_NUMBER(ran_uni)
         IF(ran_uni >= EXP(logPLratio))THEN
            objnew = obj
            logPLratio = 0._RP
            ireject = ireject + 1
         ELSE
            obj = objnew
            iaccept = iaccept + 1
         ENDIF
      ELSE
         logPLratio = 0._RP
         objnew = obj
         ireject = ireject + 1
         ioutside = 0
      ENDIF 
   ENDIF 
ENDDO
!!
!! Do Metropolis-Hastings on data-error standard deviations
!!
IF(ICOV == 1)THEN
   !! Perturb std devs with .25 probability
   CALL RANDOM_NUMBER(ran_uni_sd)
   IF(ran_uni_sd>=0.25_RP)THEN
      DO ipar = 1,NBAND
         CALL PROPOSAL_SD(obj,objnew,ipar)
         i_sdpert = 1
         IF(ioutside == 0)THEN
            CALL LOGLHOOD(objnew,dat)
!            logPLratio = (objnew%logL - obj%logL)*beta_mh
            logPLratio = (objnew%logP - obj%logP) + (objnew%logL - obj%logL)*beta_mh
            CALL RANDOM_NUMBER(ran_uni)
            IF(ran_uni >= EXP(logPLratio))THEN
               logPLratio = 0._RP
               objnew = obj
               ireject = ireject + 1
            ELSE
               obj = objnew
               iaccept = iaccept + 1
            ENDIF
         ELSE
            logPLratio = 0._RP
            objnew = obj
            ireject = ireject + 1
            ioutside = 0
         ENDIF
         i_sdpert = 0
      ENDDO
   ENDIF
ENDIF
RETURN
END SUBROUTINE EXPLORE_MHAIS
!=======================================================================

SUBROUTINE PROPOSAL_AIS(obj,objnew,dat,beta_mh,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,iwhich,iloop,ilay,ipar
TYPE(objstruc) :: obj,objnew
TYPE(datastruc):: dat
REAL(KIND=RP)  :: ran_uni, ran_gauss, beta_mh, factnew

objnew = obj

ilay = CEILING(REAL(iwhich,RP)/4._RP)
ipar = iwhich-(ilay-1)*NPL
IF(ilay > obj%k) ipar = ipar + 1

iloop = 0
!! Perturb only one parameter (determined in EXPLORE_MH call)
CALL RANDOM_NUMBER(ran_uni)
objnew%par = obj%par
IF(beta_mh > 1.e-2_RP)THEN
   objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(beta_mh)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
ELSE
   objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(1.e-2_RP)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
ENDIF
!!! WHY DID I SET IOUTSIDE=0 HERE????
!ioutside = 0
IF(((objnew%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - objnew%par(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_AIS
!=======================================================================

SUBROUTINE RESAMPLE_AIS(particles,particles_new,particles_old,kcount,kidx)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                   :: iipart,ipart,ik,idxpart
INTEGER(KIND=IB),DIMENSION(NLMX)   :: kcount
INTEGER(KIND=IB),DIMENSION(NPARTAIS,NLMX):: kidx
TYPE(objstruc),DIMENSION(NPART)    :: particles,particles_old,particles_new
REAL(KIND=RP),DIMENSION(NPARTAIS,NLMX):: cumsum     ! Cumulative distribution to pick random according to wt
REAL(KIND=RP):: ran_uni

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Find no. particles per k
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
kcount = 0
kidx   = 0 
DO ik = 1,NLMX
   iipart = 0
   DO ipart = 1,NPARTAIS
      IF(particles(ipart)%k == ik)THEN
         iipart = iipart + 1
         kcount(ik) = kcount(ik)+1
         kidx(iipart,ik) = ipart
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Resample new set of particles according to AIS weight
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
cumsum = 0._RP
DO ik = 1,NLMX
   IF(kcount(ik) >= 1)THEN
      cumsum(1,ik) = exp(particles(kidx(1,ik))%logwt)
      DO ipart = 2,kcount(ik)
         !! Compute cumulative logL distribution
!         cumsum(ipart) = LOGPLUS(cumsum(ipart-1),particles(ipart)%logwt)
         cumsum(ipart,ik) = cumsum(ipart-1,ik)+exp(particles(kidx(ipart,ik))%logwt)
      ENDDO ! PARTICLE LOOP
   ENDIF
ENDDO

!particles_old = particles    !! Now the old particles are the AIS ones
DO ipart = 1,NPART
   IF(particles(ipart)%par(1) == 0._RP)THEN
      particles(ipart)=particles(ipart-1)
      PRINT*,'Defect particle happened',ipart
   ENDIF
   !! Draw randomly according to weight using cumulative dist:
   CALL RANDOM_NUMBER(ran_uni)
   ran_uni = ran_uni * cumsum(kcount(particles(ipart)%k),particles(ipart)%k)
   idxpart = MINLOC(ABS(FLOOR(cumsum(1:kcount(particles(ipart)%k),particles(ipart)%k)-ran_uni)),DIM=1)
   iipart = kidx(idxpart,particles(ipart)%k)
   particles_new(ipart) = particles(iipart)
   particles_new(ipart)%logwt = LOG(1._RP/REAL(NPART,RP))
ENDDO
RETURN
END SUBROUTINE RESAMPLE_AIS
!=======================================================================

SUBROUTINE RESAMPLE(particles,particles_new,kcount,kidx)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                      :: iipart,ipart,ik,idxpart,idxk
INTEGER(KIND=IB)                      :: kcountidx
INTEGER(KIND=IB),DIMENSION(NLMX)      :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX):: kidx
TYPE(objstruc),DIMENSION(NPART)       :: particles,particles_new
REAL(KIND=RP),DIMENSION(NPARTAIS,NLMX):: cumsum     ! Cumulative distribution to pick random according to wt
REAL(KIND=RP),DIMENSION(NLMX)         :: cumsumk    ! Cumulative distribution to pick random according to wt
!REAL(KIND=RP),DIMENSION(NPARTAIS)     :: cumsum2    ! Cumulative distribution to pick random according to wt
REAL(KIND=RP)                         :: ran_uni,ran_unik
REAL(KIND=RP)                         :: LOGPLUS

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Find no. particles per k
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
kcount = 0
kidx   = 0
DO ik = 1,NLMX
   iipart = 0
   DO ipart = 1,NPARTAIS
      IF(particles(ipart)%k == ik)THEN
         iipart = iipart + 1
         kcount(ik) = kcount(ik)+1
         kidx(iipart,ik) = ipart
      ENDIF
   ENDDO ! PARTICLE LOOP
   WRITE(*,201) 'ik',ik,kidx(:,ik)
ENDDO
201 FORMAT(a,30I4)
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Resample new set of particles according to weight
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
cumsum  = 0._RP
DO ik = 1,NLMX
   IF(kcount(ik) >= 1)THEN
      cumsum(1,ik) = particles(kidx(1,ik))%logwt
      DO ipart = 2,kcount(ik)
         !! Compute cumulative logwt distribution
         cumsum(ipart,ik) = LOGPLUS(cumsum(ipart-1,ik),particles(kidx(ipart,ik))%logwt)
!         cumsum(ipart,ik) = cumsum(ipart-1,ik)+exp(particles(kidx(ipart,ik))%logwt)
      ENDDO ! PARTICLE LOOP
      !! Normalize to 1
      cumsum(:,ik) = cumsum(:,ik) - cumsum(kcount(ik),ik)
      cumsum(:,ik) = EXP(cumsum(:,ik))

   ENDIF
ENDDO
!stop
cumsumk = 0._RP
cumsumk(1) = kcount(1)
DO ik = 2,NLMX
   cumsumk(ik) = cumsumk(ik-1)+kcount(ik)
ENDDO
cumsumk = cumsumk/cumsumk(NLMX)

!DO ik = 1,NLMX
!   WRITE(78,*) cumsumk(ik)
!ENDDO
!CALL FLUSH(77)
!CALL FLUSH(78)
DO ik = 1,NLMX
   IF(kcount(ik) >= 1)THEN
      WRITE(*,*) 'ik,kcount',ik,kcount(ik)
      WRITE(*,*) cumsum(1:kcount(ik),ik)
      WRITE(*,*) particles(kidx(1:kcount(ik),ik))%logwt
   ENDIF
ENDDO

DO ipart = 1,NPART
   !! Draw randomly according to weight using cumulative dist:
   CALL RANDOM_NUMBER(ran_unik)
   CALL RANDOM_NUMBER(ran_uni)
   idxk = MINLOC(ABS(FLOOR(cumsumk-ran_unik)),DIM=1)
   kcountidx = kcount(particles(ipart)%k)
   idxpart = MINLOC(ABS(FLOOR(cumsum(1:kcountidx,idxk)-ran_uni)),DIM=1)
   iipart = kidx(idxpart,idxk)
   WRITE(*,*) 'k',particles(ipart)%k
   WRITE(*,*) 'kcount',kcount(particles(ipart)%k)
   WRITE(*,*) 'cumsum 1:kcount',cumsum(1:kcount(particles(ipart)%k),idxk)
   WRITE(*,*) 'cumsum 1:10',cumsum(1:10,idxk)
   WRITE(*,*) 'ran',ran_uni
   PRINT*,ipart,iipart,kidx(idxpart,idxk),idxpart,idxk
   particles_new(ipart) = particles(iipart)
!   WRITE(79,*) particles_new(ipart)%logwt
!   WRITE(80,*) particles(ipart)%logwt
!   WRITE(81,*) particles(ipart)%k
!   WRITE(82,*) iipart
!   WRITE(83,*) idxk
   particles_new(ipart)%logwt = LOG(1._RP/REAL(NPART,RP))
ENDDO
!DO ik = 1,NLMX
!   IF(kcount(ik) >= 1)THEN
!      WRITE(*,*) 'ik,kcount B',ik,kcount(ik)
!      WRITE(*,*) cumsum(1:kcount(ik),ik),particles(kidx(ipart,ik))%logwt
!   ENDIF
!ENDDO


!CALL FLUSH(79)
!CALL FLUSH(80)
!CALL FLUSH(81)
!CALL FLUSH(82)
!CALL FLUSH(83)

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!                                                                           !!
!! Resample new set of particles according to weight                         !!
!!                                                                           !!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!cumsum2 = 0._RP
!cumsum2(1) = particles(1)%logwt
!DO ipart = 2,NPARTAIS
!   !! Compute cumulative logwt distribution
!   cumsum2(ipart) = LOGPLUS(cumsum2(ipart-1),particles(ipart)%logwt)
!ENDDO ! PARTICLE LOOP
!!! Normalize to 1
!cumsum2 = cumsum2 - cumsum2(NPARTAIS)
!cumsum2 = EXP(cumsum2)
!DO ipart = 1,NPART
!   !! Draw randomly according to weight using cumulative dist:
!   CALL RANDOM_NUMBER(ran_uni)
!   ran_uni = ran_uni
!   idxpart = MINLOC(ABS(FLOOR(cumsum2-ran_uni)),DIM=1)
!   particles_new(ipart) = particles(idxpart)
!   particles_new(ipart)%logwt = LOG(1._RP/REAL(NPART,RP))
!ENDDO

RETURN
END SUBROUTINE RESAMPLE
!====================================================================

SUBROUTINE GASDEVJ(harvest)
!====================================================================
USE nrtype
USE nr
IMPLICIT NONE
REAL(DP), INTENT(OUT) :: harvest
REAL(DP) :: rsq,v1,v2
REAL(DP), SAVE :: g
LOGICAL, SAVE :: gaus_stored=.FALSE.
IF (gaus_stored) THEN
   harvest=g
   gaus_stored=.FALSE.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1=2.0_DP*v1-1.0_DP
      v2=2.0_DP*v2-1.0_DP
      rsq=v1**2+v2**2
      IF (rsq > 0.0_DP .AND. rsq < 1.0_DP) EXIT
   END DO
   rsq=SQRT(-2.0_DP*LOG(rsq)/rsq)
   harvest=v1*rsq
   g=v2*rsq
   gaus_stored=.TRUE.
END IF
END SUBROUTINE GASDEVJ

!====================================================================

      Function ASINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,ii,asinc
ii    = cmplx(0.,1.)
asinc = -ii*LOG(ii*z+SQRT(1.-z**2))
RETURN
END FUNCTION

!=======================================================================

Function COSC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,cosc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
cosc = CMPLX(COS(x)*COSH(y),-SIN(x)*SINH(y),RP)
RETURN
END FUNCTION

!=======================================================================

Function SINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,sinc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
sinc = CMPLX(SIN(x)*COSH(y),COS(x)*SINH(y),RP)
RETURN
END FUNCTION

!==============================================================================

FUNCTION CACOS(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CACOS
COMPLEX(KIND=RP) :: z
REAL(KIND=RP) :: zrp1,zrm1,zi,zizi,a1,a2,a,b

!CACOS = -CMPLX(0._RP,1._RP,RP)*LOG(z+CMPLX(0._RP,1._RP,RP)*SQRT(1._RP-z*z))
!!
!! This version from IDL; much faster than above
!!
zrp1 = REAL(z,RP)+1._RP
zrm1 = zrp1-2._RP
zi = AIMAG(z)
zizi = zi*zi
a1 = 0.5_RP*SQRT(zrp1*zrp1 + zizi)
a2 = 0.5_RP*SQRT(zrm1*zrm1 + zizi)
a = a1+a2
b = a1- a2
IF(zi >= 0._RP)THEN
   CACOS = ACOS(b) - CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ELSE
   CACOS = ACOS(b) + CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ENDIF

RETURN
END FUNCTION CACOS

!==============================================================================

FUNCTION ASINH(x)
!==============================================================================

USE DATA_TYPE
REAL(KIND=RP) :: ASINH
REAL(KIND=RP) :: x

ASINH = LOG(x+SQRT(x**2._RP+1))

RETURN
END FUNCTION ASINH

!==============================================================================
FUNCTION CSIN(z)
!==============================================================================
!! Complex sine (Jan's version)

USE DATA_TYPE
COMPLEX(KIND=RP) :: CSIN
COMPLEX(KIND=RP) :: z

CSIN =  (EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /CMPLX(0._RP,2._RP,RP)
RETURN
END FUNCTION CSIN


!==============================================================================
FUNCTION CTAN(z)
!==============================================================================
!! Complex TAN

USE DATA_TYPE
COMPLEX(KIND=RP) :: CTAN
COMPLEX(KIND=RP) :: z

CTAN =  -CMPLX(0._RP,1._RP,RP)*(EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /(EXP( CMPLX(0._RP,1._RP,RP)*z)+EXP(-CMPLX(0._RP,1._RP,RP)*z))
RETURN
END FUNCTION CTAN

!=======================================================================

FUNCTION LOGPLUS(x,y)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
REAL(KIND=RP)    :: logplus,x,y

IF(x > y)THEN
   LOGPLUS = x+LOG(1._RP+EXP(y-x))
ELSE
   LOGPLUS = y+LOG(1._RP+EXP(x-y))
ENDIF

RETURN
END FUNCTION LOGPLUS

!=======================================================================

SUBROUTINE PARALLEL_SEED()
!!
!!  Ensure unique random seed for each CPU
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds
REAL(KIND=RP)                                 :: ran_uni
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed1


CALL RANDOM_SEED
CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE( iseed1(iseedsize) )
IF(ISETSEED == 1)THEN
   iseed1 = (/2303055,     2435432,     5604058,     4289794,     3472290, &
      7717070,      765437,     3783525,     3087889,     4812786,     3028075, &
      3712062,     6345687,      436800,     7957708,     2047897,     1944360, &
      1222992,     7347775,     7876874,     7588112,     7590343,     1426393, &
      1753301,     7680986,     2842400,     4411488,     4804010,      497639, &
      4978920,     4675495,      754842,     7360599,     5816102/)
   CALL RANDOM_SEED(PUT=iseed1)
ELSE
   CALL RANDOM_SEED(GET=iseed1)
!   IF(rank==src)WRITE(6,*) 'Master seed:',iseed1
ENDIF

ALLOCATE( iseed(iseedsize), rseeds(iseedsize,NTHREAD), iseeds(iseedsize,NTHREAD) )
iseed = 0
rseeds = 0._RP
iseeds = 0
IF(rank == src)THEN
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1000000._RP)
ENDIF
DO i = 1,iseedsize
   CALL MPI_BCAST( iseeds(i,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
ENDDO
iseed = iseeds(:,rank+1)

!!
!! Write seeds to seed logfile:
!!
IF(rank == src)THEN
   OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
   ACTION='WRITE',POSITION='REWIND',RECL=1024)
   WRITE(50,*) 'Rank: ',rank
   WRITE(50,201) iseed
   WRITE(50,*) ''
   CLOSE(50)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
DO i = 1,NTHREAD-1
   IF(rank == i)THEN
      OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
      ACTION='WRITE',POSITION='APPEND',RECL=1024)
      WRITE(50,*) 'Rank: ',rank
      WRITE(50,201) iseed
      WRITE(50,*) ''
      CLOSE(50)
   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
ENDDO
CALL RANDOM_SEED(PUT=iseed)
DO i = 1,10000
   CALL RANDOM_NUMBER(ran_uni)
ENDDO

201   FORMAT(50I10)
END SUBROUTINE PARALLEL_SEED

!==============================================================================
FUNCTION RANDPERM(num)
!==============================================================================
USE data_type, ONLY : IB, RP
IMPLICIT NONE
INTEGER(KIND=IB), INTENT(IN) :: num
INTEGER(KIND=IB) :: number, i, j, k
INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
REAL(KIND=RP), DIMENSION(num) :: rand2
INTRINSIC RANDOM_NUMBER
CALL RANDOM_NUMBER(rand2)
DO i=1,num
   number=1
   DO j=1,num
      IF (rand2(i) > rand2(j)) THEN
           number=number+1
      END IF
   END DO
   DO k=1,i-1
      IF (rand2(i) <= rand2(k) .AND. rand2(i) >= rand2(k)) THEN
           number=number+1
      END IF
   END DO
   RANDPERM(i)=number
END DO
RETURN
END FUNCTION RANDPERM
!=======================================================================
LOGICAL FUNCTION ISNAN(a)
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) a
IF (a.NE.a) THEN
ISNAN = .TRUE.
ELSE
ISNAN = .FALSE.
END IF
RETURN
END
!=======================================================================
LOGICAL FUNCTION ISINF(a)
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) a
IF ((a*0).NE.0) THEN
ISINF = .TRUE.
ELSE
ISINF = .FALSE.
END IF
RETURN
END
!=======================================================================
! This is the end my fiend...
! EOF
