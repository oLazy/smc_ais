!==============================================================================
!
!                Sequential Reversible Jump MCMC Sampling 
!                   with AIS bridging distributions
!           for relfection coefficient inversion along track
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, September 18 2010
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: September 18 2010
!
!  Sequential particle filter based on Gilks and Berzuini 2001
!  RJMCMC based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
!  Agostinetti Malinverno 2010
!
!==============================================================================
!!$
!!$MODULE DATA_TYPE
!!$   IMPLICIT NONE
!!$   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.0D0), SP=KIND(0.0)
!!$   REAL(KIND=RP),   PARAMETER :: PI  = 3.141592653589793238462643383279502884197_RP
!!$END MODULE DATA_TYPE
!!$
!!$!==============================================================================
!!$MODULE RJMCMC_COM
!!$   USE MPI
!!$   USE DATA_TYPE
!!$   IMPLICIT NONE
!!$
!!$!
!!$! General switches
!!$!
!!$   INTEGER(KIND=IB) :: IMAP       !! WRITE REPLICA AND EXIT
!!$   INTEGER(KIND=IB) :: ICOV       !! 1 = Sample over sigma
!!$                                  !! 2 = Use sigma from ping ave
!!$                                  !! 3 = Use sigma from ping ave but scale by factor (free parameter)
!!$   INTEGER(KIND=IB) :: ENOS       !! 0 = uniform prior on k
!!$                                  !! 1 = Apply even-numbered order statistics as prior on k (avoids very thin layers)
!!$   INTEGER(KIND=IB) :: IPOIPR     !! 1 = Apply Poisson prior on k
!!$   INTEGER(KIND=IB) :: IAIS       !! 0 = Apply rjMCMC on small No. particles; 1 = Apply AIS; 
!!$   INTEGER(KIND=IB) :: IdB        !! 1=Carry out computation in dB
!!$   INTEGER(KIND=IB) :: ISETSEED   !! Fix the random seed 
!!$   INTEGER(KIND=IB) :: INODATA    !! Sample prior without data
!!$   INTEGER(KIND=IB) :: ISPHER     !! Spherical refl. coeff. 
!!$   INTEGER(KIND=IB) :: IHAM       !! Hamilton prior
!!$   INTEGER(KIND=IB) :: ISIM       !! WRITE REPLICA WITH NOISE AND EXIT
!!$   INTEGER(KIND=IB) :: subsmp     !! Subsample Sommerfeld integral plane wave part
!!$   INTEGER(KIND=IB) :: NPAVE      !! Skip pings
!!$   INTEGER(KIND=IB) :: IGA        !! Type of averaging (0=intensity, 1=Gaussian, 2=1/3 octave intensity)
!!$
!!$!!
!!$!!  AUV DATA RJMCMC trial ping
!!$!!
!!$   INTEGER(KIND=IB):: NPL         !! Number parameters per layer
!!$   INTEGER(KIND=IB):: NPING       !! Number of pings
!!$   INTEGER(KIND=IB):: NPART       !! Number of particles
!!$   INTEGER(KIND=IB):: NPARTAIS    !! Number of particles to perform AIS on
!!$   INTEGER(KIND=IB):: NANG    
!!$   INTEGER(KIND=IB):: NLMX    
!!$   INTEGER(KIND=IB):: NBAND   
!!$   INTEGER(KIND=IB):: NAVEF        !! No. freq per band
!!$   INTEGER(KIND=IB):: icheckpoint  !! Checkpointing variable (needed on some computer clusters)
!!$
!!$   INTEGER(KIND=IB):: ipingst      !! Ping # to start at (read from ./checkpoint/status.txt)
!!$   REAL(KIND=RP)   :: frbw         !! Frac. bandwidth for freq ave.
!!$   REAL(KIND=RP),ALLOCATABLE, DIMENSION(:):: bands
!!$   REAL(KIND=RP),ALLOCATABLE, DIMENSION(:):: sdtrgt
!!$   REAL(KIND=RP) :: FBW
!!$
!!$   INTEGER(KIND=IB)                  :: filebaselen, filebaselengl
!!$   CHARACTER(len=64)                 :: filebasegl, particle_init_file
!!$   CHARACTER(len=64),ALLOCATABLE,DIMENSION(:):: pingfilebase, particlefile
!!$
!!$!!
!!$!!  Prior variables and good seeding model
!!$!!
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pk       ! Poisson prior on k
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlim
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlim
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpert
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsd
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsc
!!$   INTEGER(KIND=IB)            :: kmin               ! Min number of layers
!!$   INTEGER(KIND=IB)            :: kmax               ! Max number of layers
!!$   REAL(KIND=RP)               :: lambda             ! Lambda parameter for Poisson prior on k
!!$   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness
!!$   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: fr,fstep   ! Total frequency array
!!$   REAL(KIND=RP)               :: z_t,cw,rw,hmx
!!$   REAL(KIND=RP)               :: logLtrgt
!!$   INTEGER(KIND=IB)            :: kmintmp,kmaxtmp
!!$   CHARACTER(LEN=64) :: logfile
!!$   CHARACTER(LEN=64) :: seedfile
!!$   CHARACTER(len=64) :: mapfile
!!$   CHARACTER(len=64) :: repfile
!!$   INTEGER(kind=IB) :: ulog
!!$!!
!!$!!  Standard deviation prior variables:
!!$!!
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsd
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlimsd
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpertsd
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsd
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdsc
!!$!!
!!$!!  Sampling specific parameters
!!$!!
!!$   INTEGER(KIND=IB)           :: NFPMX                             !! = (NLMX * NPL) + (NPL-1)
!!$   INTEGER(KIND=IB)           :: len_snd,len_rcv 
!!$   INTEGER(KIND=IB)           :: ioutside = 0
!!$   INTEGER(KIND=IB)           :: ireject = 0, iaccept = 0
!!$   INTEGER(KIND=IB)           :: ireject_bd = 0, iaccept_bd = 0
!!$   INTEGER(KIND=IB)           :: i_bd                              !! Birth-Death track (0=MCMC, 1=birth, 2=death)
!!$   INTEGER(KIND=IB)           :: i_zpert                           !! Z perturb track (0=nothing, 1=birth-death, 2=perturb 1 z)
!!$   INTEGER(KIND=IB)           :: ipartfail 
!!$   REAL(KIND=RP)              :: acc_ratio
!!$
!!$!!
!!$!!  Convergence parameters
!!$!!
!!$   INTEGER(KIND=IB)       :: iskipais = 0                    !! Skip AIS if initial rjMCMC and resampling successful
!!$
!!$!!
!!$!! RJMCMC parameters
!!$!!
!!$   INTEGER(KIND=IB):: MCMC_BI     !! # balance steps in traget region burn-in (exits if target reached)
!!$   INTEGER(KIND=IB):: MCMC_STEPS  !! # balance steps in traget region (always does these)
!!$   INTEGER(KIND=IB):: NTEMP       !! # temperature steps for initial AIS (if fail, increases by factor 3)
!!$   INTEGER(KIND=IB):: NTEMP1        
!!$   INTEGER(KIND=IB):: NBAL        !! # steps to balance at each AIS temp (if fail, increases by factor 1.3)
!!$   INTEGER(KIND=IB):: NBAL1         
!!$   REAL(KIND=IB)   :: npercnt     !! % particles that must be in logL target region
!!$   INTEGER(KIND=IB):: NCOOLTRAJ   !! # cooling trajectories to try before giving up
!!$   INTEGER(KIND=IB):: NAP         !! Misc parameters in sample (for bookeeping)
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: beta             !! Inverse Temprerature
!!$   REAL(KIND=RP)                         :: beta1gl          !! Global inverse T to define annealing schedule (start)
!!$   REAL(KIND=RP)                         :: beta4gl          !! Global inverse T to define annealing schedule (end for 3 exp legs)
!!$
!!$!!
!!$!!  Structures for objects and data 
!!$!!
!!$   TYPE :: objstruc
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     !! Forward parameters
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdpar   !! Std dev parameters
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: z       !! Depth partition
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: h       !! Depth partition
!!$      INTEGER(KIND=IB)                      :: k       !! Layer dimension
!!$      INTEGER(KIND=IB)                      :: NFP     !! Number forward parameters
!!$      INTEGER(KIND=IB)                      :: nfail   !! Counter to keep track of how often particle failed
!!$      REAL(KIND=RP)                         :: logL    !! log likelihood
!!$      REAL(KIND=RP)                         :: logP    !! log likelihood
!!$      REAL(KIND=RP)                         :: logwt   !! log wt (likelihood ratio when moving from one ping to next)
!!$      REAL(KIND=RP)                         :: logPr   !! log Prior probability ratio
!!$      REAL(KIND=RP)                         :: lognorm !! Data covariance matrices
!!$   END TYPE objstruc
!!$
!!$   TYPE :: datastruc
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: Robs    !! Observed data
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: angobs  !! Observed angles
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: Rrep    !! Replica data for trial models
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: res     !! Replica data for trial models
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: sigma   !! Index for bad points/data gaps
!!$                                                           !! meaningful for simulations)
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: lognorm !! Data lognorm
!!$      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: logdet  !! Data log-determinant
!!$      INTEGER(KIND=IB),ALLOCATABLE, DIMENSION(:):: NDPF    !! No. data per freq
!!$      INTEGER(KIND=IB)                          :: NANG    !! No. data/angles (struc copy)
!!$      INTEGER(KIND=IB)                          :: NBAND   !! No. freq bands (struc copy)
!!$   END TYPE datastruc
!!$
!!$   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: sdev        !! Standard devs
!!$   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: icount
!!$
!!$!!
!!$!!  Global variables for spherical reflection coeff computation
!!$!!
!!$   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE    :: N,Nr
!!$   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE       :: drTh
!!$   COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE    :: diTh
!!$   COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE  :: rTh,RPs
!!$   COMPLEX(KIND=RP),DIMENSION(:,:,:),ALLOCATABLE:: btR
!!$
!!$!!
!!$!!  MPI global variables
!!$!!
!!$  INTEGER(KIND=IB)            :: rank,NTHREAD,ierr
!!$  INTEGER(KIND=IB), PARAMETER :: src = 0_IB
!!$  INTEGER                     :: to,from,tag,COMM
!!$  INTEGER                     :: status(MPI_STATUS_SIZE)
!!$  INTEGER(KIND=IB)            :: isize1,isize2,isize3
!!$
!!$  REAL(KIND=RP)               :: tcmp,tstartcmp,tendcmp    !! Timing
!!$  REAL(KIND=RP)               :: tstartsnd, tendsnd        !! Communication time
!!$  REAL(KIND=RP)               :: tstart, tend              !! Overall time 
!!$  REAL(KIND=RP)               :: tiping1, tiping2          !! Overall time 
!!$  REAL(KIND=RP)               :: tstart2, tend2            !! Time for one forward model computation
!!$
!!$  INTERFACE
!!$    FUNCTION RANDPERM(num)
!!$       USE data_type, ONLY : IB
!!$       IMPLICIT NONE
!!$       INTEGER(KIND=IB), INTENT(IN) :: num
!!$       INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
!!$    END FUNCTION RANDPERM
!!$  END INTERFACE
!!$  CONTAINS
!!$  !!==============================================================================
!!$  INTEGER FUNCTION NEWUNIT(unit)
!!$  !!==============================================================================
!!$    integer, intent(out), optional :: unit
!!$    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
!!$    logical :: opened
!!$    integer :: lun
!!$    newunit=-1
!!$    do lun=LUN_MIN,LUN_MAX
!!$      inquire(unit=lun,opened=opened)
!!$      if (.not. opened) then
!!$        newunit=lun
!!$        exit
!!$      end if
!!$    end do
!!$    if (present(unit)) unit=newunit
!!$  END FUNCTION NEWUNIT
!!$  !!=======================================================================
!!$
!!$  SUBROUTINE ALLOC_STRUC(dat,particles,particles_new,obj,objnew)
!!$  !!
!!$  !! Allocates all derived types. 
!!$  !!
!!$  !!=======================================================================
!!$  USE MPI
!!$  !USE RJMCMC_COM
!!$  IMPLICIT NONE
!!$  TYPE(objstruc),ALLOCATABLE,DIMENSION(:):: particles,particles_new
!!$  TYPE(objstruc):: obj,objnew
!!$  TYPE(datastruc),ALLOCATABLE,DIMENSION(:):: dat
!!$  INTEGER(KIND=IB):: ip
!!$
!!$  ALLOCATE(dat(NPING))
!!$  DO ip = 1,NPING
!!$    ALLOCATE(dat(ip)%Robs(NBAND,NANG))
!!$    ALLOCATE(dat(ip)%angobs(NANG))
!!$    ALLOCATE(dat(ip)%Rrep(NBAND,NANG))
!!$    ALLOCATE(dat(ip)%res(NBAND,NANG))
!!$    ALLOCATE(dat(ip)%sigma(NBAND,NANG))
!!$    ALLOCATE(dat(ip)%lognorm(NBAND))
!!$    ALLOCATE(dat(ip)%logdet(NBAND))
!!$    ALLOCATE(dat(ip)%NDPF(NBAND))
!!$    dat(ip)%NANG  = NANG  ! # data/angles (struc copy)
!!$    dat(ip)%NBAND = NBAND ! # freq bands (struc copy)
!!$  ENDDO
!!$
!!$  ALLOCATE(particles(NPART),particles_new(NPART))
!!$  DO ip = 1,NPART
!!$    ALLOCATE(particles(ip)%par((NLMX*NPL)+(NPL-1)))
!!$    ALLOCATE(particles(ip)%sdpar(NBAND))
!!$    ALLOCATE(particles(ip)%z(NLMX))
!!$    ALLOCATE(particles(ip)%h(NLMX))
!!$    ALLOCATE(particles_new(ip)%par((NLMX*NPL)+(NPL-1)))
!!$    ALLOCATE(particles_new(ip)%sdpar(NBAND))
!!$    ALLOCATE(particles_new(ip)%z(NLMX))
!!$    ALLOCATE(particles_new(ip)%h(NLMX))
!!$  ENDDO
!!$  ALLOCATE(obj%par((NLMX*NPL)+(NPL-1)))
!!$  ALLOCATE(obj%sdpar(NBAND))
!!$  ALLOCATE(obj%z(NLMX))
!!$  ALLOCATE(obj%h(NLMX))
!!$  ALLOCATE(objnew%par((NLMX*NPL)+(NPL-1)))
!!$  ALLOCATE(objnew%sdpar(NBAND))
!!$  ALLOCATE(objnew%z(NLMX))
!!$  ALLOCATE(objnew%h(NLMX))
!!$  RETURN
!!$  END SUBROUTINE ALLOC_STRUC
!!$
!!$END MODULE RJMCMC_COM
!!$!=======================================================================

PROGRAM  SRJMCMC_PLANE

!=======================================================================
USE MPI
USE RJMCMC_COM
USE M_VALMED
USE UTILS
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: iping,iiping,ipart,ipartsnd,ipartrcv,iloop,ithread,idest,imcmc
INTEGER(KIND=IB)  :: idone,icountmcmc,ithin,ifail,ik,iipart,itrgtcnt
INTEGER(KIND=IB)  :: isource,nfailtmp,irealloc
INTEGER(KIND=IB)  :: ifr,ustat,upartcl,itry
INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:)  :: idxfail,idxpart
INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:)  :: kcount
INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: kidx

TYPE(objstruc)                             :: obj           ! Objects in RJMH chain
TYPE(objstruc)                             :: objnew        ! Objects in RJMH chain
TYPE(objstruc),ALLOCATABLE,DIMENSION(:)    :: particles     ! All particles (current ping)
TYPE(objstruc),ALLOCATABLE,DIMENSION(:)    :: particles_new ! All particles (current ping)
TYPE(datastruc),ALLOCATABLE,DIMENSION(:)   :: dat           ! Data

REAL(KIND=RP)                              :: ran_uni    ! Uniform random number
REAL(KIND=RP)                              :: beta1,beta4! Inverse Ts to pass to MAKE_BETA
REAL(KIND=RP)                              :: logLG      ! Global best log likelihood
REAL(KIND=RP)                              :: LOGPLUS    ! Function: Addition carried out in log space

REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sample_p01 ! PPD for ping 01 (to get initial particles

INTERFACE
   FUNCTION LOGFACTORIAL(n)
     USE DATA_TYPE
     REAL(KIND=RP) :: LOGFACTORIAL
     REAL(KIND=RP),INTENT(IN):: n
   END FUNCTION LOGFACTORIAL
END INTERFACE
!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds
!!---------------------------------------------------------------------!
!!     EXECUTABLE STATEMENTS START HERE                                !
!!---------------------------------------------------------------------!
!!
CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )

CALL READPARFILE()
ALLOCATE(idxfail(NPART),idxpart(NPART))
ALLOCATE(kcount(NLMX))
ALLOCATE(kidx(NPART,NLMX))
CALL ALLOC_STRUC(dat,particles,particles_new,obj,objnew)

!!
!! Do checkpointing if applicable
!!

CALL GETCHECKPOINT(particles,particle_init_file)
obj = particles(1)

!!
!! All threads READ ALL DATA for all pings
!!
!iiping = ipingst
iiping = 1
DO iping = 1,NPING
   WRITE(pingfilebase(iping),'(a,i3.3,a)')'p',iiping,'_1000_2400'
   filebaselen = 14
!   WRITE(pingfilebase(iping),'(a,i3.3,a)')'p',iiping,'_pave_3_1000_2400'
!   filebaselen = 21
   !! Read data for new ping
   CALL READ_DATA(dat(iping),pingfilebase(iping),particlefile(iping))
   iiping = iiping+NPAVE
ENDDO

IF(rank == src)WRITE(ulog,208),'Done reading data from files ',pingfilebase(1),' to ',pingfilebase(NPING),'.'

!!
!!  Prior bounds
!!
minlim(1) = hmin
maxlim(1) = hmx

kmintmp = kmin
kmaxtmp = kmax
maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc

CALL WRITE_LOG()      !! Write sampling parameters to log
CALL PARALLEL_SEED()  !! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)

ALLOCATE( icount(NTHREAD) )
icount = 0

tstart = MPI_WTIME()

CALL INIT_FREQ(dat(1))  !! Compute freq for freq-averaging, Initialize Sommerfeld integral for Spherical R

IF(IMAP == 1)THEN
  IF(rank == src)WRITE(ulog,*) 'IMAP activated, exiting after computing replica for MAP.'
  CALL COMPUTE_MAP(obj,dat(1))
  CALL FLUSH(ulog)
  STOP
ELSE
  IF(rank == src)WRITE(ulog,*) 'Starting model:'
  IF(rank == src)CALL PRINTPAR(obj)
  tstart2 = MPI_WTIME()
  IF(INODATA == 0)THEN
    CALL LOGLHOOD(obj,dat(1))
  ELSE
    CALL LOGLHOOD2(obj)
  ENDIF
  icount(rank+1) = icount(rank+1) + 1
  tend2 = MPI_WTIME()
  IF(rank == src)WRITE(ulog,*) 'logL = ',obj%logL
  IF(rank == src)WRITE(ulog,*) 'time = ',tend2-tstart2
ENDIF

!------------------------------------------------------------------------
!
!       ************ Sequential RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
logLtrgt = 0._RP
ifail = 0

iiping = ipingst+NPAVE
DO iping = 2,NPING
  IF(rank==src)THEN
    WRITE(ulog,*)   ''
    WRITE(ulog,*)   ''
    WRITE(ulog,207) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(ulog,207) '~~                                                                                             ~~'
    WRITE(ulog,206) '~~                           Working on ping ',iiping,'                                             ~~'
    WRITE(ulog,207) '~~                                                                                             ~~'
    WRITE(ulog,207) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    CALL FLUSH(ulog)
    tiping1 = MPI_WTIME()
  ENDIF

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !! Initialize filter: Compute logL for dat(iiping) & resample (all on master)
  !!                    Resample from NPART->NPARTAIS
  !! Eric: Should probably use OpenMP to parallelize the logL comp loop...
  !! (so that it is at least parallel on the node the master is on)
  !!
  CALL FILTER_INIT(particles(1:NPART),particles_new(1:NPARTAIS),obj,dat(iiping),&
                   iiping,beta1,beta4)
  particles(1:NPARTAIS) = particles_new(1:NPARTAIS)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

  !! rjMCMC for subset of particles (NPARTAIS)
  CALL FILTER_rjMCMC(particles(1:NPARTAIS),particles_new(1:NPARTAIS),obj,dat(iiping),iiping,MCMC_BI/4, & 
                     MCMC_STEPS/4,NPARTAIS)
  particles(1:NPARTAIS) = particles_new(1:NPARTAIS)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF(rank==src)THEN
    CALL CHECK_TRGT(particles(1:NPARTAIS),NPARTAIS)
  ENDIF
  CALL MPI_BCAST(iskipais,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )

  IF(IAIS == 1)THEN
  IF(iskipais == 0)THEN
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !! Carry out AIS filtering step for sub-set of particles
    !!
    DO itry = 1,NCOOLTRAJ
      CALL FILTER_AIS(particles(1:NPARTAIS),particles_new(1:NPARTAIS),obj,dat(iiping),&
                      iiping,NPARTAIS)
      particles = particles_new
      !! Resample step (does not change No. particles)
      IF(rank==src)THEN
        CALL RESAMPLE_AIS(particles(1:NPARTAIS),particles_new(1:NPARTAIS),NPARTAIS,NPARTAIS)
        WRITE(ulog,*)'Before resampling:'
        CALL CHECK_MINMAX(particles(1:NPARTAIS),NPARTAIS)
        WRITE(ulog,*)'After resampling:'
        CALL CHECK_MINMAX(particles_new(1:NPARTAIS),NPARTAIS)
        particles(1:NPARTAIS) = particles_new(1:NPARTAIS)
      ENDIF
      CALL CHECK_TRGT_AIS(particles_new(1:NPARTAIS),NPARTAIS,beta1,beta4,iiping,itry)
      IF(iskipais == 1)EXIT
    ENDDO
    particles = particles_new
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  ENDIF;ENDIF
  !! Resample back up from NPARTAIS->NPART (to full number of particles
  IF(rank==src)THEN
    CALL RESAMPLE_ALL_K(particles(1:NPARTAIS),particles_new(1:NPART),NPARTAIS,NPART)
    WRITE(ulog,*)'After resampling:'
    CALL CHECK_MINMAX(particles_new(1:NPART),NPART)
    particles = particles_new
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !! Carry out rjMCMC re-balancing step (full No. particles)
  !!
  CALL FILTER_rjMCMC(particles(1:NPART),particles_new(1:NPART),obj,dat(iiping), & 
                     iiping,MCMC_BI,MCMC_STEPS,NPART)
  particles = particles_new
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!  Master saves posterior particle sample
  !!
  IF(rank==src)THEN
    OPEN(NEWUNIT(upartcl),FILE=particlefile(iiping),FORM='formatted',STATUS='REPLACE', &
    ACTION='WRITE',RECL=4096)
    WRITE(ulog,201) 'Writing to file ',particlefile(iiping)
    DO ipart = 1,NPART
      WRITE(upartcl,202) particles_new(ipart)%logL, particles_new(ipart)%logPr, particles_new(ipart)%lognorm, &
                    REAL(particles_new(ipart)%k,RP),particles_new(ipart)%par,particles_new(ipart)%sdpar
    ENDDO
    CLOSE(upartcl)
    !!
    !! Write checkpointing file with ping ID
    !!
    IF(NTHREAD > 1)THEN
      icheckpoint = iiping
      WRITE(ulog,*) 'MASTER updating checkpoint/status.txt with ', icheckpoint
      OPEN(NEWUNIT(ustat),FILE='checkpoint/status.txt',FORM='formatted',STATUS='REPLACE',ACTION='WRITE')
      WRITE(ustat,211) icheckpoint
      CLOSE(ustat)
    ENDIF
    tiping2 = MPI_WTIME()
    WRITE(ulog,*) '# particles failed at this ping:',SUM(idxfail)
    WRITE(ulog,*) 'Total time for ping:            ',tiping2-tiping1
    CALL FLUSH(ulog)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  iiping = iiping+NPAVE
ENDDO !! PING LOOP

IF(rank==src)WRITE(ulog,*) 'SRJMH run finished'
IF(rank==src)CLOSE(ulog)  !! Closing log file

201 FORMAT(a,a)
202 FORMAT(500ES18.8)
204 FORMAT(I6,A4,F12.6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
206 FORMAT(A45,I5,A46)
207 FORMAT(A96)
208 FORMAT(A29,A15,A5,A15,A1)
209 FORMAT(I8,I8,I9,4F13.4,I7,F12.4)
210 FORMAT(A96)
211 FORMAT(I8)

CALL MPI_FINALIZE( ierr ) 
END PROGRAM SRJMCMC_PLANE
!=======================================================================

SUBROUTINE FILTER_INIT(particles,particles_new,obj,dat,iiping,beta1,beta4)
!=======================================================================
USE MPI
USE RJMCMC_COM
USE M_VALMED
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: iiping,ipart,ipartsnd,ipartrcv,iloop,ithread,idest,imcmc
INTEGER(KIND=IB)  :: idone,icountmcmc,ithin,ifail,ik,iipart,itrgtcnt
INTEGER(KIND=IB)  :: isource,nfailtmp,irealloc
INTEGER(KIND=IB)  :: ifr
INTEGER(KIND=IB),DIMENSION(NLMX)      :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX):: kidx

TYPE(objstruc)                         :: obj           !! Objects in RJMH chain
TYPE(objstruc),DIMENSION(NPART)        :: particles     !! All particles (current ping)
TYPE(objstruc),DIMENSION(NPARTAIS)     :: particles_new !! All particles (current ping)
TYPE(datastruc)                        :: dat           !! Data
REAL(KIND=RP)                          :: beta1,beta4   !! Inverse Ts to pass to MAKE_BETA
REAL(KIND=RP),DIMENSION(NBAND)         :: logLtrgttmp
REAL(KIND=RP),DIMENSION(NPARTAIS,NBAND):: sdtmp
REAL(KIND=RP),DIMENSION(NBAND)         :: sdmead
REAL(KIND=RP)                          :: logLmin,logLmax

!! Make cooling schedule
NTEMP1 = NTEMP
NBAL1  = NBAL
IF(allocated(beta) .EQV. .TRUE.)DEALLOCATE( beta )
ALLOCATE( beta(NTEMP1) )
beta1 = beta1gl
beta4 = beta4gl
CALL MAKE_BETA(beta1,beta4)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
irealloc = 0

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !! Randomize order of particles since NPARTAIS << NPART
  !! Not really needed?? I am resampling to NPARTAIS anyway...
  !idxpart = RANDPERM(NPART)
  !DO ipart = 1,NPART
  !  particles_old(ipart) = particles(idxpart(ipart))
  !ENDDO
  !particles = particles_old

  kmintmp = NLMX
  kmaxtmp = 0
PARTICLE_LOOP:  DO ipart = 1,NPART
    !! Compute likelihoods for new ping for all particles
    IF(INODATA == 0)THEN
      CALL LOGLHOOD(particles(ipart),dat)
    ELSE
      CALL LOGLHOOD2(particles(ipart))
    ENDIF
    !! Compute particle log-weights (given by logL)
    particles(ipart)%logwt = particles(ipart)%logL
    kmintmp = MIN(kmintmp,particles(ipart)%k)
    kmaxtmp = MAX(kmaxtmp,particles(ipart)%k)
 END DO PARTICLE_LOOP 

  IF(kmintmp > 1)    kmintmp = kmintmp - 1
  IF(kmaxtmp < NLMX) kmaxtmp = kmaxtmp + 1

  WRITE(ulog,*) 'kmin for ping',iiping,'= ',kmintmp
  WRITE(ulog,*) 'kmax for ping',iiping,'= ',kmaxtmp
  WRITE(ulog,*) 'NTEMP1 = ',NTEMP1
  WRITE(ulog,*) 'NBAL1 = ',NBAL1
  WRITE(ulog,*) 'beta1 = ',beta1
  WRITE(ulog,*) 'beta4 = ',beta4

  !!
  !! Resample new set of particles according to weight
  !! resampled array is stored in particles_new
  !!
  CALL RESAMPLE_ALL_K(particles(1:NPART),particles_new(1:NPARTAIS),NPART,NPARTAIS)

  !!
  !! Estimate target logL region
  !!
  logLmin = HUGE(1._RP)
  logLmax = -HUGE(1._RP)
  DO ipart = 1,NPARTAIS
    IF(INODATA == 0)THEN
      CALL LOGLHOOD(particles_new(ipart),dat)
    ELSE
      CALL LOGLHOOD2(particles_new(ipart))
    ENDIF
    IF(particles_new(ipart)%logL > logLmax)logLmax = particles_new(ipart)%logL
    IF(particles_new(ipart)%logL < logLmin)logLmin = particles_new(ipart)%logL
    DO ifr = 1,NBAND
      sdtmp(ipart,ifr) = particles_new(ipart)%sdpar(ifr)
    ENDDO
  ENDDO
  DO ifr = 1,NBAND
    sdmead(ifr) = VALMED(sdtmp(:,ifr))
    logLtrgttmp(ifr) = -(REAL(NANG,RP)/2._RP)*LOG(2._RP*PI) -(REAL(NANG,RP)*LOG(sdtrgt(ifr))+REAL(NANG,RP)/2._RP)
  ENDDO
  logLtrgt = SUM(logLtrgttmp)
  WRITE(ulog,213) 'median sigma(ifr) = ',sdmead
  WRITE(ulog,213) 'logLtrgt(ifr)     = ',logLtrgttmp
  WRITE(ulog,213) 'logLtrgt          = ',logLtrgt
  CALL FLUSH(ulog)

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!    NO SLAVE PART
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ENDIF

!! Broadcast target logL:
CALL MPI_BCAST( logLtrgt,1, MPI_DOUBLE_PRECISION, src, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( kmintmp,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( kmaxtmp,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )

209 FORMAT(I8,I8,I9,4F13.4,I7,F12.4)
210 FORMAT(A96)
213 FORMAT(A,10F10.4)
RETURN
END SUBROUTINE FILTER_INIT
!=======================================================================

SUBROUTINE FILTER_AIS(particles,particles_new,obj,dat,iiping,N1)
!=======================================================================
USE MPI
USE RJMCMC_COM
USE M_VALMED
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: iiping,ipart,ipartsnd,ipartrcv,iloop,ithread,idest,imcmc
INTEGER(KIND=IB)  :: idone,icountmcmc,ithin,ifail,ik,iipart,itrgtcnt
INTEGER(KIND=IB)  :: isource,nfailtmp,irealloc,aiscount
INTEGER(KIND=IB)  :: ifr,N1
INTEGER(KIND=IB),DIMENSION(N1):: idxfail
TYPE(objstruc)                :: obj           ! Objects in RJMH chain
TYPE(objstruc),DIMENSION(N1)  :: particles     ! All particles (current ping)
TYPE(objstruc),DIMENSION(N1)  :: particles_new ! All particles (current ping)
TYPE(datastruc)               :: dat           ! Data

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
   ipartfail = 0
   aiscount = 1

   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   !!
   !! Compute AIS weights
   !!
   !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   WRITE(ulog,*) 'AIS:'
   WRITE(ulog,210) '    ping,  iloop,particle,    logL new,    logL old,   logwt new,   logwt old,source, time(comp)'
   WRITE(ulog,210) '------------------------------------------------------------------------------------------------'
   CALL FLUSH(ulog)
   ipartsnd = 0
   idxfail = 0
   DO ithread = 1,MIN(NTHREAD-1,N1)
      ipartsnd = ipartsnd + 1
      idest = ipartsnd
      idone = 0
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,irealloc)
   ENDDO ! PARTICLE LOOP
   iloop = 0
   ifail = 0
   DO 
      iloop = iloop + 1
      CALL RCV_PARTICLE(obj,isource,ipartrcv,ifail)
      particles_new(ipartrcv) = obj
      idest = isource
      IF(MOD(iloop,100)==0)THEN
         WRITE(ulog,209) iiping,iloop,ipartrcv,particles_new(ipartrcv)%logL,particles(ipartrcv)%logL,&
                      particles_new(ipartrcv)%logwt,particles(ipartrcv)%logwt,isource,tcmp
         CALL FLUSH(ulog)
      ENDIF
      IF(ipartsnd < N1)THEN
         tendsnd = MPI_WTIME()
         ipartsnd = ipartsnd + 1
         idone = 0
         CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,irealloc)
      ENDIF

      IF(iloop == N1)EXIT
   ENDDO ! PARTICLE LOOP
   !!
   !! Send all slaveѕ the done if AIS was success.
   !!
   DO ithread = 1,MIN(NTHREAD-1,N1)
      idone = 1
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,ithread,idone,irealloc)
   ENDDO
   !!
   !! Resample step
   !!

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    SLAVE PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE
   !!
   !! Compute AIS weights for resampling
   !!
   DO
      ifail = 0
      idone = 0
      tstartcmp = MPI_WTIME()
      CALL SND_PARTICLE(obj,ipart,idest,idone,irealloc)
      IF(idone == 1) EXIT

      CALL AIS(obj,dat)
      tendcmp = MPI_WTIME()
      tcmp = tendcmp-tstartcmp
      !! Get updated particle back to master
      CALL RCV_PARTICLE(obj,isource,ipart,ifail)
   ENDDO

ENDIF

209 FORMAT(I8,I8,I9,4F13.4,I7,F12.4)
210 FORMAT(A96)
RETURN
END SUBROUTINE FILTER_AIS
!=======================================================================

SUBROUTINE FILTER_rjMCMC(particles,particles_new,obj,dat,iiping,MCMC1,& 
                         MCMC2,N1)
!=======================================================================
USE MPI
USE RJMCMC_COM
USE M_VALMED
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: iiping,ipart,ipartsnd,ipartrcv,iloop,ithread,idest,imcmc
INTEGER(KIND=IB)  :: idone,icountmcmc,ithin,ifail,ik,iipart,itrgtcnt
INTEGER(KIND=IB)  :: isource,nfailtmp,irealloc
INTEGER(KIND=IB)  :: ifr,MCMC1,MCMC2,N1
INTEGER(KIND=IB),DIMENSION(N1):: idxfail

TYPE(objstruc)                :: obj           ! Objects in RJMH chain
TYPE(objstruc),DIMENSION(N1)  :: particles     ! All particles (current ping)
TYPE(objstruc),DIMENSION(N1)  :: particles_new ! All particles (current ping)
TYPE(datastruc)               :: dat           ! Data

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
  DO ipart = 1,N1
    particles_new(ipart)%nfail = 0
  ENDDO ! PARTICLE LOOP
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!                                                                                !!
  !! Send new particles to slaves, then receive answers (auto-load balancing):      !!
  !!                                                                                !!
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
   
  WRITE(ulog,202) 'rjMCMC for',MCMC1,' burn-in and ',MCMC2,' retained steps:'
  WRITE(ulog,203) '    ping,  iloop,particle,idone,        logL,  k,iacc_bd, acc_ratio, time(prtcle),source, time(comp),ifail,nfail'
  WRITE(ulog,203) '----------------------------------------------------------------------------------------------------------------'
  CALL FLUSH(ulog)
  ipartsnd = 0
  idxfail = 0
  DO iloop = 1,MIN(NTHREAD-1,N1)
    ipartsnd = ipartsnd + 1
    idest = ipartsnd
    idone = 0
    CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,irealloc)
  ENDDO ! PARTICLE LOOP
  iloop = 0
  !! Auto-load balancing particle loop.
  DO 
    iloop = iloop + 1
    tstartsnd = MPI_WTIME()
    CALL RCV_PARTICLE(obj,isource,ipartrcv,ifail)
    tendsnd = MPI_WTIME()
    particles_new(ipartrcv) = obj
    idest = isource

    IF(MOD(iloop,100)==0)THEN
       WRITE(ulog,201) iiping,iloop,ipartrcv,ipartsnd,particles_new(ipartrcv)%logL,particles_new(ipartrcv)%k,&
                    iaccept_bd,acc_ratio,tendsnd-tstartsnd,isource,tcmp,ifail,idxfail(ipartrcv)
      CALL FLUSH(ulog)
    ENDIF
    IF(ipartsnd < N1)THEN
      tendsnd = MPI_WTIME()
      ipartsnd = ipartsnd + 1
      idone = 0
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,irealloc)
    ELSE
      idone = 1
      CALL SND_PARTICLE(particles(ipartsnd),ipartsnd,idest,idone,irealloc)
    ENDIF

    IF(iloop == N1)EXIT
  ENDDO ! PARTICLE LOOP

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    SLAVE PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE
  !!
  !! REBALANCE WITH rjMCMC
  !!

  !! Keep working until master says you're done
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
    CALL SND_PARTICLE(obj,ipart,idest,idone,irealloc)

    IF(idone == 1) EXIT
    imcmc = 0
    DO 
      CALL EXPLORE_MH(obj,dat,1._RP)
      imcmc = imcmc + 1
      IF(obj%logL > logLtrgt)EXIT
      IF(imcmc >= MCMC1)EXIT
    ENDDO
    DO imcmc = 1,MCMC2
      CALL EXPLORE_MH(obj,dat,1._RP)
    ENDDO
    tendcmp = MPI_WTIME()
   
    tcmp = tendcmp-tstartcmp
    !! Get updated particle back to master
    CALL RCV_PARTICLE(obj,isource,ipart,ifail)
  ENDDO
ENDIF !! MPI ENDIF

201 FORMAT(I8,I8,I9,I6,F13.4,I4,I8,F11.4,F14.4,I7,F12.4,I6,I5)
202 FORMAT(a,I4,a,I4,a)
203 FORMAT(A112)
RETURN
END SUBROUTINE FILTER_rjMCMC
!!=======================================================================

SUBROUTINE CHECK_MINMAX(particles,N1)
!! Called from master only!
!!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB):: ipart,N1,itrgtcnt
TYPE(objstruc),DIMENSION(N1):: particles     ! All particles (current ping)
REAL(KIND=RP)               :: logLmin,logLmax,logwtmin,logwtmax

!!
!! Check if at least some particles are in target region. If not, repeat with AIS three times as many interp. dists.
!!
itrgtcnt = 0
logLmin = HUGE(1._RP)
logLmax = -HUGE(1._RP)
logwtmin = HUGE(1._RP)
logwtmax = -HUGE(1._RP)
DO ipart = 1,N1
   IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
   IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
   IF(particles(ipart)%logwt > logwtmax)logwtmax = particles(ipart)%logwt
   IF(particles(ipart)%logwt < logwtmin)logwtmin = particles(ipart)%logwt
ENDDO
WRITE(ulog,*) 'logLmin  = ',logLmin, 'logLmax  = ',logLmax
WRITE(ulog,*) 'logwtmin = ',logwtmin,'logwtmax = ',logwtmax
CALL FLUSH(ulog)
RETURN
END SUBROUTINE CHECK_MINMAX
!!=======================================================================

SUBROUTINE CHECK_TRGT(particles,N1)
!! Called from master only!
!!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB):: ipart,N1,itrgtcnt
TYPE(objstruc),DIMENSION(N1):: particles     ! All particles (current ping)
REAL(KIND=RP)               :: logLmin,logLmax

!!
!! Check if at least some particles are in target region. If not, repeat with AIS three times as many interp. dists.
!!
itrgtcnt = 0
logLmin = HUGE(1._RP)
logLmax = -HUGE(1._RP)
DO ipart = 1,N1
   IF(particles(ipart)%logL > logLtrgt)itrgtcnt = itrgtcnt + 1
   IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
   IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
ENDDO
IF(itrgtcnt < NINT(npercnt*REAL(N1,RP)/100._RP))THEN
  iskipais = 0
  WRITE(ulog,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
  WRITE(ulog,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
  WRITE(ulog,*) itrgtcnt,'particles are in target region. Proceeding with AIS.'
ELSE
  iskipais = 1
  WRITE(ulog,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
  WRITE(ulog,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
  WRITE(ulog,*) itrgtcnt,'particles are in target region. Skipping AIS.'
ENDIF

RETURN
END SUBROUTINE CHECK_TRGT
!!=======================================================================

SUBROUTINE CHECK_TRGT_AIS(particles,N1,beta1,beta4,iiping,itry)
!! Called from all threads!
!!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB)            :: ipart,N1,itrgtcnt,iiping,itry
TYPE(objstruc),DIMENSION(N1):: particles     ! All particles (current ping)
REAL(KIND=RP)               :: beta1,beta4! Inverse Ts to pass to MAKE_BETA
REAL(KIND=RP)               :: logLmin,logLmax

IF(rank==src)THEN
  !!
  !! Check if at least some particles are in target region. If not, repeat with AIS three times as many interp. dists.
  !!
  itrgtcnt = 0
  logLmin = HUGE(1._RP)
  logLmax = -HUGE(1._RP)
  DO ipart = 1,N1
     IF(particles(ipart)%logL > logLtrgt)itrgtcnt = itrgtcnt + 1
     IF(particles(ipart)%logL > logLmax)logLmax = particles(ipart)%logL
     IF(particles(ipart)%logL < logLmin)logLmin = particles(ipart)%logL
  ENDDO
  iskipais = 1
  IF(itrgtcnt < NINT(npercnt*REAL(N1,RP)/100._RP))iskipais = 0
  WRITE(ulog,*) 'logL target:',logLtrgt,'# in target:',itrgtcnt
  WRITE(ulog,*) 'logLmin = ',logLmin,'logLmax = ',logLmax
  WRITE(ulog,*) itrgtcnt,'particles are in target region.'
  CALL FLUSH(ulog)
ENDIF
IF(itry==NCOOLTRAJ)iskipais = 1
CALL MPI_BCAST(iskipais,1, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
IF(iskipais==0)THEN
  beta1 = beta1/1.5_RP
  beta4 = beta4gl
  NBAL1 = NBAL1 + 1
  DEALLOCATE(beta)
  NTEMP1 = 3_IB*NTEMP1
  ALLOCATE(beta(NTEMP1))
  CALL MAKE_BETA(beta1,beta4)
  IF(rank==src)THEN
    WRITE(ulog,*) 'Repeating AIS for ping',iiping
    WRITE(ulog,*) 'NTEMP1 = ',NTEMP1
    WRITE(ulog,*) 'NBAL1 = ',NBAL1
    WRITE(ulog,*) 'beta1 = ',beta1
    WRITE(ulog,*) 'beta4 = ',beta4
  ENDIF
ELSE
  IF(rank==src)WRITE(ulog,*) 'Proceeding with rjMCMC.'
ENDIF

RETURN
END SUBROUTINE CHECK_TRGT_AIS
!!=======================================================================

SUBROUTINE LOGLHOOD(obj,dat)
!!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB) :: ifreq,ifr,iang,ncra,ipar,iidx
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP),DIMENSION(obj%NFP)            :: m_inr
REAL(KIND=RP),DIMENSION(NAVEF,obj%k+1)      :: vp,alf1dB
REAL(KIND=RP),DIMENSION(NAVEF,NANG)         :: Rpltry
REAL(KIND=RP),DIMENSION(NBAND)              :: Etmp
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+1)          :: idxc,idxa
REAL(KIND=RP)                               :: fref,cref,alfrdB,y
REAL(KIND=RP), DIMENSION(NFPMX)             :: mtmp
INTEGER(KIND=IB),DIMENSION(obj%k)           :: idxh

!!
!!  Compute plane wave refl. coeff. (band averaged)
!!
mtmp = 0._RP
mtmp(1:obj%NFP) = obj%par(1:obj%NFP)
idxh  = 0
idxh  = (/1:obj%k*NPL:NPL/)
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
                           Rpltry,NAVEF,obj%NFP,ifreq)
      ELSE
         CALL SPH_REF_NLAY(dat%angobs,mtmp(1:obj%NFP),fr(ifreq),&
                           Rpltry,NAVEF,obj%NFP,ifreq)
      ENDIF
   ENDIF

   IF(NAVEF > 1)THEN
      IF(IGA == 1)THEN
         DO iang = 1,NANG
            !!
            !! GAUSSIAN AVERAGE
            !!
            dat%Rrep(ifreq,iang) = SUM(Rpltry(:,iang)*EXP(-(fr((ifreq-1)*NAVEF+1:ifreq*NAVEF)-bands(ifreq))**2/ &
                                   (frbw*bands(ifreq))**2)*fstep(ifreq),1)/ &
                                   SUM(EXP(-(fr-bands(ifreq))**2._RP/ &
                                  (frbw*bands(ifreq))**2)*fstep(ifreq),1);
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

!!
!!  Compute log likelihood
!!
dat%res = (dat%Robs-dat%Rrep)!*dat%Rex
IF(ICOV == 1)THEN
   !!
   !! Sample over sigma (one for all freqs)
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = -(REAL(NANG,RP)/2._RP)*LOG((2._RP*PI)) &
          -(SUM(dat%res(ifreq,:)**2)/(2._RP*obj%sdpar(ifreq)**2)+REAL(NANG,RP)*LOG(obj%sdpar(ifreq)))
   ENDDO
ELSEIF(ICOV == 2)THEN
   !!
   !! Keep sigma fixed at ping avergaed estimate
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = dat%lognorm(ifreq) &
          -(SUM(dat%res(ifreq,:)**2/(2._RP*dat%sigma(ifreq,:)**2)))
   ENDDO
ELSEIF(ICOV == 3)THEN
   !!
   !! Sample over scaling (one for all freqs) for ping averaged sigma estimate
   !!
   DO ifreq = 1,NBAND
      Etmp(ifreq) = -(REAL(NANG,RP)/2._RP)*LOG((2._RP*PI)) &
          -(SUM(dat%res(ifreq,:)**2/(2._RP*(obj%sdpar(1)*dat%sigma(ifreq,:))**2))&
          +REAL(NANG,RP)*LOG(obj%sdpar(ifreq))+SUM(LOG(dat%sigma(ifreq,:))))
   ENDDO
ENDIF
obj%logL = SUM(Etmp)
obj%logP = LOG(PRODUCT(1._RP/maxpert))

RETURN
END SUBROUTINE LOGLHOOD
!!==============================================================================

SUBROUTINE LOGLHOOD2(obj)
!!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj

obj%logL = 1._RP

END SUBROUTINE LOGLHOOD2
!!==============================================================================

SUBROUTINE EXPLORE_MH(obj,dat,beta_mh)
!!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,idel,ipar,ipar2,ilay,idxz,i_comp_pred
INTEGER(KIND=IB)                            :: ncra
TYPE(objstruc)                              :: obj,objnew1
TYPE (datastruc)                            :: dat
REAL(KIND=RP)                               :: logPLratio,logPratio,logy,logq1_1,logq1_2
REAL(KIND=RP)                               :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                               :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp
REAL(KIND=RP)                               :: logarp
INTEGER(KIND=IB)                            :: arptype

objnew1 = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)
i_bd = 0

!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(kmin /= kmax)THEN
   !! Perturbing k:
   CALL RANDOM_NUMBER(ran_unik)
   i_bd = 0
   IF(obj%k == kmax)THEN  !! If k == kmax, no birth allowed
      IF(ran_unik<=0.3333_RP) i_bd = 2
   ELSEIF(obj%k == kmin)THEN  !! If k == kmin, no death allowed
      IF(ran_unik<=0.3333_RP) i_bd = 1
   ELSE
      IF(ran_unik <= 0.3333_RP) i_bd = 1
      IF(ran_unik > 0.6666_RP) i_bd = 2
   ENDIF
   IF(i_bd == 1)THEN
      CALL BIRTH(obj,objnew1)
   ELSEIF(i_bd == 2)THEN
      CALL DEATH(obj,objnew1)
   ENDIF
ENDIF

!!
!! Do Metropolis-Hastings
!!
IF(obj%k /= objnew1%k)THEN
   !!
   !! If k changed, check BD acceptance
   !!
   !obj%ipropose_bd = obj%ipropose_bd + 1
   CALL CHECKBOUNDS(objnew1)
   IF(ioutside == 0)THEN
      IF(INODATA == 0)THEN
        i_comp_pred = 1
        CALL LOGLHOOD(objnew1,dat)
      ELSE
        CALL LOGLHOOD2(objnew1)
      ENDIF
      logPratio  = objnew1%logPr
      logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh
      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
         objnew1 = obj
      ELSE
         obj = objnew1
         !obj%iaccept_bd = obj%iaccept_bd + 1
      ENDIF
   ELSE
      objnew1 = obj
      ioutside = 0
   ENDIF
ELSE  ! k-change if
!IF(1==2)THEN
   !!
   !! If k did not change, carry out MH sweep
   !!
   idxrand = 0
   idxrand(1:obj%NFP) = RANDPERM(obj%NFP)
   !!
   !! Do Metropolis-Hastings update on c, rho, alpha
   !!
   DO ipar = 1,obj%NFP
     iwhich = idxrand(ipar)
     ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
     ipar2  = iwhich-(ilay-1)*NPL
     IF(ilay > obj%k) ipar2 = ipar2 + 1

     !ipropose(ipar2) = ipropose(ipar2) + 1
     CALL PROPOSAL(obj,objnew1,iwhich,1._RP)
     CALL CHECKBOUNDS(objnew1)
     IF(ioutside == 0)THEN
       IF(INODATA == 0)THEN
         i_comp_pred = 1
         CALL LOGLHOOD(objnew1,dat)
       ELSE
         CALL LOGLHOOD2(objnew1)
       ENDIF
       logPratio  = objnew1%logPr
       logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh
       !logPLratio = (objnew1%logL - obj%logL)*beta_mh
       CALL RANDOM_NUMBER(ran_uni)
       IF(ran_uni >= EXP(logPLratio))THEN
         objnew1 = obj
       ELSE
         obj = objnew1
         !iaccept(ipar2) = iaccept(ipar2) + 1
       ENDIF
     ELSE !! outside before delayed rejection
       objnew1 = obj
       ioutside = 0
     ENDIF
   ENDDO
!   ipropose_buf(:,ibuf,ic) = ipropose
!   iaccept_buf(:,ibuf,ic)  = iaccept
!   IF(ibuf == NBUF) ibuf = 0
!   ibuf = ibuf + 1

   !!
   !! Do Metropolis-Hastings on data-error standard deviations
   !!
   IF(ICOV == 1)THEN
     DO ipar = 1,NBAND
        !ipropose_sd = ipropose_sd + 1
        CALL PROPOSAL_SD(obj,objnew1,ipar)
        IF(ioutside == 0)THEN
           IF(INODATA == 0)THEN
             i_comp_pred = 0
             CALL LOGLHOOD(objnew1,dat)
           ELSE
             CALL LOGLHOOD2(objnew1)
           ENDIF
           logPLratio = (objnew1%logL - obj%logL)*beta_mh
           CALL RANDOM_NUMBER(ran_uni)
           IF(ran_uni >= EXP(logPLratio))THEN
              objnew1 = obj
           ELSE
              obj = objnew1
              !iaccept_sd = iaccept_sd + 1
           ENDIF
        ELSE
           objnew1 = obj
           ioutside = 0
        ENDIF
     ENDDO
   ENDIF ! ICOV if
ENDIF ! k-change if
END SUBROUTINE EXPLORE_MH
!!=======================================================================

SUBROUTINE CHECKBOUNDS(obj)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: ip,iwhich,ipar,ilay,ncra
INTEGER(KIND=IB)                            :: ih,OS,ihamrej
TYPE(objstruc)                              :: obj
REAL(KIND=RP)                               :: vspmin,vspmax,zhere
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp

IF(obj%par((obj%k-1)*NPL+1) > maxlim(1))ioutside = 1
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
   IF(((obj%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - obj%par(iwhich)) < 0._RP))THEN
     ioutside = 1
     IF(rank == src)WRITE(ulog,*) 'ipar=',ipar,iwhich,minlim(ipar),maxlim(ipar),obj%par(iwhich)
   ENDIF
ENDDO

ihamrej = 0
IF(IHAM == 1)CALL COUPLE_HAM(obj,ihamrej)
IF(ihamrej == 1) ioutside = 1

RETURN
END SUBROUTINE CHECKBOUNDS
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,iwhich,factor)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,ipar,ipar2,ilay,ivel,ihamrej
TYPE(objstruc) :: obj,objnew,objtmp
REAL(KIND=RP)  :: ran_uni, ran_nor, factor
REAL(KIND=RP)  :: zp,zjm1,zj,zjp1

objnew = obj
ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
ipar  = iwhich-(ilay-1)*NPL
IF(ilay > obj%k) ipar = ipar + 1
!! Initialize prior as zero (if ENOS is applied, need to compute). 
objnew%logPr = 0._RP

IF(ipar /= 1)THEN
  !! CAUCHY proposal
  CALL RANDOM_NUMBER(ran_uni)
  objnew%par(iwhich) = obj%par(iwhich) + fact/factor*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
ELSE
  !!
  !! Compute prior for even-numbered order statistics (Green 1995) where 
  !! znew is perturbed interface, zj (j) is original interface, 
  !! zjp1 is interface below (j+1) and zjm1 is interface above.
  !!
  IF(ENOS == 0)THEN
    !! CAUCHY proposal
    CALL RANDOM_NUMBER(ran_uni)
    objnew%par(iwhich) = obj%par(iwhich) + fact/factor*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
    objnew%logPr = 0._RP
  ELSE
    CALL RANDOM_NUMBER(ran_uni)
    zj = obj%par(iwhich)
    IF(ilay == 1)THEN
      zjm1 = 0._RP
    ELSE
      zjm1 = obj%par(iwhich-NPL)
    ENDIF
    IF(ilay == obj%k)THEN
      zjp1 = hmx
    ELSE
      zjp1 = obj%par(iwhich+NPL)
    ENDIF
    !! sample uniform:
    zp = zjm1+ran_uni*(zjp1-zjm1)
    objnew%par(iwhich) = zp
    !! Apply even-numbered order statistics in prior:
    objnew%logPr = LOG(zjp1-zp)+LOG(zp-zjm1)-LOG(zjp1-zj)-LOG(zj-zjm1)
  ENDIF
ENDIF
CALL CALCH(objnew)

RETURN
END SUBROUTINE PROPOSAL
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
REAL(KIND=RP)  :: zp,zjm1,zj,zjp1

objnew = obj
!! Initialize prior as zero (if ENOS is applied, need to compute). 
objnew%logPr = 0._RP

ilay = CEILING(REAL(iwhich,RP)/4._RP)
ipar = iwhich-(ilay-1)*NPL
IF(ilay > obj%k) ipar = ipar + 1

iloop = 0
!! Perturb only one parameter (determined in EXPLORE_MH call)
CALL RANDOM_NUMBER(ran_uni)


IF(ipar /= 1)THEN
  !! CAUCHY proposal
  CALL RANDOM_NUMBER(ran_uni)
  IF(beta_mh > 1.e-2_RP)THEN
    objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(beta_mh)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
  ELSE
    objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(1.e-2_RP)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
  ENDIF
  IF(((objnew%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - objnew%par(iwhich)) < 0._RP))ioutside = 1
ELSE
  !!
  !! Compute prior for even-numbered order statistics (Green 1995) where 
  !! znew is perturbed interface, zj (j) is original interface, 
  !! zjp1 is interface below (j+1) and zjm1 is interface above.
  !!
  IF(ENOS == 0)THEN
    !! CAUCHY proposal
    CALL RANDOM_NUMBER(ran_uni)
    IF(beta_mh > 1.e-2_RP)THEN
      objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(beta_mh)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
    ELSE
      objnew%par(iwhich) = obj%par(iwhich)+1._RP/SQRT(1.e-2_RP)*pertsd(ipar)*TAN(PI*(ran_uni-0.5_RP))
    ENDIF
    IF(((objnew%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - objnew%par(iwhich)) < 0._RP))ioutside = 1
    objnew%logPr = 0._RP
  ELSE
    CALL RANDOM_NUMBER(ran_uni)
    zj = obj%par(iwhich)
    IF(ilay == 1)THEN
      zjm1 = 0._RP
    ELSE
      zjm1 = obj%par(iwhich-NPL)
    ENDIF
    IF(ilay == obj%k)THEN
      zjp1 = hmx
    ELSE
      zjp1 = obj%par(iwhich+NPL)
    ENDIF
    !! sample uniform:
    zp = zjm1+ran_uni*(zjp1-zjm1)
    objnew%par(iwhich) = zp
    !! Apply even-numbered order statistics in prior:
    objnew%logPr = LOG(zjp1-zp)+LOG(zp-zjm1)-LOG(zjp1-zj)-LOG(zj-zjm1)
  ENDIF
ENDIF
CALL CALCH(objnew)
RETURN
END SUBROUTINE PROPOSAL_AIS
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
REAL(KIND=RP)                    :: zdel,zj,zjp1
REAL(KIND=RP),DIMENSION(NPL-1)   :: cra

objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)
!! Pick random layer:
idxdeath = 0
idxdeath(1:obj%k) = RANDPERM(obj%k)
idel = idxdeath(1)

objnew%par = 0._RP  ! Take care here that there are no old "left overs" from other dimensions
cra = 0._RP
IF(idel == 1)THEN
   !! Killed interface is first interface:
   cra = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
   objnew%par(1:objnew%NFP) = obj%par(NPL+1:obj%NFP)
   objnew%par(2:NPL) = cra
   zdel = obj%par(1)
   zj   = 0._RP
   zjp1 = obj%par(1+NPL)
ELSEIF(idel == obj%k)THEN
   !! Killed interface is last interface:
   cra = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
   objnew%par(1:objnew%NFP-(NPL-1)) = (/ obj%par(1:obj%NFP-(NPL*2-1)) /)
   objnew%par(objnew%NFP-(NPL-2):objnew%NFP) = cra
   zdel = obj%par(obj%NFP-(NPL-2)-NPL)
   zj   = obj%par(obj%NFP-(NPL-2)-2*NPL)
   zjp1 = hmx-hmin
ELSE
   !! Killed interface is in normal layer stack:
   cra = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
   objnew%par(1:objnew%NFP) = (/ obj%par(1:(idel-1)*NPL),obj%par(idel*NPL+1:obj%NFP) /)
   objnew%par((idel-1)*NPL+2:(idel-1)*NPL+NPL) = cra
   zdel = obj%par(((idel-1)*NPL)+1)
   zj   = obj%par((idel-1)*NPL-(NPL-1))
   zjp1 = obj%par(idel*NPL+1)
ENDIF
CALL CALCH(objnew)

!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
!! Apply only even-numbered order statistics in prior:
  objnew%logPr = 2._RP*LOG(hmx-hmin)-LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
!! Apply only Poisson prior:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
!! Apply Poisson prior with ENOS:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+2._RP*LOG(hmx-hmin)- &
                 LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ENDIF
END SUBROUTINE DEATH
!=======================================================================

SUBROUTINE BIRTH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)              :: i, iznew, ipert, iwhich, ipar
INTEGER(KIND=IB)              :: ihamrej
TYPE(objstruc)                :: obj,objnew,objnew2
REAL(KIND=RP)                 :: ran_uni,ran_nor
REAL(KIND=RP)                 :: znew,zj,zjp1
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
   objnew%par(1:objnew%NFP) = (/ obj%par(1:NPL),obj%par /)
   objnew%par(1)            = znew
   zj                       = 0._RP
   zjp1                     = objnew%par(1+NPL)
ELSEIF(iznew > obj%k)THEN
   !! New layer is created in half-space:
   objnew%par(1:objnew%NFP) = (/ obj%par(1:obj%NFP-(NPL-1)),znew, &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP), &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP) /)
   objnew%par(objnew%NFP-((NPL*2)-2)) = znew
   zj                                 = objnew%par((objnew%NFP)-(2*NPL)-(NPL-2))
   zjp1                               = (hmx-hmin)
ELSE
   !! New layer is within layer stack:
   objnew%par(1:objnew%NFP) = (/ obj%par(1:iznew*NPL), &
                                 obj%par(((iznew-1)*NPL)+1:iznew*NPL), &
                                 obj%par(iznew*NPL+1:obj%NFP) /)
   objnew%par(((iznew-1)*NPL)+1) = znew
   zj                            = objnew%par((iznew-1)*NPL-(NPL-1))
   zjp1                          = objnew%par(iznew*NPL+1)
ENDIF
CALL CALCH(objnew)

!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
  objnew%logPr = LOG(2._RP*obj%k+2._RP)+LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+LOG(2._RP*obj%k+2._RP)+ &
                 LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ENDIF

ipert = iznew
objnew2 = objnew
DO
  DO ipar = 2,NPL
    iwhich = (ipert-1)*NPL+ipar
    IF(ipert > objnew%k) iwhich = iwhich - 1
    !!
    !! Proposal from prior
    !!
    CALL RANDOM_NUMBER(ran_uni)
    objnew2%par(iwhich) = minlim(ipar) + maxpert(ipar)*ran_uni
  ENDDO
  ihamrej = 0
  iwhich = (ipert-1)*NPL+2
  IF(IHAM == 1)CALL COUPLE_HAM_BIRTH(objnew2,iwhich,ihamrej)
  IF(ihamrej == 0)EXIT
ENDDO
objnew = objnew2
END SUBROUTINE BIRTH
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

IF(obj%k == 0)THEN
  obj%h = 0._RP
  obj%z = 0._RP
ELSE
  idxz = 0
  idxz = (/1:obj%k*NPL:NPL/)
  obj%z = 0._RP
  obj%z(1:obj%k) = obj%par(idxz)
  obj%h = 0._RP
  obj%h(1) = obj%par(1)
  DO i = 2,obj%k
    obj%h(i) = obj%par(idxz(i))-obj%par(idxz(i-1))
  ENDDO
ENDIF
END SUBROUTINE CALCH
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

IF(rank == src)WRITE(ulog,*)'Reading particles from file ',samplefiletmp

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
      obj(ipart)%sdpar(1:NBAND) = 0.033_RP
   ENDIF
   CALL CALCH(obj(ipart))   !! Compute z-partition for object obj
ENDDO
CLOSE(33)
RETURN
END SUBROUTINE READ_PPD_P01
!!=======================================================================

SUBROUTINE READPARFILE()
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
CHARACTER(len=64):: filebasefile,parfile
INTEGER(KIND=IB) :: ubase,upar
INTERFACE
   FUNCTION LOGFACTORIAL(n)
     USE DATA_TYPE
     REAL(KIND=RP) :: LOGFACTORIAL
     REAL(KIND=RP),INTENT(IN):: n
   END FUNCTION LOGFACTORIAL
END INTERFACE
INTEGER :: ik

!!
!! Read parameter file
!!
filebasefile = 'filebase.txt'
OPEN(NEWUNIT(ubase),FILE=filebasefile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(ubase,*) filebasegl
READ(ubase,*) filebaselengl
CLOSE(ubase)

!!
!! Read parameter file
!!
parfile = filebasegl(1:filebaselengl) // '_parameters.dat'
OPEN(NEWUNIT(upar),FILE=parfile,FORM='formatted',STATUS='OLD',ACTION='READ')

READ(upar,*) particle_init_file
READ(upar,*) IMAP
READ(upar,*) ICOV
READ(upar,*) ENOS
READ(upar,*) IPOIPR
READ(upar,*) IAIS
READ(upar,*) IdB
READ(upar,*) ISETSEED
READ(upar,*) INODATA
READ(upar,*) ISPHER
READ(upar,*) IHAM
READ(upar,*) ISIM
READ(upar,*) subsmp
READ(upar,*) NPAVE
READ(upar,*) NBAL

READ(upar,*) NPL
READ(upar,*) NPING
READ(upar,*) NPART
READ(upar,*) NPARTAIS

READ(upar,*) MCMC_BI
READ(upar,*) MCMC_STEPS
READ(upar,*) NTEMP
READ(upar,*) NCOOLTRAJ
READ(upar,*) npercnt
READ(upar,*) NAP
READ(upar,*) NANG
READ(upar,*) NLMX
READ(upar,*) NBAND
READ(upar,*) NAVEF

READ(upar,*) IGA
READ(upar,*) frbw
frbw = 1._RP/frbw
READ(upar,*) FBW
ALLOCATE(bands(NBAND),sdtrgt(NBAND))
READ(upar,*) bands(1:NBAND)
READ(upar,*) sdtrgt(1:NBAND)

READ(upar,*) hmin
READ(upar,*) beta1gl
READ(upar,*) beta4gl
READ(upar,*) lambda

!!                   h      c         rt       alpha
!! minlim(1:4) = (/ hmin, 1450._RP, 1.25_RP, 0.00001_RP /)
!! maxlim(1:4) = (/ hmx,  1650._RP, 1.90_RP, 1.00000_RP /)
ALLOCATE(minlim(NPL),maxlim(NPL),maxpert(NPL),pertsd(NPL),pertsdsc(NPL))
minlim = 0._RP; maxlim = 0._RP; maxpert = 0._RP; pertsd = 0._RP
pertsdsc = 30._RP
READ(upar,*) minlim(2:NPL)
READ(upar,*) maxlim(2:NPL)
pertsdsc = 40._RP

!! Set prior and proposal for noise std dev.:
ALLOCATE(minlimsd(NBAND), maxlimsd(NBAND), maxpertsd(NBAND), &
         pertsdsd(NBAND), pertsdsdsc(NBAND))
READ(upar,*) minlimsd
READ(upar,*) maxlimsd
pertsdsdsc = 30._RP
maxpertsd  = maxlimsd-minlimsd
pertsdsd   = maxpertsd/pertsdsdsc
pertsdsdsc= 18._RP


!READ(20,*) NAVEF
!READ(20,*) NBAND
!NFREQ = NBAND*NAVEF
!ALLOCATE(bands(NBAND))
!READ(20,*) bands
!READ(20,*) NANG
CLOSE(upar)

!!
!! Some array dimensions:
!!
!NFPMX  = (NLMX * NPL) + (NPL-1)
!NLMX2  = NLMX
!NFPMX2 = (NLMX2 * NPL) + (NPL-1) !! Transition layer case
!NSDEVM = ((NLMX*NPL)+NPL-1)*NLMX

ALLOCATE(sdev(NBAND))         !! Standard devs
sdev = 0._RP
ALLOCATE(pingfilebase(NPING), particlefile(NPING))
kmin = 1
kmax = NLMX
NFPMX = (NLMX * NPL) + (NPL-1)
!kmin = 8
!kmax = 8

!! Poisson Prior on k:
ALLOCATE(pk(kmax))
DO ik = kmin,kmax
  pk(ik)  = EXP(-lambda)*lambda**ik/EXP(LOGFACTORIAL(REAL(ik,RP)))
ENDDO

RETURN
END SUBROUTINE READPARFILE
!=======================================================================

SUBROUTINE GETCHECKPOINT(particles,samplefiletmp)
!=======================================================================
!! Only master should read the sample. Particles are sent to slaves.
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: ipart,ustatus,usmp,baselen
REAL(KIND=RP),DIMENSION(NAP+NFPMX+NBAND):: sample_p01
TYPE(objstruc),DIMENSION(NPART):: particles
CHARACTER(len=64):: samplefiletmp
CHARACTER(len=64):: pingfilebasetmp
CHARACTER(len=64):: filebasechckpt

logfile         = filebasegl(1:filebaselengl) // '_SRJMH.log'
seedfile        = filebasegl(1:filebaselengl) // '_seeds.log'
mapfile         = filebasegl(1:filebaselengl) // '_map.dat'
repfile         = filebasegl(1:filebaselengl) // '_rep.dat'

!! Read checkpoint status from last run (gives ping No. to start with)
OPEN(NEWUNIT(ustatus),FILE='checkpoint/status.txt',FORM='formatted',STATUS='OLD',ACTION='READ')
READ(ustatus,*) icheckpoint
CLOSE(ustatus)
ipingst = icheckpoint
IF(rank == src)WRITE(ulog,*)'icheckpoint=',icheckpoint

IF(icheckpoint == 1)THEN
!! Set starting ping to last ping considered
  IF(rank == src)THEN
    OPEN(NEWUNIT(ulog),FILE=logfile,FORM='formatted',STATUS='REPLACE', &
    ACTION='WRITE',RECL=2048)
    WRITE(ulog,*) 'Starting Sequential RJMCMC sampling...'
    WRITE(ulog,*) '...running on ',NTHREAD,' cores'
  ENDIF
ELSE
  !! Set starting ping to last ping considered
  IF(rank == src)THEN
    OPEN(NEWUNIT(ulog),FILE=logfile,FORM='formatted',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
    IF(rank == src)WRITE(ulog,*) 'RJMCMC restarted after checkpoint on',NTHREAD,' cores for ping ',ipingst
  ENDIF
ENDIF

IF(icheckpoint == 1)THEN
  samplefiletmp = particle_init_file
ELSE
  WRITE(filebasechckpt,'(a,i3.3,a)')'p',ipingst,'_1000_2400'
  baselen = 14
  samplefiletmp = filebasechckpt(1:baselen) // '_particles.txt'
ENDIF
OPEN(UNIT=NEWUNIT(usmp),FILE=samplefiletmp,FORM='formatted',ACTION='READ')
IF(rank == src)WRITE(ulog,*)'Reading particles from file ',samplefiletmp,ipingst

DO ipart = 1,NPART
   READ(usmp,*) sample_p01
   particles(ipart)%k   = INT(sample_p01(4))
   particles(ipart)%logL= sample_p01(1)
   particles(ipart)%NFP = (particles(ipart)%k * NPL) + (NPL-1)
   particles(ipart)%par = 0._RP
   particles(ipart)%par(1:particles(ipart)%NFP) = sample_p01(5:4+particles(ipart)%NFP)
   particles(ipart)%sdpar = 0._RP
   IF(ICOV == 1)THEN
      particles(ipart)%sdpar(1:NBAND) = sample_p01(4+NFPMX+1:4+NFPMX+NBAND)
      particles(ipart)%sdpar(1:NBAND) = 0.033_RP
   ENDIF
   CALL CALCH(particles(ipart))   !! Compute z-partition for particlesect particles
ENDDO
CLOSE(usmp)
RETURN

END SUBROUTINE GETCHECKPOINT
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

!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
!IF(rank == src)WRITE(ulog,209) 'Loading data from file...',infile
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
len_snd  = 6+NFPMX+NBAND+2*NLMX+1
len_rcv  = 6+NFPMX+NBAND+2*NLMX+6

209 FORMAT(A26,A40)
RETURN
END SUBROUTINE READ_DATA
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
   CALL SOMMER_INTEGRAND(dat%angobs,fr,NBAND*NAVEF)
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
REAL(KIND=RP):: t1,t2,ran_nor

t1 = MPI_WTIME()
IF(INODATA == 0)THEN
  CALL LOGLHOOD(obj,dat)
ELSE
  CALL LOGLHOOD2(obj)
ENDIF
t2 = MPI_WTIME()
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
IF(rank == src)WRITE(ulog,*) 'time = ',t2-t1
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
   WRITE(ulog,201) obj%par((i-1)*NPL+1:i*NPL)
ENDDO
WRITE(ulog,202) '            ',obj%par(obj%k*NPL+1:obj%k*NPL+(NPL-1))
IF(ICOV == 1)THEN
   WRITE(ulog,*) 'SD parameters:'
   WRITE(ulog,203) obj%sdpar
ENDIF

201 FORMAT(4F12.4)
202 FORMAT(A12,4F12.4)
203 FORMAT(20F12.4)
END SUBROUTINE PRINTPAR
!=======================================================================

SUBROUTINE COUPLE_HAM(obj,ihamrej)
!=======================================================================
!!
!! Constyrain to Hamilton bounds
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: i,j,ncra
INTEGER(KIND=IB)                            :: ihamrej
TYPE (objstruc)                             :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+1)          :: idxc,idxr,idxa
REAL(KIND=RP),DIMENSION(kmax+1)             :: c,r,a
REAL(KIND=RP)                               :: cl,ch,al,ah

idxcra= 0
idxc  = 0
idxr  = 0
CALL GETIDXCRA(obj,idxcra,ncra)
idxc = idxcra(1:ncra:(NPL-1))
idxr = idxcra(2:ncra:(NPL-1))
idxa = idxcra(3:ncra:(NPL-1))

c = 0
r = 0
c(1:obj%k+1) = obj%par(idxc(1:obj%k+1))
r(1:obj%k+1) = obj%par(idxr(1:obj%k+1))
a(1:obj%k+1) = obj%par(idxa(1:obj%k+1))

!!
!! ioutside must not be set to 0 here!! Would destroy the 
!! perturbation setting from earlier (in PROPOSAL)
!!
DO i=1,obj%k+1
   cl=(1.54_RP-0.907_RP*r(i)+0.3695_RP*exp(1.88_RP*log(r(i))))*1.5004_RP*1000._RP
   ch=(1.6_RP-0.907_RP*r(i)+0.3695_RP*exp(2.01*log(r(i))))*1.5014_RP*1000._RP
   IF(cl < minlim(2))cl = minlim(2)
   IF(ch > maxlim(2))ch = maxlim(2)
   IF(c(i) > ch)THEN
     !ioutside = 1
     ihamrej = 1
     !IF(rank == src)WRITE(ulog,201) i,'ch',ioutside,c(i),ch
     CALL FLUSH(ulog)
   ENDIF
   IF(c(i) < cl)THEN
      !ioutside = 1
      ihamrej = 1
      !IF(rank == src)WRITE(ulog,201) i,'cl',ioutside,c(i),cl
   ENDIF
ENDDO

201 FORMAT(I3,a,I2,2F10.4)
END SUBROUTINE COUPLE_HAM
!=======================================================================

SUBROUTINE COUPLE_HAM_BIRTH(obj,iwhich,ihamrej)
!=======================================================================
!!
!! Constyrain to Hamilton bounds
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: i,j,ncra
INTEGER(KIND=IB)                            :: ihamrej,iwhich
TYPE (objstruc)                             :: obj
REAL(KIND=RP)                               :: c,r,a
REAL(KIND=RP)                               :: cl,ch,al,ah

c = obj%par(iwhich)
r = obj%par(iwhich+1)
a = obj%par(iwhich+2)
!!
!! ioutside must not be set to 0 here!! Would destroy the 
!! perturbation setting from earlier (in PROPOSAL)
!!
cl=(1.54_RP-0.907_RP*r+0.3695_RP*exp(1.88_RP*log(r)))*1500.4_RP!*1000._RP
ch=(1.6_RP-0.907_RP*r+0.3695_RP*exp(2.01_RP*log(r)))*1501.4_RP!*1000._RP
IF(cl < minlim(2))cl = minlim(2)
IF(ch > maxlim(2))ch = maxlim(2)
IF(c > ch .OR. c < cl)THEN
  ihamrej = 1
ENDIF

201 FORMAT(I3,a,I2,2F10.4)
END SUBROUTINE COUPLE_HAM_BIRTH
!=======================================================================

SUBROUTINE PRIOR(obj,dat)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=IB) :: i
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP)    :: ran_uni

DO i=1,obj%NFP
   CALL RANDOM_NUMBER(ran_uni)
   obj%par(i) = minlim(i) + maxpert(i) * ran_uni
ENDDO
IF(INODATA == 0)THEN
  CALL LOGLHOOD(obj,dat)
ELSE
  CALL LOGLHOOD2(obj)
ENDIF
RETURN
END SUBROUTINE PRIOR

!=======================================================================

SUBROUTINE SOMMER_INTEGRAND(thd,freq,nfrq)
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
INTEGER(KIND=IB) :: nfrq, Nmx, Nrmx, Nr_fr
REAL(KIND=RP) :: w,exp_arg,umax
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
!imx = MAXLOC(N,DIM=1)
!Nrmx = Nr(imx)
!N = Nmx
!Nr = Nrmx
!drTh = drTh(imx)
!diTh = diTh(imx)
ALLOCATE( btR(Nmx,NANG,nfrq),rTh(Nmx,nfrq) )
btR = 0._RP
rTh = 0._RP

DO kf=1,nfrq
   ALLOCATE( tmp2(N(kf),2),exp_termR(N(kf)),kr_sinRTh(N(kf),NANG) )
   tmp2 = 0._RP
   exp_termR = 0._RP
   kr_sinRTh = 0._RP
   rTh(1:Nr(kf),kf) = (/ ((i*drTh(kf))-drTh(kf) ,i=1,Nr(kf)) /)
   rTh(Nr(kf)+1:N(kf),kf) = PI/2._RP-CMPLX(0._RP,1._RP,RP)*(/ (i*drTh(kf),i=1,N(kf)-Nr(kf)) /)
   rTh(1:N(kf),kf) = rTh(1:N(kf),kf) + TINY(1._RP)

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

   tmp2(:,1) = SIN(rTh(1:N(kf),kf))
   tmp2(:,2) = 0._RP            !! Had to add second line here to allow for use of MATMUL...
   kr_sinRTh=k(kf)*MATMUL(tmp2,rr)

   ! Bessel function first kind, zeroth order (far field approx. and then exact near field)
   btR(1:N(kf),:,kf) = SQRT(2._RP/(PI*kr_sinRTh))*COS(kr_sinRTh-0.25_RP*PI) !large arg approx >4;
   !!
   !!  ??????BESSEL FUNCTION FOR COMPLEX VARIABLE??????
   !!
   DO i = 1,NANG
   DO j = 1,N(kf)
      IF(ABS(kr_sinRTh(j,i)) < 4._RP)THEN 
         btR(j,i,kf)=dbesj0(REAL(kr_sinRTh(j,i),RP))
      ENDIF
   ENDDO
   ENDDO

   exp_termR=EXP(CMPLX(0._RP,1._RP,RP)*k(kf)*(z_t)*COS(rTh(1:N(kf),kf)))*SIN(rTh(1:N(kf),kf))
   tmp = m*exp_termR(1:Nr(kf))
   tmpi = mi*exp_termR(Nr(kf):N(kf))

   DO i = 1,NANG
      btR(1:Nr(kf),i,kf)=btR(1:Nr(kf),i,kf)*tmp
      btR(Nr(kf):N(kf),i,kf)=btR(Nr(kf):N(kf),i,kf)*tmpi
   ENDDO

   DEALLOCATE( tmp2,m,mi,tmp,tmpi,kr_sinRTh,exp_termR )
ENDDO
RETURN
END SUBROUTINE SOMMER_INTEGRAND
!=======================================================================

SUBROUTINE SPH_REF_NLAY(thd,m_rg,freq,ref_sph,nfrq,NFP,idx)
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
INTEGER(KIND=IB) :: nfrq, NFP,N_tmp
REAL(KIND=RP)    :: dbesj0
REAL(KIND=RP)    :: ASINH
REAL(KIND=RP), DIMENSION(nfrq) :: freq
REAL(KIND=RP), DIMENSION(NFP)  :: m_rg
REAL(KIND=RP), DIMENSION(NANG) :: thd
COMPLEX(KIND=RP),DIMENSION(NANG) :: rRP,iRP
REAL(KIND=RP),DIMENSION(nfrq,NANG) :: ref_sph
COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE :: ref,ref2,rTh_tmp

DO kf=1,nfrq
   lf = (idx-1)*nfrq+kf
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
         rTh_tmp(j) = rTh(i,lf)
         j = j + 1
      ENDDO
      rTh_tmp(N_tmp-Ni:N_tmp) = rTh(N(lf)-Ni:N(lf),lf)

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
      CALL REF_NLAY4(90._RP-rTh(1:N(lf),lf)*180._RP/PI,m_rg,freq(kf),ref2,cw,rw,1,N(lf),NFP)
   ENDIF

   rRP = drTh(lf)/3._RP*MATMUL(ref2(1:Nr(lf)),btR(1:Nr(lf),:,lf))
   iRP = diTh(lf)/3._RP*MATMUL(ref2(Nr(lf):N(lf)),btR(Nr(lf):N(lf),:,lf))

   ref_sph(kf,:) = ABS(CMPLX(0._RP,1._RP,RP)*(rRP+iRP)/RPs(:,lf))
   DEALLOCATE( ref,ref2,rTh_tmp )
ENDDO
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
!REAL(KIND=RP):: t1,t2

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

SUBROUTINE REF_NLAY44(thd,m_rg,freq,ref,cw,rw,nfrq,NANG,NFP)
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
REAL(KIND=RP),DIMENSION(nfrq)            :: freq
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: freq2
REAL(KIND=RP),DIMENSION(NFP)  :: m_rg
COMPLEX(KIND=RP), DIMENSION(NANG) :: thd,th1,th1tmp
COMPLEX(KIND=RP), DIMENSION(nfrq,NANG) :: ref
COMPLEX(KIND=RP), DIMENSION(NANG,nfrq) :: reftmp,z1
REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: c,alf,r,d
COMPLEX(KIND=RP), DIMENSION(:), ALLOCATABLE :: v
COMPLEX(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: zz,th,v2,k
COMPLEX(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: z
COMPLEX(KIND=RP), DIMENSION(NANG,nfrq) :: znlay,znlay1,zin,zm1
COMPLEX(KIND=RP), DIMENSION(nfrq,NANG) :: tankdtmp,tankd,thtmp,ktmp
COMPLEX(KIND=RP) :: CACOS,CTAN
LOGICAL :: ISNAN

NLAY = ((SIZE(m_rg,1)-3)/4)+1
NLAY2 = NLAY+1

ALLOCATE(c(NLAY2),alf(NLAY2),r(NLAY2),d(NLAY-1))
ALLOCATE(v(NLAY2),v2(NLAY2,nfrq),zz(NLAY2,NANG),th(NLAY2,NANG))
ALLOCATE(z(NLAY2,NANG,nfrq),freq2(NLAY2,nfrq),k(nfrq,NLAY2))

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

!ref = ABS(TRANSPOSE(reftmp))
ref = TRANSPOSE(reftmp)
!at 0 degrees, the reflection coefficeint is -1
DO j=1,nfrq
   DO i=1,NANG
      IF(ISNAN(ABS(ref(j,i))) == .TRUE.) ref(j,i)=-1._RP
   ENDDO
ENDDO
DEALLOCATE(c,alf,r,d,v,v2,zz,th,z,freq2,k)
RETURN
END SUBROUTINE REF_NLAY44

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
ELSE
    vp=1._RP/(1._RP/cr + Q2*alfr/(2._RP*PI*fr)**y)
    alf1dB= alfrdB*(f1/fr)**(y-1._RP)
ENDIF

RETURN
END SUBROUTINE KK_WATERS_dBmkHz
!=======================================================================

SUBROUTINE SND_PARTICLE(obj,ipart,idest,idone,irealloc)
!=======================================================================
!!
!! Exchanging particles
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB):: ipart,idest,idone,irealloc
REAL(KIND=RP):: t1,t2
TYPE(objstruc):: obj
REAL(KIND=RP),DIMENSION(len_snd):: obj_tmp

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
         CALL MPI_SSEND( beta, NTEMP1, MPI_DOUBLE_PRECISION, idest, &
                         ipart, MPI_COMM_WORLD, ierr )
      ENDIF
      CALL MPI_SSEND( logLtrgt, 1, MPI_DOUBLE_PRECISION, idest, &
                      ipart, MPI_COMM_WORLD, ierr )
      obj_tmp = 0._RP
      obj_tmp =  (/ obj%logL, obj%logwt, obj%logPr, obj%lognorm, REAL(obj%k,RP), REAL(obj%NFP,RP), &
                    obj%par, obj%sdpar, obj%z, obj%h, REAL(obj%nfail,RP) /)
      CALL MPI_SSEND( obj_tmp, len_snd, MPI_DOUBLE_PRECISION, idest, &
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
      CALL MPI_RECV(obj_tmp, len_snd, MPI_DOUBLE_PRECISION, src,&
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
      obj%h          = 0._RP
      obj%h          = obj_tmp(7+NFPMX+NBAND+NLMX:7+NFPMX+NBAND+2*NLMX-1)
      obj%nfail      = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX))
   ENDIF
ENDIF !  MPI
RETURN
END SUBROUTINE SND_PARTICLE
!=======================================================================

SUBROUTINE RCV_PARTICLE(obj,isource,ipart,ifail)
!=======================================================================
!!
!! Exchanging particles
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB):: isource,ipart,ifail
REAL(KIND=RP):: logLG
REAL(KIND=RP):: LOGPLUS
REAL(KIND=RP):: t1,t2
TYPE(objstruc):: obj
REAL(KIND=RP),DIMENSION(len_rcv):: obj_tmp

!!
!!  Sending samples to master
!!
IF(rank == src)THEN
   obj_tmp = 0._RP
   CALL MPI_RECV(obj_tmp, len_rcv, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
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
   obj%h          = 0._RP
   obj%h          = obj_tmp(7+NFPMX+NBAND+NLMX:7+NFPMX+NBAND+2*NLMX-1)
   acc_ratio      = obj_tmp(7+NFPMX+NBAND+2*NLMX)
   iaccept_bd     = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX+1))
   ireject_bd     = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX+2))
   i_bd           = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX+3))
   i_zpert        = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX+4))
   obj%nfail      = INT(obj_tmp(7+NFPMX+NBAND+2*NLMX+5))
ELSE
   obj_tmp = 0._RP
   obj_tmp =  (/ obj%logL, obj%logwt, obj%logPr, obj%lognorm, REAL(obj%k,RP), REAL(obj%NFP,RP), &
                 obj%par , obj%sdpar , obj%z, obj%h, REAL(iaccept,RP)/REAL(ireject,RP),REAL(iaccept_bd,RP),&
                 REAL(ireject_bd,RP),REAL(i_bd,RP),REAL(i_zpert,RP),REAL(rank,RP),&
                 REAL(obj%nfail,RP) /)
   CALL MPI_SSEND( obj_tmp, len_rcv, MPI_DOUBLE_PRECISION, src, &
                   ipart, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( tcmp, 1, MPI_DOUBLE_PRECISION, src, &
                   ipart, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( ifail, 1, MPI_INTEGER, src, &
                   ipart, MPI_COMM_WORLD, ierr )
ENDIF !  MPI
RETURN
END SUBROUTINE RCV_PARTICLE
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

WRITE(ulog,*) 'Global best model:'
CALL PRINTPAR(obj)
WRITE(ulog,*) 'Global best logL = ',obj%logL
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
   WRITE(ulog,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(ulog,*) '~~~                                                        ~~~  '
   WRITE(ulog,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
   WRITE(ulog,*) '~~~                                                        ~~~  '
   WRITE(ulog,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(ulog,*) '...running on ',NTHREAD,' cores'
   WRITE(ulog,*) ''
   WRITE(ulog,210) 'ICOV                 :   ',ICOV
   WRITE(ulog,210) 'ENOS                 :   ',ICOV
   WRITE(ulog,210) 'IPOIPR               :   ',IPOIPR
   WRITE(ulog,210) 'Number of angles     :   ',NANG
   WRITE(ulog,210) 'Number of frequencies:   ',NBAND
   WRITE(ulog,210) 'NAVEF                :   ',NAVEF
   WRITE(ulog,210) 'IGA                  :   ',IGA
   WRITE(ulog,210) 'ISPHER               :   ',ISPHER
   WRITE(ulog,210) 'IAIS                 :   ',IAIS
   WRITE(ulog,210) 'INODATA              :   ',INODATA
   WRITE(ulog,210) 'ISETSEED             :   ',ISETSEED
   WRITE(ulog,210) 'subsmp               :   ',subsmp
   WRITE(ulog,210) 'NPAVE                :   ',NPAVE
   !WRITE(ulog,209) 'Particle files:           ',particlefile
   WRITE(ulog,*) ''
   WRITE(ulog,*) 'NPART    = ',NPART
   WRITE(ulog,*) 'NPARTAIS = ',NPARTAIS
   WRITE(ulog,*) 'MCMC_BI  = ',MCMC_BI
   WRITE(ulog,*) 'MCMC_STEP= ',MCMC_STEPS
   WRITE(ulog,*) 'NTEMP    = ',NTEMP
   WRITE(ulog,*) 'NBAL     = ',NBAL 
   WRITE(ulog,*) 'NFPMX    = ',NFPMX
   WRITE(ulog,*) 'kmin     = ',kmin
   WRITE(ulog,*) 'kmax     = ',kmax
   WRITE(ulog,*) 'hmin     = ',hmin
   WRITE(ulog,*) 'minlim:  '
   WRITE(ulog,201) minlim
   WRITE(ulog,*) 'maxlim:  '
   WRITE(ulog,201) maxlim
   WRITE(ulog,*) 'pertsdsc:'
   WRITE(ulog,201) pertsdsc
   WRITE(ulog,*) 'subsmp = ',subsmp
   WRITE(ulog,*) ''
   WRITE(ulog,*) 'Done reading data.'
   WRITE(ulog,*) ''
   IF(IGA == 0)THEN
      WRITE(ulog,*)' !!!   NARROW BAND INTENSITY FREQ AVERAGING   !!!'
      WRITE(ulog,*)' FBW = ',FBW
   ELSEIF(IGA == 1)THEN
      WRITE(ulog,*)' !!!   GAUSSIAN FREQ AVERAGING   !!!'
      WRITE(ulog,*)' frbw = ',frbw
   ELSEIF(IGA == 2)THEN
      WRITE(ulog,*)' !!!   1/3 octave INTENSITY FREQ AVERAGING   !!!'
   ENDIF
   WRITE(ulog,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(ulog,*) ''
   CALL FLUSH(ulog)
ENDIF
RETURN
201 FORMAT(200F12.4)
209 FORMAT(A30,A64)
210 FORMAT(A26,I4)
END SUBROUTINE WRITE_LOG
!=======================================================================

SUBROUTINE MAKE_BETA(beta1,beta4)
!!
!! Make cooling schedule (three geometrically spaced legs)
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: itemp
REAL(KIND=RP)                  :: beta1,beta2,beta3,beta4
REAL(KIND=RP),DIMENSION(NTEMP1):: logbeta

beta2 = beta1+0.2_RP*ABS(beta1-beta4)
beta3 = beta1+0.75_RP*ABS(beta1-beta4)
!IF(rank == src)THEN
!   WRITE(ulog,*)'Computing cooling schedule (inverse temperature,geometrical spacing):'
!   WRITE(ulog,*)'beta 1 =',beta1,'1 to ',NTEMP/20
!   WRITE(ulog,*)'beta 2 =',beta2,NTEMP/20+1,' to ',NTEMP/4
!   WRITE(ulog,*)'beta 3 =',beta3,NTEMP/4+1,' to ',NTEMP
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
!     WRITE(ulog0,*) beta(itemp)
!   ENDDO 
!   CALL FLUSH(60)
!   CLOSE(60)
!ENDIF
!STOP

209 FORMAT(50ES16.8)
RETURN
END SUBROUTINE MAKE_BETA
!==============================================================================

SUBROUTINE AIS(obj,dat)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,itemp,imcmc
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
!      idxrand(1:obj%NFP) = RANDPERM(obj%NFP)
!      ipar = 1
      CALL EXPLORE_MHAIS(obj,dat,beta(itemp),logPLratio)
      logwt = logwt - logPLratio
   ENDDO

ENDDO
obj%logwt = logwt + obj%logP + obj%logL
CALL CALCH(obj)   !! Compute z-partition for AIS object obj
RETURN
END SUBROUTINE AIS
!==============================================================================

SUBROUTINE EXPLORE_MHAIS(obj,dat,beta_mh,logPLratio)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,idel,ipar,ipar2,ilay,idxz,i_comp_pred
INTEGER(KIND=IB)                            :: ncra
TYPE(objstruc)                              :: obj,objnew1
TYPE (datastruc)                      :: dat
REAL(KIND=RP)                               :: logPLratio,logPratio,logy,logq1_1,logq1_2
REAL(KIND=RP)                               :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                               :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp
REAL(KIND=RP)                               :: logarp
INTEGER(KIND=IB)                            :: arptype

objnew1 = obj
!!
!! Do Metropolis-Hastings
!!
!IF(1==2)THEN
   !!
   !! If k did not change, carry out MH sweep
   !!
   idxrand = 0
   idxrand(1:obj%NFP) = RANDPERM(obj%NFP)
   !!
   !! Do Metropolis-Hastings update on c, rho, alpha
   !!
   DO ipar = 1,obj%NFP
     iwhich = idxrand(ipar)
     ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
     ipar2  = iwhich-(ilay-1)*NPL
     IF(ilay > obj%k) ipar2 = ipar2 + 1

     !ipropose(ipar2) = ipropose(ipar2) + 1
     !CALL PROPOSAL(obj,objnew1,iwhich,1._RP)
     CALL PROPOSAL_AIS(obj,objnew1,dat,beta_mh,iwhich)
     CALL CHECKBOUNDS(objnew1)
     IF(ioutside == 0)THEN
       IF(INODATA == 0)THEN
         i_comp_pred = 1
         CALL LOGLHOOD(objnew1,dat)
       ELSE
         CALL LOGLHOOD2(objnew1)
       ENDIF
       logPratio  = objnew1%logPr
       logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh
       !logPLratio = (objnew1%logL - obj%logL)*beta_mh
       CALL RANDOM_NUMBER(ran_uni)
       IF(ran_uni >= EXP(logPLratio))THEN
         objnew1 = obj
       ELSE
         obj = objnew1
         !iaccept(ipar2) = iaccept(ipar2) + 1
       ENDIF
     ELSE !! outside before delayed rejection
       objnew1 = obj
       ioutside = 0
     ENDIF
   ENDDO
!   ipropose_buf(:,ibuf,ic) = ipropose
!   iaccept_buf(:,ibuf,ic)  = iaccept
!   IF(ibuf == NBUF) ibuf = 0
!   ibuf = ibuf + 1

   !!
   !! Do Metropolis-Hastings on data-error standard deviations
   !!
   IF(ICOV == 1)THEN
     DO ipar = 1,NBAND
        !ipropose_sd = ipropose_sd + 1
        CALL PROPOSAL_SD(obj,objnew1,ipar)
        IF(ioutside == 0)THEN
           IF(INODATA == 0)THEN
             i_comp_pred = 0
             CALL LOGLHOOD(objnew1,dat)
           ELSE
             CALL LOGLHOOD2(objnew1)
           ENDIF
           logPLratio = (objnew1%logL - obj%logL)*beta_mh
           CALL RANDOM_NUMBER(ran_uni)
           IF(ran_uni >= EXP(logPLratio))THEN
              objnew1 = obj
           ELSE
              obj = objnew1
              !iaccept_sd = iaccept_sd + 1
           ENDIF
        ELSE
           objnew1 = obj
           ioutside = 0
        ENDIF
     ENDDO
   ENDIF ! ICOV if
END SUBROUTINE EXPLORE_MHAIS

!==============================================================================

!SUBROUTINE EXPLORE_MHAIS(obj,logPLratio,beta_mh,dat)
!==============================================================================
!USE MPI
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: iwhich,NFPnew,idel,ipar,ilay,idxz,infloop
!INTEGER(KIND=IB) :: ncra
!TYPE(objstruc) :: obj,objnew
!TYPE(datastruc):: dat
!REAL(KIND=RP)  :: logPLratio,ran_uni,ran_uni_Z,ran_uni_sd
!REAL(KIND=RP)  :: znew,beta_mh
!INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
!INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
!REAL(KIND=RP)   ,DIMENSION(obj%k)           :: ztmp
!
!objnew = obj
!!! Draw uniform Birth-Death probability
!CALL RANDOM_NUMBER(ran_uni_Z)
!
!!! Perturb z with 0.5 probability
!IF((ran_uni_Z > 0.5_RP))THEN
!   i_zpert = 2
!   i_bd = 0
!   idxrand = 0
!   idxrand(1:obj%k) = RANDPERM(obj%k)
!   idxz = idxrand(1)
!   iwhich = ((idxz-1)*NPL)+1
!   infloop = 0
!   CALL PROPOSALZ(obj,objnew,idxz,iwhich)
!ENDIF

!!
!! Do Metropolis-Hastings
!!
!! Perturb only c, rho, alpha, and sigma
!CALL GETIDXCRA(obj,idxcra,ncra)
!idxrand = 0
!idxrand(1:ncra) = RANDPERM(ncra)
!DO ipar = 1,ncra
!   iwhich = idxcra(idxrand(ipar))
!   CALL PROPOSAL_AIS(obj,objnew,dat,beta_mh,iwhich)
!   IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)
!   CALL CHECKBOUNDS(objnew)
!   IF(ioutside == 0)THEN
!      CALL LOGLHOOD(objnew,dat)
!!      logPLratio = (objnew%logL - obj%logL)*beta_mh
!      logPLratio = (objnew%logP - obj%logP) + (objnew%logL - obj%logL)*beta_mh
!      CALL RANDOM_NUMBER(ran_uni)
!      IF(ran_uni >= EXP(logPLratio))THEN
!         logPLratio = 0._RP
!         objnew = obj
!         ireject = ireject + 1
!      ELSE
!         obj = objnew
!         iaccept = iaccept + 1
!      ENDIF
!   ELSE
!      logPLratio = 0._RP
!      objnew = obj
!      ireject = ireject + 1
!      ioutside = 0
!   ENDIF 
!ENDDO
!!
!! Do Metropolis-Hastings on data-error standard deviations
!!
!IF(ICOV == 1)THEN
!   !! Perturb std devs with .25 probability
!   CALL RANDOM_NUMBER(ran_uni_sd)
!   IF(ran_uni_sd>=0.25_RP)THEN
!      DO ipar = 1,NBAND
!         CALL PROPOSAL_SD(obj,objnew,ipar)
!         IF(ioutside == 0)THEN
!            CALL LOGLHOOD(objnew,dat)
!!            logPLratio = (objnew%logL - obj%logL)*beta_mh
!            logPLratio = (objnew%logP - obj%logP) + (objnew%logL - obj%logL)*beta_mh
!            CALL RANDOM_NUMBER(ran_uni)
!            IF(ran_uni >= EXP(logPLratio))THEN
!               logPLratio = 0._RP
!               objnew = obj
!               ireject = ireject + 1
!            ELSE
!               obj = objnew
!               iaccept = iaccept + 1
!            ENDIF
!         ELSE
!            logPLratio = 0._RP
!            objnew = obj
!            ireject = ireject + 1
!            ioutside = 0
!         ENDIF
!      ENDDO
!   ENDIF
!ENDIF
!RETURN
!END SUBROUTINE EXPLORE_MHAIS
!=======================================================================

SUBROUTINE RESAMPLE_AIS(particles,particles_new,N1,N2)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)             :: ik,ipart,idxpart,N1,N2
TYPE(objstruc),DIMENSION(N1) :: particles
TYPE(objstruc),DIMENSION(N2) :: particles_new
REAL(KIND=RP),DIMENSION(N1)  :: cumsum     ! Cumulative distribution to pick random according to wt
REAL(KIND=RP)                :: ran_uni
REAL(KIND=RP)                :: LOGPLUS
INTEGER(KIND=IB),DIMENSION(NLMX) :: kcount
REAL(KIND=IB),DIMENSION(NLMX) :: logPok       ! Posterior on k

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Find no. particles per k
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
kcount = 0
DO ik = 1,NLMX
   DO ipart = 1,N1
      IF(particles(ipart)%k == ik)THEN
         kcount(ik) = kcount(ik)+1
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO
!! Noralized log posterior marginal on k
logPok = LOG(REAL(kcount,RP))-LOG(REAL(N1))

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Resample new set of particles according to weight
!!
!! Eric: Something is wrong with the weights in AIS! Changing resampling to logL 
!! here for now!
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
cumsum  = 0._RP
!cumsum(1) = (particles(1)%logwt+logPok(particles(1)%k))
cumsum(1) = (particles(1)%logL+logPok(particles(1)%k))
DO ipart = 2,N1
  !! Compute cumulative logwt distribution
  !cumsum(ipart) = LOGPLUS(cumsum(ipart-1),particles(ipart)%logwt+logPok(particles(ipart)%k))
  cumsum(ipart) = LOGPLUS(cumsum(ipart-1),particles(ipart)%logL+logPok(particles(ipart)%k))
ENDDO ! PARTICLE LOOP
!! Normalize to 1
cumsum = cumsum - cumsum(N1)
cumsum = EXP(cumsum)

DO ipart = 1,N2
  !! Draw randomly according to weight using cumulative dist:
  CALL RANDOM_NUMBER(ran_uni)
  idxpart = MINLOC(ABS(cumsum-ran_uni),DIM=1)
  particles_new(ipart) = particles(idxpart)
  particles_new(ipart)%logwt = LOG(1._RP/REAL(N2,RP))
ENDDO

RETURN
END SUBROUTINE RESAMPLE_AIS
!=======================================================================

SUBROUTINE RESAMPLE_PER_K(particles,particles_new,kcount,kidx)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                      :: iipart,ipart,ik,idxpart,idxk
INTEGER(KIND=IB),DIMENSION(NLMX)      :: kcount
INTEGER(KIND=IB),DIMENSION(NPART,NLMX):: kidx
TYPE(objstruc),DIMENSION(NPART)       :: particles,particles_new
REAL(KIND=RP),DIMENSION(NPARTAIS,NLMX):: cumsum     ! Cumulative distribution to pick random according to wt
REAL(KIND=RP),DIMENSION(NLMX)         :: cumsumk    ! Cumulative distribution to pick random according to wt
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
ENDDO
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Resample new set of particles according to weight
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
cumsum  = 0._RP
DO ik = 1,NLMX
   IF(kcount(ik) >= 1)THEN
      cumsum(1,ik) = (particles(kidx(1,ik))%logwt)
      DO ipart = 2,kcount(ik)
         !! Compute cumulative logwt distribution
         cumsum(ipart,ik) = LOGPLUS(cumsum(ipart-1,ik),particles(kidx(ipart,ik))%logwt)
         !cumsum(ipart,ik) = cumsum(ipart-1,ik)+exp(particles(kidx(ipart,ik))%logwt)
      ENDDO ! PARTICLE LOOP
      !! Normalize to 1
      cumsum(:,ik) = cumsum(:,ik) - cumsum(kcount(ik),ik)
      cumsum(:,ik) = EXP(cumsum(:,ik))
      !cumsum(:,ik) = cumsum(:,ik)/cumsum(kcount(ik),ik)
      WRITE(81,*) ik,kcount(ik)
      DO ipart = 1,kcount(ik)
         WRITE(81,*) cumsum(ipart,ik)
      ENDDO

   ENDIF
ENDDO
cumsumk = 0._RP
cumsumk(1) = kcount(1)
DO ik = 2,NLMX
   cumsumk(ik) = cumsumk(ik-1)+kcount(ik)
ENDDO
cumsumk = cumsumk/cumsumk(NLMX)

DO ipart = 1,NPART
   !! Draw randomly according to weight using cumulative dist:
   CALL RANDOM_NUMBER(ran_unik)
   CALL RANDOM_NUMBER(ran_uni)
   idxk = MINLOC(ABS(FLOOR(cumsumk-ran_unik)),DIM=1)
   idxpart = MINLOC(ABS(cumsum(1:kcount(particles(ipart)%k),idxk)-ran_uni),DIM=1)
   iipart = kidx(idxpart,idxk)
   particles_new(ipart) = particles(iipart)
   particles_new(ipart)%logwt = LOG(1._RP/REAL(NPART,RP))
ENDDO

RETURN
END SUBROUTINE RESAMPLE_PER_K
!!=======================================================================

SUBROUTINE RESAMPLE_ALL_K(particles,particles_new,N1,N2)
!! Resample from N1 to N2 particles according to logwt
!!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB)             :: ik,ipart,idxpart,N1,N2
TYPE(objstruc),DIMENSION(N1) :: particles
TYPE(objstruc),DIMENSION(N2) :: particles_new
REAL(KIND=RP),DIMENSION(N1)  :: cumsum     ! Cumulative distribution to pick random according to wt
REAL(KIND=RP)                :: ran_uni
REAL(KIND=RP)                :: LOGPLUS
INTEGER(KIND=IB),DIMENSION(NLMX) :: kcount
REAL(KIND=IB),DIMENSION(NLMX) :: logPok       ! Posterior on k

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Find no. particles per k
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
kcount = 0
DO ik = 1,NLMX
   DO ipart = 1,N1
      IF(particles(ipart)%k == ik)THEN
         kcount(ik) = kcount(ik)+1
      ENDIF
   ENDDO ! PARTICLE LOOP
ENDDO
!! Noralized log posterior marginal on k
logPok = LOG(REAL(kcount,RP))-LOG(REAL(N1))

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!! Resample new set of particles according to weight
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
cumsum  = 0._RP
!PRINT*,'rank=',rank,'k=',particles(1)%k,iping,N1,N2
!PRINT*,'k=',particles_new(1)%k
cumsum(1) = (particles(1)%logwt+logPok(particles(1)%k))
DO ipart = 2,N1
  !! Compute cumulative logwt distribution
  cumsum(ipart) = LOGPLUS(cumsum(ipart-1),particles(ipart)%logwt+logPok(particles(ipart)%k))
ENDDO ! PARTICLE LOOP
!! Normalize to 1
cumsum = cumsum - cumsum(N1)
cumsum = EXP(cumsum)

!DO ipart = 1,N1
!  WRITE(77,*) particles(ipart)%logL,logPok(particles(ipart)%k)
!  WRITE(81,*) cumsum(ipart)
!ENDDO
!STOP

DO ipart = 1,N2
  !! Draw randomly according to weight using cumulative dist:
  CALL RANDOM_NUMBER(ran_uni)
  idxpart = MINLOC(ABS(cumsum-ran_uni),DIM=1)
  particles_new(ipart) = particles(idxpart)
  particles_new(ipart)%logwt = LOG(1._RP/REAL(N2,RP))
ENDDO

RETURN
END SUBROUTINE RESAMPLE_ALL_K
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

ASINH = LOG(x+SQRT(x**2+1))

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

Function LOGPLUS(x,y)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
REAL(KIND=RP)    :: logplus,x,y

IF(x > y)THEN
   LOGPLUS = x+LOG(1._RP+EXP(y-x))
ELSE
   LOGPLUS = y+log(1._RP+EXP(x-y))
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
!   IF(rank==src)WRITE(ulog,*) 'Master seed:',iseed1
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
!!=======================================================================
RECURSIVE FUNCTION LOGFACTORIAL(n)  RESULT(fact)
!-----Factorial------------------------------------------------------
!!=======================================================================

USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: fact
REAL(KIND=RP), INTENT(IN) :: n

IF (n == 0) THEN
   fact = 0
ELSE
   fact = LOG(n) + LOGFACTORIAL(n-1)
END IF

END FUNCTION LOGFACTORIAL
!=======================================================================
! This is the end my fiend...
! EOF
