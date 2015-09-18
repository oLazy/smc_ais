MODULE RJMCMC_COM
  USE MPI
  USE DATA_TYPE
  IMPLICIT NONE
  
  !
  ! General switches
  !
  INTEGER(KIND=IB) :: IMAP       !! WRITE REPLICA AND EXIT
  INTEGER(KIND=IB) :: ICOV       !! 1 = Sample over sigma
  !! 2 = Use sigma from ping ave
  !! 3 = Use sigma from ping ave but scale by factor (free parameter)
  INTEGER(KIND=IB) :: ENOS       !! 0 = uniform prior on k
  !! 1 = Apply even-numbered order statistics as prior on k (avoids very thin layers)
  INTEGER(KIND=IB) :: IPOIPR     !! 1 = Apply Poisson prior on k
  INTEGER(KIND=IB) :: IAIS       !! 0 = Apply rjMCMC on small No. particles; 1 = Apply AIS; 
  INTEGER(KIND=IB) :: IdB        !! 1=Carry out computation in dB
  INTEGER(KIND=IB) :: ISETSEED   !! Fix the random seed 
  INTEGER(KIND=IB) :: INODATA    !! Sample prior without data
  INTEGER(KIND=IB) :: ISPHER     !! Spherical refl. coeff. 
  INTEGER(KIND=IB) :: IHAM       !! Hamilton prior
  INTEGER(KIND=IB) :: ISIM       !! WRITE REPLICA WITH NOISE AND EXIT
  INTEGER(KIND=IB) :: subsmp     !! Subsample Sommerfeld integral plane wave part
  INTEGER(KIND=IB) :: NPAVE      !! Skip pings
  INTEGER(KIND=IB) :: IGA        !! Type of averaging (0=intensity, 1=Gaussian, 2=1/3 octave intensity)

!!
!!  AUV DATA RJMCMC trial ping
!!
   INTEGER(KIND=IB):: NPL         !! Number parameters per layer
   INTEGER(KIND=IB):: NPING       !! Number of pings
   INTEGER(KIND=IB):: NPART       !! Number of particles
   INTEGER(KIND=IB):: NPARTAIS    !! Number of particles to perform AIS on
   INTEGER(KIND=IB):: NANG    
   INTEGER(KIND=IB):: NLMX    
   INTEGER(KIND=IB):: NBAND   
   INTEGER(KIND=IB):: NAVEF        !! No. freq per band
   INTEGER(KIND=IB):: icheckpoint  !! Checkpointing variable (needed on some computer clusters)

   INTEGER(KIND=IB):: ipingst      !! Ping # to start at (read from ./checkpoint/status.txt)
   REAL(KIND=RP)   :: frbw         !! Frac. bandwidth for freq ave.
   REAL(KIND=RP),ALLOCATABLE, DIMENSION(:):: bands
   REAL(KIND=RP),ALLOCATABLE, DIMENSION(:):: sdtrgt
   REAL(KIND=RP) :: FBW

   INTEGER(KIND=IB)                  :: filebaselen, filebaselengl
   CHARACTER(len=64)                 :: filebasegl, particle_init_file
   CHARACTER(len=64),ALLOCATABLE,DIMENSION(:):: pingfilebase, particlefile

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pk       ! Poisson prior on k
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlim
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlim
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpert
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsc
   INTEGER(KIND=IB)            :: kmin               ! Min number of layers
   INTEGER(KIND=IB)            :: kmax               ! Max number of layers
   REAL(KIND=RP)               :: lambda             ! Lambda parameter for Poisson prior on k
   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness
   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: fr,fstep   ! Total frequency array
   REAL(KIND=RP)               :: z_t,cw,rw,hmx
   REAL(KIND=RP)               :: logLtrgt
   INTEGER(KIND=IB)            :: kmintmp,kmaxtmp
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(len=64) :: repfile
   INTEGER(kind=IB) :: ulog
!!
!!  Standard deviation prior variables:
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlimsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdsc
!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX                             !! = (NLMX * NPL) + (NPL-1)
   INTEGER(KIND=IB)           :: len_snd,len_rcv 
   INTEGER(KIND=IB)           :: ioutside = 0
   INTEGER(KIND=IB)           :: ireject = 0, iaccept = 0
   INTEGER(KIND=IB)           :: ireject_bd = 0, iaccept_bd = 0
   INTEGER(KIND=IB)           :: i_bd                              !! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_zpert                           !! Z perturb track (0=nothing, 1=birth-death, 2=perturb 1 z)
   INTEGER(KIND=IB)           :: ipartfail 
   REAL(KIND=RP)              :: acc_ratio

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iskipais = 0                    !! Skip AIS if initial rjMCMC and resampling successful

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB):: MCMC_BI     !! # balance steps in traget region burn-in (exits if target reached)
   INTEGER(KIND=IB):: MCMC_STEPS  !! # balance steps in traget region (always does these)
   INTEGER(KIND=IB):: NTEMP       !! # temperature steps for initial AIS (if fail, increases by factor 3)
   INTEGER(KIND=IB):: NTEMP1        
   INTEGER(KIND=IB):: NBAL        !! # steps to balance at each AIS temp (if fail, increases by factor 1.3)
   INTEGER(KIND=IB):: NBAL1         
   REAL(KIND=IB)   :: npercnt     !! % particles that must be in logL target region
   INTEGER(KIND=IB):: NCOOLTRAJ   !! # cooling trajectories to try before giving up
   INTEGER(KIND=IB):: NAP         !! Misc parameters in sample (for bookeeping)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: beta             !! Inverse Temprerature
   REAL(KIND=RP)                         :: beta1gl          !! Global inverse T to define annealing schedule (start)
   REAL(KIND=RP)                         :: beta4gl          !! Global inverse T to define annealing schedule (end for 3 exp legs)

!!
!!  Structures for objects and data 
!!
   TYPE objstruc
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     !! Forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdpar   !! Std dev parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: z       !! Depth partition
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: h       !! Depth partition
      INTEGER(KIND=IB)                      :: k       !! Layer dimension
      INTEGER(KIND=IB)                      :: NFP     !! Number forward parameters
      INTEGER(KIND=IB)                      :: nfail   !! Counter to keep track of how often particle failed
      REAL(KIND=RP)                         :: logL    !! log likelihood
      REAL(KIND=RP)                         :: logP    !! log likelihood
      REAL(KIND=RP)                         :: logwt   !! log wt (likelihood ratio when moving from one ping to next)
      REAL(KIND=RP)                         :: logPr   !! log Prior probability ratio
      REAL(KIND=RP)                         :: lognorm !! Data covariance matrices
   END TYPE objstruc

   TYPE :: datastruc
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: Robs    !! Observed data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: angobs  !! Observed angles
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: Rrep    !! Replica data for trial models
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: res     !! Replica data for trial models
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)  :: sigma   !! Index for bad points/data gaps
                                                           !! meaningful for simulations)
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: lognorm !! Data lognorm
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)    :: logdet  !! Data log-determinant
      INTEGER(KIND=IB),ALLOCATABLE, DIMENSION(:):: NDPF    !! No. data per freq
      INTEGER(KIND=IB)                          :: NANG    !! No. data/angles (struc copy)
      INTEGER(KIND=IB)                          :: NBAND   !! No. freq bands (struc copy)
   END TYPE datastruc

   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: sdev        !! Standard devs
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: icount

!!
!!  Global variables for spherical reflection coeff computation
!!
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE    :: N,Nr
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE       :: drTh
   COMPLEX(KIND=RP),DIMENSION(:),ALLOCATABLE    :: diTh
   COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE  :: rTh,RPs
   COMPLEX(KIND=RP),DIMENSION(:,:,:),ALLOCATABLE:: btR

!!
!!  MPI global variables
!!
  INTEGER(KIND=IB)            :: rank,NTHREAD,ierr
  INTEGER(KIND=IB), PARAMETER :: src = 0_IB
  INTEGER                     :: to,from,tag,COMM
  INTEGER                     :: status(MPI_STATUS_SIZE)
  INTEGER(KIND=IB)            :: isize1,isize2,isize3

  REAL(KIND=RP)               :: tcmp,tstartcmp,tendcmp    !! Timing
  REAL(KIND=RP)               :: tstartsnd, tendsnd        !! Communication time
  REAL(KIND=RP)               :: tstart, tend              !! Overall time 
  REAL(KIND=RP)               :: tiping1, tiping2          !! Overall time 
  REAL(KIND=RP)               :: tstart2, tend2            !! Time for one forward model computation

  INTERFACE
    FUNCTION RANDPERM(num)
       USE data_type, ONLY : IB
       IMPLICIT NONE
       INTEGER(KIND=IB), INTENT(IN) :: num
       INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
    END FUNCTION RANDPERM
  END INTERFACE
  CONTAINS
  !!==============================================================================
  INTEGER FUNCTION NEWUNIT(unit)
  !!==============================================================================
    integer, intent(out), optional :: unit
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
        newunit=lun
        exit
      end if
    end do
    if (present(unit)) unit=newunit
  END FUNCTION NEWUNIT
  !!=======================================================================

  SUBROUTINE ALLOC_STRUC(dat,particles,particles_new,obj,objnew)
  !!
  !! Allocates all derived types. 
  !!
  !!=======================================================================
  USE MPI
  !USE RJMCMC_COM
  IMPLICIT NONE
  TYPE(objstruc),ALLOCATABLE,DIMENSION(:):: particles,particles_new
  TYPE(objstruc):: obj,objnew
  TYPE(datastruc),ALLOCATABLE,DIMENSION(:):: dat
  INTEGER(KIND=IB):: ip

  ALLOCATE(dat(NPING))
  DO ip = 1,NPING
    ALLOCATE(dat(ip)%Robs(NBAND,NANG))
    ALLOCATE(dat(ip)%angobs(NANG))
    ALLOCATE(dat(ip)%Rrep(NBAND,NANG))
    ALLOCATE(dat(ip)%res(NBAND,NANG))
    ALLOCATE(dat(ip)%sigma(NBAND,NANG))
    ALLOCATE(dat(ip)%lognorm(NBAND))
    ALLOCATE(dat(ip)%logdet(NBAND))
    ALLOCATE(dat(ip)%NDPF(NBAND))
    dat(ip)%NANG  = NANG  ! # data/angles (struc copy)
    dat(ip)%NBAND = NBAND ! # freq bands (struc copy)
  ENDDO

  ALLOCATE(particles(NPART),particles_new(NPART))
  DO ip = 1,NPART
    ALLOCATE(particles(ip)%par((NLMX*NPL)+(NPL-1)))
    ALLOCATE(particles(ip)%sdpar(NBAND))
    ALLOCATE(particles(ip)%z(NLMX))
    ALLOCATE(particles(ip)%h(NLMX))
    ALLOCATE(particles_new(ip)%par((NLMX*NPL)+(NPL-1)))
    ALLOCATE(particles_new(ip)%sdpar(NBAND))
    ALLOCATE(particles_new(ip)%z(NLMX))
    ALLOCATE(particles_new(ip)%h(NLMX))
  ENDDO
  ALLOCATE(obj%par((NLMX*NPL)+(NPL-1)))
  ALLOCATE(obj%sdpar(NBAND))
  ALLOCATE(obj%z(NLMX))
  ALLOCATE(obj%h(NLMX))
  ALLOCATE(objnew%par((NLMX*NPL)+(NPL-1)))
  ALLOCATE(objnew%sdpar(NBAND))
  ALLOCATE(objnew%z(NLMX))
  ALLOCATE(objnew%h(NLMX))
  RETURN
  END SUBROUTINE ALLOC_STRUC

END MODULE RJMCMC_COM
