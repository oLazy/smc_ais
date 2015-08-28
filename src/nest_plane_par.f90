!==============================================================================
!
!  Nested Sampling for plane wave reflection coefficient inversion
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, July 12 2008
!  jand@uvic.ca                       (250) 472 4342
!  http://web.uvic.ca~/jand/
!  Last change: July 12 2008
!
!  Based on "nested sampling algorithm" (John Skilling 2004, 2005, 2006)
!
!==============================================================================

MODULE DATA_TYPE
   IMPLICIT NONE
   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.0D0), SP=KIND(0.0)
   REAL(KIND=RP),   PARAMETER :: PI  = 3.141592653589793238462643383279502884197_RP
END MODULE DATA_TYPE

!==============================================================================
MODULE NEST_COM
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE

!
! General switches
!
   INTEGER(KIND=IB), PARAMETER :: ICOV    = 0            ! Data cov. switch
   INTEGER(KIND=IB), PARAMETER :: IPREJAC = 0            ! Data cov. switch
   INTEGER(KIND=IB), PARAMETER :: ISEEDST = 1            ! 1: Seed good starting model
   INTEGER(KIND=IB), PARAMETER :: NAVEF   = 12           ! # freq per band
   REAL(KIND=RP),    PARAMETER :: frbw    = 1._RP/20._RP ! Frac. bandwidth for freq ave.

!!
!!  SIM A
!!
!   INTEGER(KIND=IB), PARAMETER :: NANG = 131
!   INTEGER(KIND=IB), PARAMETER :: NFP  = 7 
!   REAL(KIND=RP),    PARAMETER :: cw   = 1511.0_RP
!   REAL(KIND=RP),    PARAMETER :: rw   = 1.029_RP
!   CHARACTER(len=64):: infile1       = 'sim_A_1_3_16_1lay.txt'
!   CHARACTER(len=64):: covfile       = 'sim_A_1_3_16_1lay_cov.txt'
!   CHARACTER(len=64):: mapfile       = 'sim_A_1_3_16_1lay_map.dat'
!   CHARACTER(len=64):: sdfile        = 'sim_A_1_3_16_exstd.txt'
!   CHARACTER(len=64):: samplefile    = 'sim_A_1_3_16_1lay_sample.txt'
!   CHARACTER(len=64):: jacfile       = 'sim_A_1_3_16_1lay_preJac.txt'
!   CHARACTER(len=64):: lincovfile    = 'sim_A_1_3_16_1lay_lincov.txt'
!   CHARACTER(len=64):: nonlincovfile = 'sim_A_1_3_16_1lay_nonlincov.txt'
!
!   INTEGER(KIND=IB), PARAMETER :: NBAND  = 8
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 300._RP, 400._RP,  504._RP,  635._RP,  800._RP, 1008._RP, 1270._RP, & 
!                                                 1600._RP /)

!!
!!  SIM B
!!
!   INTEGER(KIND=IB), PARAMETER :: NANG   = 131
!   INTEGER(KIND=IB), PARAMETER :: NFP    = 19
!   REAL(KIND=RP),    PARAMETER :: cw     = 1511.0_RP
!   REAL(KIND=RP),    PARAMETER :: rw     = 1.029_RP
!   CHARACTER(len=64) :: infile1    = 'sim_B_1_10_50_4lay.txt'
!   CHARACTER(LEN=64) :: logfile    = 'sim_B_1_10_50_4lay_FGS.log'
!   CHARACTER(len=64) :: covfile    = 'sim_B_1_10_50_4lay_cov.txt'
!   CHARACTER(len=64) :: mapfile    = 'sim_B_1_10_50_4lay_map.dat'
!   CHARACTER(len=64) :: sdfile     = 'sim_B_1_10_50_exstd.txt'
!   CHARACTER(len=64) :: samplefile = 'sim_B_1_10_50_4lay_sample.txt'
!
!   INTEGER(KIND=IB), PARAMETER :: NBAND  = 14
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 1000._RP, 1200._RP, 1400._RP, 1683._RP, 1889._RP, 2121._RP, &
!                                                  2381._RP, 2672._RP, 3000._RP, 3367._RP, 3779._RP, 4242._RP, &
!                                                  4762._RP, 5345._RP /)

!!
!!  Site 01
!!
!   INTEGER(KIND=IB), PARAMETER :: NANG    = 131
!   INTEGER(KIND=IB), PARAMETER :: NFP     = 11
!   REAL(KIND=RP),    PARAMETER :: cw     = 1511.18_RP
!   REAL(KIND=RP),    PARAMETER :: rw   = 1.029_RP
!   CHARACTER(len=64) :: infile1    = 'x_s01_1_3_20_2layb.txt'
!   CHARACTER(len=64) :: covfile    = 'x_s01_1_3_20_2lay_cov.txt'
!   CHARACTER(len=64) :: mapfile    = 'x_s01_1_3_20_2laybcov_map.dat'
!   CHARACTER(len=64) :: sdfile     = 'x_s01_1_3_20_exstdb.txt'
!   CHARACTER(len=64) :: samplefile = 'x_s01_1_3_20_2laybcov_sample.txt'
!
!   INTEGER(KIND=IB), PARAMETER :: NBAND  = 9
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 300._RP,  400._RP,  504._RP,  635._RP,  800._RP, 1008._RP, 1270._RP, &
!                                                 1600._RP, 2000._RP /)

!!
!!  Site 02
!!
   INTEGER(KIND=IB), PARAMETER :: NANG = 131
   INTEGER(KIND=IB), PARAMETER :: NFP  = 3
   REAL(KIND=RP),    PARAMETER :: cw   = 1510.78_RP
   REAL(KIND=RP),    PARAMETER :: rw   = 1.029_RP
   CHARACTER(len=64) :: infile1    = 'x_s02_1_3_16_7lay.txt'
   CHARACTER(LEN=64) :: logfile    = 'x_s02_1_3_16_7lay_FGS.log'
   CHARACTER(len=64) :: covfile    = 'x_s02_1_3_16_4lay_cov.txt'
   CHARACTER(len=64) :: mapfile    = 'x_s02_1_3_16_7laycov_map.dat'
   CHARACTER(len=64) :: sdfile     = 'x_s02_1_3_16_exstd.txt'
   CHARACTER(len=64) :: samplefile = 'x_s02_1_3_16_7lay_sample.txt'
   CHARACTER(len=64):: jacfile       = 'sim_A_1_3_16_1lay_preJac.txt'
   CHARACTER(len=64):: lincovfile    = 'sim_A_1_3_16_1lay_lincov.txt'
   CHARACTER(len=64):: nonlincovfile = 'sim_A_1_3_16_1lay_nonlincov.txt'

   INTEGER(KIND=IB), PARAMETER :: NBAND  = 8
   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 300._RP, 400._RP,  504._RP,  635._RP,  800._RP, 1008._RP, 1270._RP, & 
                                                 1600._RP /)

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP), DIMENSION(NFP):: minlim   = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxlim   = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxpert  = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: w        = 0._RP 
   REAL(KIND=RP), DIMENSION(NFP):: minlimp  = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxlimp  = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxpertp = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: wp       = 0._RP 
   REAL(KIND=RP), DIMENSION(NFP):: wscale   = 10._RP 

!!
!!  Nested sampling specific parameters
!!
   INTEGER(KIND=IB)           :: N     = NFP*NFP*200_IB ! # objects in box
!   INTEGER(KIND=IB)           :: N     = NFP*NFP*40_IB ! # objects in box
!   INTEGER(KIND=IB)           :: N     = NFP*NFP*10_IB ! # objects in box
!   INTEGER(KIND=IB)           :: N     = 80_IB         ! # objects in box
   INTEGER(KIND=IB),PARAMETER :: JMAX = 20_IB           ! Max # MCMC iterations
   INTEGER(KIND=IB)           :: NCPU                  ! # objects per CPU (calculated from N)

!!
!!  Convergence parameters
!!
   REAL(KIND=RP),PARAMETER:: endfac   = 1.5_RP  ! Termination factor
   REAL(KIND=RP),PARAMETER:: logtol   = 1E-2_RP ! Termination tolerance (log space)
   REAL(KIND=RP),PARAMETER:: corconv1 = 0.8_RP  ! Start using non-linear estimate in burn-in
   REAL(KIND=RP),PARAMETER:: corconv2 = 0.5_RP  ! Finish burn-in
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch nested sampling
   INTEGER(KIND=IB)       :: ilast    = 0       ! Last call to save in nested sampling
!!
!! Bookkeeping parameters
!!
   INTEGER(KIND=IB),PARAMETER:: IMAX   = 1E6_IB   ! # iterations (# nested boxes)
   INTEGER(KIND=IB),PARAMETER:: NKEEP  = 100
   INTEGER(KIND=IB),PARAMETER:: NPRINT = 100
   INTEGER(KIND=IB),PARAMETER:: NCONV  = 100
   INTEGER(KIND=IB),PARAMETER:: NCOPY  = 3
   INTEGER(KIND=IB),PARAMETER:: NDM    = 60
   INTEGER(KIND=IB),PARAMETER:: NAP    = 5

!!
!!  Structures for objects and data 
!!
   TYPE :: objstruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NFP):: par   = 0._RP   ! Forward parameters
      REAL(KIND=RP),DIMENSION(NFP):: prior = 0._RP   ! Prior parameters drawn from [0,1]
      REAL(KIND=RP)               :: logwt           ! log weights
      REAL(KIND=RP)               :: logL            ! log likelihood
   END TYPE objstruc

   TYPE :: datastruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Robs   = 0._RP ! Observed data
      REAL(KIND=RP),DIMENSION(NANG)           :: angobs = 0._RP ! Observed angles
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rrep   = 0._RP ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: res    = 0._RP ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rex    = 0._RP ! Index for bad points/data gaps
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cdi    = 0._RP ! Inverse data covariance matrices
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cd     = 0._RP ! Data covariance matrices
      REAL(KIND=RP),DIMENSION(NFP)            :: start  = 0._RP ! Good starting object
      REAL(KIND=RP),DIMENSION(NFP)            :: true   = 0._RP ! True object (only
                                                                ! meaningful for simulations)
      REAL(KIND=RP),DIMENSION(NFP,NFP)        :: Cov0   = 0._RP ! Linear model covariance
      REAL(KIND=RP),DIMENSION(NFP,NFP)        :: Cov    = 0._RP ! Non-linear model covariance
      REAL(KIND=RP),DIMENSION(NFP,NFP)        :: VV     = 0._RP ! Non-linear rotation matrix
      REAL(KIND=RP),DIMENSION(NFP)            :: sdevm  = 0._RP ! Non-linear rotation matrix
      REAL(KIND=RP)                           :: lognorm=0._RP  ! Data covariance matrices
      INTEGER(KIND=IB), DIMENSION(NBAND)      :: NDPF           ! # data per freq
      INTEGER(KIND=IB)                        :: NANG   = NANG  ! # data/angles (struc copy)
      INTEGER(KIND=IB)                        :: NBAND  = NBAND ! # freq bands (struc copy)
      INTEGER(KIND=IB),DIMENSION(NBAND)       :: NDAT   = 0     ! Observed angles
   END TYPE datastruc
   REAL(KIND=RP),DIMENSION(NBAND):: sdev = 0._RP                ! Standard devs
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: icount
   INTEGER(KIND=IB)                         :: idxdim
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: idxr,idxc
   REAL(KIND=RP),   DIMENSION(:),ALLOCATABLE:: dc,dr

   TYPE :: covstruc
      SEQUENCE
      INTEGER(KIND=IB)            :: jsamp   = 0_IB
      INTEGER(KIND=IB)            :: ksamp   = 0_IB
      INTEGER(KIND=IB)            :: iburnin = 0_IB
      REAL(KIND=RP)               :: dcor    = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: msum1   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: msum2   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: msum3   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcsum1  = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcsum2  = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcsum3  = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: mave1   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: mave2   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP):: mave3   = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcross1 = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcross2 = 0._RP   
      REAL(KIND=RP),DIMENSION(NFP,NFP):: mcross3 = 0._RP   
   END TYPE covstruc

!   TYPE :: smpstruc
!      SEQUENCE
!      REAL(KIND=RP)               :: logL            ! log likelihood (worst)
!      REAL(KIND=RP)               :: logLb           ! log likelihood (best)
!      REAL(KIND=RP)               :: logwidth        ! logwidth
!      REAL(KIND=RP)               :: logwt           ! log weights
!      REAL(KIND=RP)               :: H               ! information
!      REAL(KIND=RP)               :: logZ            ! log evidence
!      REAL(KIND=RP),DIMENSION(NFP):: par   = 0._RP   ! Forward parameters
!   END TYPE smpstruc

!!
!!  Global variables
!!
   REAL(KIND=RP)               :: logZG, HG

!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ncount,ierr
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,tag,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   REAL(KIND=RP)               :: tstart, tend 
   INTEGER(KIND=IB)            :: CHUNK_BG,CHUNK_SM,REST,CHUNK_SL
   INTEGER(KIND=IB)            :: CHUNK_BEG,CHUNK_END

END MODULE NEST_COM
  
!=======================================================================

PROGRAM  NESTED_PLANE

!=======================================================================
USE MPI
USE NEST_COM
USE NR
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ifreq,inest,ikeep,iprint1,iprint2,icopy2,iconvtest
!INTEGER(KIND=IB), DIMENSION(2), PARAMETER :: iseed = 5678

TYPE (objstruc), DIMENSION(:),ALLOCATABLE  :: obj     ! Objects in likelihood box
TYPE (datastruc)                           :: dat     ! Data

REAL(KIND=RP),   DIMENSION(:,:),ALLOCATABLE:: sample
REAL(KIND=RP)                              :: ran_uni ! Uniform random number

INTEGER(KIND=IB)  :: iworst,ibest,icopy
REAL(KIND=RP)     :: logLstar      ! Hard likelihood constraint
REAL(KIND=RP)     :: logwidth      ! Log width of integration interval
REAL(KIND=RP)     :: logwidth2     ! Log width of inwardstep interval
REAL(KIND=RP)     :: logZ,logZnew  ! Running sum for log evidence
REAL(KIND=RP)     :: H             ! Running sum for neg. entropy/information
REAL(KIND=RP)     :: LOGPLUS       ! Function: Addition carried out in log space

!!
!!  Specific to linear model covariance matrix estimate
!!
REAL(KIND=RP),DIMENSION(NANG,NANG)      :: VVV  = 0._RP
REAL(KIND=RP),DIMENSION(NANG,NANG)      :: WWWi = 0._RP
REAL(KIND=RP),DIMENSION(NANG)           :: WWW  = 0._RP
REAL(KIND=RP)                           :: Wmax = 0._RP
REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: CdC  = 0._RP
REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: CdCi = 0._RP
INTEGER(KIND=IB),DIMENSION(NANG)        :: indx = 0_IB
REAL(KIND=SP)                           :: d    = 0._SP
REAL(KIND=RP)                           :: logdet = 0._RP

!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

INTEGER(KIND=IB), DIMENSION(2), PARAMETER :: iseed1 = 6788

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )

NCPU = CEILING(REAL(N,RP)/REAL(NTHREAD,RP))
N = NCPU*NTHREAD

ALLOCATE( obj(NCPU),sample(NKEEP+NCPU,NAP+NFP))
ALLOCATE( icount(NTHREAD) )
icount = 0

CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!DO i = 1,NTHREAD
!   IF(i-1 .EQ. rank)THEN
!      WRITE(6,*) 'Process ', rank, ' of ', NTHREAD, ' is alive'
!      CALL FLUSH(6)
!   ENDIF
!   CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!ENDDO

tstart = MPI_WTIME()

!!------------------------------------------------------------------------
!!
!!  Print nested sampling parameters to screen for logging
!!
IF(rank == src)THEN
   WRITE(6,*) 'Max # iterations:',IMAX
   WRITE(6,*) 'Max # MCMC steps:',JMAX
   WRITE(6,*) '# objects in box:',N
   WRITE(6,*) 'Convergence when inest > N*H*',endfac

   WRITE(6,*) infile1
   WRITE(6,*) covfile
   WRITE(6,*) sdfile
   WRITE(6,*) samplefile
   WRITE(6,*) ''
ENDIF

!!
!! Find indexes for velocity-density coupling
!!
idxdim = ((NFP-3)/4)+1
ALLOCATE(idxc(idxdim),idxr(idxdim),dc(idxdim-1),dr(idxdim-1))
idxc(1) = 2
idxr(1) = 3
IF(idxdim > 2)THEN
   DO i = 2,idxdim-1
      idxc(i) = idxc(i-1) + 4
      idxr(i) = idxr(i-1) + 4
   ENDDO
ENDIF
idxc(idxdim) = NFP-2
idxr(idxdim) = NFP-1

IF(rank == src) WRITE(6,*) 'idxc',idxc
IF(rank == src) WRITE(6,*) 'idxr',idxr

!------------------------------------------------------------------------
!  Read in prior info
!------------------------------------------------------------------------
IF(rank == src)WRITE(6,*) 'Loading data...'
OPEN(UNIT=20,FILE=infile1,FORM='formatted',STATUS='OLD',ACTION='READ')

READ(20,*) (minlim(j),j=1,NFP)
READ(20,*) (maxlim(j),j=1,NFP)
READ(20,*) (dat%true(j),j=1,NFP)

!!
!! Interval size maxpert and initial slice sampling step size w
!!
maxpert = maxlim - minlim
w = (maxlim-minlim)/wscale

!!
!!  Set rotated bounds to normal bounds initially and rotation matrix to 
!!  identity matrix. After Cov is calculated, this gets updated.
!!
maxlimp  = maxlim
minlimp  = minlim
maxpertp = maxpertp
wp = w
DO i = 1,NFP
   dat%VV(i,i) = 1._RP
ENDDO


!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
DO ifreq = 1,NBAND
    READ(20,*) (dat%Robs(ifreq,j),j=1,NANG)
ENDDO

READ(20,*) (dat%angobs(j),j=1,NANG)

DO ifreq = 1,NBAND
    READ(20,*) (dat%Rex(ifreq,j),j=1,NANG)
ENDDO
CLOSE(20)

IF(rank == src)THEN
   WRITE(6,*) 'Number of angles:         ',NANG
   WRITE(6,*) 'Read data.'
   CALL FLUSH(6)
ENDIF

IF(ICOV .EQ. 1)THEN
   IF(rank == src) WRITE(*,*) 'Using data covariance matrix estimate.'
   OPEN(UNIT=20,FILE=covfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   DO ifreq = 1,NBAND
      DO i = 1,NANG
         READ(20,*) (dat%Cdi(i,j,ifreq),j=1,NANG)
      ENDDO
   ENDDO
   DO ifreq = 1,NBAND
      DO i = 1,NANG
         READ(20,*) (dat%Cd(i,j,ifreq),j=1,NANG)
      ENDDO
   ENDDO
   CLOSE(20)
ELSE
   IF(rank == src) WRITE(*,*) 'Using experimental error (one stdev per band) estimate.'
   OPEN(UNIT=20,FILE=sdfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   READ(20,*) (sdev(ifreq),ifreq=1,NBAND)
   CLOSE(20)
   DO ifreq = 1,NBAND
      DO i = 1,NANG
         dat%Cdi(i,i,ifreq) = 1._RP/sdev(ifreq)/sdev(ifreq)
         dat%Cd(i,i,ifreq)  = sdev(ifreq)*sdev(ifreq)
      ENDDO
   ENDDO
ENDIF
DO ifreq = 1,NBAND

   CALL CHOLESKY(dat%Cd(:,:,ifreq),NANG,NANG,CdC(:,:,ifreq))

   CdCi(:,:,ifreq) = REAL(CdC(:,:,ifreq),KIND=SP)
   CALL SVDCMP(CdCi(:,:,ifreq),WWW,VVV)
   Wmax = MAXVAL(WWW,1)
   WWWi = 0.  
   DO i=1,NANG
       IF (WWW(i)/Wmax .GT. 1.e-7) THEN
             WWWi(i,i) = 1./WWW(i)
       ELSE
             WRITE(6,*) 'Warning: Small singular value'
             WRITE(6,*) WWW/Wmax
       ENDIF
   ENDDO
   CdCi(:,:,ifreq) = MATMUL(MATMUL(VVV,WWWi),TRANSPOSE(VVV))
ENDDO

!!
!!  Compute log(norm) for Gaussian data errors
!!  (computed log determinant by summing log of main diag of LU dcmp)
!!
dat%lognorm = 0._RP
DO ifreq = 1,NBAND
   VVV = REAL(dat%Cd(:,:,ifreq),KIND=SP)
   CALL LUDCMP(VVV,indx,d)
   logdet = 0._RP
   DO i = 1,NANG
      logdet = logdet + LOG(REAL(ABS(VVV(i,i)),KIND=RP))
   ENDDO
   dat%NDAT(ifreq) = COUNT(dat%Rex(ifreq,:)==1)
   dat%lognorm = dat%lognorm + (-0.5_RP*REAL(dat%NDAT(ifreq),KIND=RP)*LOG(2._RP*PI)- 0.5_RP*logdet)
ENDDO
IF(rank == src)WRITE(*,*) 'log(|C_d|) = ',logdet
IF(rank == src)WRITE(*,*) 'lognorm = ',dat%lognorm


!!
!!  Use good starting object (if expensive optimization was carried out)
!!
IF(ISEEDST == 1)THEN
   IF(rank == src) WRITE(6,*)'Reading map starting model:',mapfile
   OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   READ(20,*) (dat%start(j),j=1,NFP)
   CLOSE(20)

!!
!!  Constrain good starting object to within bounds 
!!
   DO i=1,NFP
     IF (dat%start(i) < minlim(i)) THEN
        IF(rank == src) WRITE(6,*) 'Starting Model Outside Bounds',i,dat%start(i),minlim(i),maxlim(i)
        dat%start(i) = minlim(i) + ABS(dat%start(i)-minlim(i))
     ENDIF
     IF (dat%start(i) > maxlim(i)) THEN
        IF(rank == src) WRITE(6,*) 'Starting Model Outside Bounds',i,dat%start(i),minlim(i),maxlim(i)
        dat%start(i) = maxlim(i) - ABS(dat%start(i)-maxlim(i))
     ENDIF
   ENDDO
ENDIF

!!
!!  Initialize variables for running sums:
!!
logZ  = -1.0E308_RP      ! Initializes Z to zero
H     = 0._RP            ! Initializes H to zero

!!
!!  Ensure unique random seed for each CPU
!!
CALL RANDOM_SEED
CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE( iseed(iseedsize), rseeds(iseedsize,NTHREAD), iseeds(iseedsize,NTHREAD) )
IF(rank == src)THEN
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1000000._RP)
ENDIF
CALL MPI_BCAST( iseeds(1,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
CALL MPI_BCAST( iseeds(2,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
iseed = iseeds(:,rank+1)
CALL RANDOM_SEED(PUT=iseed)
WRITE(*,*) 'Starting seed:',rank,iseed
CALL FLUSH(6)
CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

!!
!! Compute linear covariance estimate
!!
CALL LINROT(dat,CdCi)

obj(1)%par   = dat%start
obj(1)%prior = (obj(1)%par - minlim) / maxpert
CALL LOGLHOOD(obj(1),dat)
icount(rank+1) = icount(rank+1) + 1
!!
!! Compute non-linear covariance estimate by sampling (burn-in)
!!
CALL NONLINROT(obj(1),dat)

!------------------------------------------------------------------------
!
!           ************Nested Sampling************
!
! -----------------------------------------------------------------------
IF(rank == src)THEN
   WRITE(6,*) 'Starting nested sampling...'
   !!
   !!  File to save posterior sample
   !!
   OPEN(UNIT=40,FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
ENDIF

!!
!!  Draw objects (obj) randomly from prior
!!
DO i=1,NCPU

   CALL PRIOR(obj(i),dat)

ENDDO

!!
!!  Seed likely starting object
!!
IF(ISEEDST == 1)THEN
   IF(rank == src)THEN
      obj(1)%par   = dat%start
      obj(1)%prior = (obj(1)%par - minlim) / maxpert
      CALL LOGLHOOD(obj(1),dat)
      icount(rank+1) = icount(rank+1) + 1
   ENDIF
ENDIF

logwidth  = LOG(1.0_RP - EXP(-1.0_RP / N))    ! Outermost interval of prior mass
logwidth2 = LOG(1.0_RP - EXP(-1.0_RP / NCPU)) ! Outermost interval of prior mass

ikeep     = 1
iprint1   = 0
iprint2   = 0
icopy2    = 0
iconvtest = 1

!!
!!  Start nested sampling loop 
!!
DO inest = 1,IMAX

   CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
   iworst = 1
   ibest  = 1
   DO i=2,NCPU
      IF(obj(i)%logL < obj(iworst)%logL) iworst = i
      IF(obj(i)%logL > obj(ibest)%logL)  ibest  = i
   ENDDO
   obj(iworst)%logwt = logwidth2 + obj(iworst)%logL - LOG(REAL(NTHREAD,RP))

   !!
   !! Update evidence Z and information H
   !!
   logZnew = LOGPLUS(logZ, obj(iworst)%logwt)
   H = EXP(obj(iworst)%logwt - logZnew) * obj(iworst)%logL &
     + EXP(logZ - logZnew) * (H + logZ) - logZnew
   logZ = logZnew

   !!
   !! Update global evidence and test for convergence
   !!
   IF(inest >= (iconvtest*NCONV))THEN
      CALL CONVERGENCE(logZ,H,obj(ibest),inest)
      iconvtest = iconvtest + 1
   ENDIF

   !!
   !! Save sample
   !!
   sample(ikeep,:) = (/ obj(iworst)%logL,logwidth, obj(iworst)%logwt, H, logZG, obj(iworst)%par /)
   IF(ikeep == NKEEP .OR. ilast == 1)THEN
      CALL SAVESAMPLE(sample,ikeep)
      ikeep  = 0
      sample = 0._RP
   ENDIF
   ikeep  = ikeep + 1

   !!
   !! Write to screen
   !!
   IF(inest >= (iprint1*20*NPRINT)+1)THEN
      IF (rank == src) THEN
         tend = MPI_WTIME()
         WRITE(*,*) 'Time: ',tend-tstart,'s on',NTHREAD,'CPUs'
         WRITE(*,*) ''
         WRITE(6,206) ' inest,        N*H*fac,          logZ,              H,     logL(worst),     logL(best)'
         WRITE(6,206) '--------------------------------------------------------------------------------------'
      ENDIF
      iprint1 = iprint1 + 1 
   ENDIF
   IF(inest >= (iprint2*NPRINT))THEN
      IF (rank == src) THEN
         WRITE(6,205) inest,N*H*endfac,logZG,H,obj(iworst)%logL,obj(ibest)%logL
      ENDIF
      iprint2 = iprint2 + 1 
   ENDIF

   !!
   !! Evolve worst object into better one
   !!
   CALL RANDOM_NUMBER(ran_uni)
   logLstar  = obj(iworst)%logL
   IF(inest >= icopy2*NCOPY)THEN
     icopy = INT(ran_uni*NCPU)+1
     obj(iworst) = obj(icopy)
     icopy2 = icopy2 + 1
   ENDIF

   DO j = 1,JMAX
      CALL EXPLORE_SLICE_ROT(obj(iworst),dat,logLstar)
   ENDDO     !  MCMC iteration loop to JMAX

   !
   ! Shrink interval for prior mass
   !
   logwidth  = logwidth  - 1.0_RP/REAL(N,RP)
   logwidth2 = logwidth2 - 1.0_RP/REAL(NCPU,RP)   

   !!
   !! Convergence
   !!
   IF(iconv == 1)EXIT 

ENDDO
ikeep = ikeep - 1
!!
!!  End nested sampling loop 
!!

!!
!! Final correction:
!!
logwidth  = -REAL(inest,RP)/REAL(N,RP) - LOG(REAL(N,RP))
logwidth2 = -REAL(inest,RP)/REAL(NCPU,RP) - LOG(REAL(NCPU,RP))
DO i = 1,NCPU
   obj(i)%logwt = logwidth2 + obj(i)%logL - LOG(REAL(NTHREAD,RP))
   logZnew = LOGPLUS(logZ, obj(i)%logwt)
   H = EXP(obj(i)%logwt - logZnew) * obj(i)%logL &
     + EXP(logZ - logZnew) * (H + logZ) - logZnew
   logZ = logZnew

ENDDO

!!
!! Add the live models to the count
!!
ilast = 1
CALL CONVERGENCE(logZ,H,obj(ibest),inest)

!!
!! Additional last NCPU samples
!!
DO i = 1,NCPU
   sample(ikeep+i,:) = (/ obj(i)%logL,logwidth, obj(i)%logwt, H, logZG, obj(i)%par /)
ENDDO
ikeep = ikeep + NCPU
CALL SAVESAMPLE(sample,ikeep)
CLOSE(40)

!!
!! Print results:
!!
IF(rank == src)THEN
   WRITE(*,*) 'FINAL RESULT:'
   WRITE(*,202) logZG, SQRT(H/N),inest
   tend = MPI_WTIME()
   WRITE(6,*)'Number of likelihood evals: ',SUM(icount),' on ',NTHREAD,' CPUs'
   WRITE(6,*)'Total time: ',tend-tstart,'s per CPU'
ENDIF

201 FORMAT(2x, 'Thread ',I3,' results:')
202 FORMAT(2x, 'Evidence ln(Z) = ',F12.6,'+-', F12.6,' No. iterates per CPU = ',I6)
203 FORMAT(2x, 'Information H  = ',F12.6,'nats = ', F12.6,' bits')
205 FORMAT(I6,5(f16.4))
206 FORMAT(A86)
207 FORMAT(50F16.4)

DEALLOCATE( iseed,rseeds,iseeds )
CALL MPI_FINALIZE( ierr ) 
END PROGRAM NESTED_PLANE

!=======================================================================

SUBROUTINE SAVESAMPLE(sample,ikeep)
!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE MPI
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,ikeep

REAL(KIND=RP),DIMENSION(NKEEP+NCPU,NAP+NFP):: sample
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sampleM
REAL(KIND=RP),DIMENSION(:),  ALLOCATABLE   :: buf2
INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE  :: ikeepM

!!
!! Exchanging number of samples
!!
IF(rank == src)THEN
   ALLOCATE( ikeepM(NTHREAD) )
   ikeepM = 0_IB

   ikeepM(1) = ikeep
   DO i = 1,NTHREAD-1
      from  = i
      tag   = MPI_ANY_TAG
      CALL MPI_RECV(ikeepM(i+1), 1, MPI_INTEGER, from, &
                  tag, MPI_COMM_WORLD, status, ierr )
   ENDDO
ELSE
   CALL MPI_SSEND( ikeep, 1, MPI_INTEGER, src, & 
                   rank, MPI_COMM_WORLD, ierr )
ENDIF

!!
!!  Sending samples to master
!!
IF(rank == src)THEN
   ALLOCATE( sampleM(NTHREAD*(ikeep),NAP+NFP) )
   sampleM(1:ikeep,:) = sample
   ikeep = ikeepM(1)
   DO i = 1,NTHREAD-1
      from  = i
      tag   = MPI_ANY_TAG
      ncount = (ikeepM(i+1))*(NFP+NAP)
      ALLOCATE( buf2(ncount) )
      CALL MPI_RECV(buf2, ncount, MPI_DOUBLE_PRECISION, from,&
                    tag, MPI_COMM_WORLD, status, ierr )
      sampleM(ikeep+1:ikeep+ikeepM(i+1),:) = RESHAPE(buf2,(/ ikeepM(i+1),NFP+NAP /))
      ikeep = ikeep + ikeepM(i+1)
      DEALLOCATE( buf2 )
   ENDDO
   ikeep = ikeepM(1)
   !!
   !! Master writes sample to file
   !!
   DO i=1,NTHREAD*(ikeep)
      WRITE(40,207) sampleM(i,:)
   ENDDO
   DEALLOCATE( sampleM,ikeepM )
ELSE
   ncount = ikeep*(NFP+NAP)
   ALLOCATE( buf2(ncount) )
   buf2 = PACK(sample(1:ikeep,:),sample(1:ikeep,:) /= SQRT(-1._RP))
   CALL MPI_SSEND( buf2, ncount, MPI_DOUBLE_PRECISION, src, &
                   rank, MPI_COMM_WORLD, ierr )
   DEALLOCATE( buf2 )
ENDIF !  MPI

207   FORMAT(50F16.4)
RETURN
END SUBROUTINE SAVESAMPLE

!=======================================================================

SUBROUTINE CONVERGENCE(logZ,H,obj,inest)
!=======================================================================
USE MPI
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB)                      :: i,inest
REAL(KIND=RP)                         :: LOGPLUS  ! Function: Addition in log space
TYPE (objstruc)                       :: obj      ! Best object

INTEGER(KIND=IB)                      :: ibestG
REAL(KIND=RP)                         :: logZ,H
REAL(KIND=RP),DIMENSION(:),ALLOCATABLE:: logZM,HM,logLM,logwtM

ALLOCATE( logZM(NTHREAD),HM(NTHREAD),logLM(NTHREAD),logwtM(NTHREAD) )
logZM  = 0._RP
logLM  = 0._RP
logwtM = 0._RP
HM     = 0._RP

!!
!! Exchanging logZ and H
!!
IF(rank == src)THEN
   logZM(1) = logZ
   logLM(1) = obj%logL
   logwtM(1) = obj%logwt
   HM(1) = H
   DO i = 1,NTHREAD-1
      from  = i
      tag   = MPI_ANY_TAG
      CALL MPI_RECV(logZM(i+1), 1, MPI_DOUBLE_PRECISION, from, &
                  tag, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV(logLM(i+1), 1, MPI_DOUBLE_PRECISION, from, &
                  tag, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV(logwtM(i+1), 1, MPI_DOUBLE_PRECISION, from, &
                  tag, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV(HM(i+1), 1, MPI_DOUBLE_PRECISION, from, &
                  tag, MPI_COMM_WORLD, status, ierr )
   ENDDO
ELSE
   CALL MPI_SSEND( logZ, 1, MPI_DOUBLE_PRECISION, src, & 
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( obj%logL, 1, MPI_DOUBLE_PRECISION, src, & 
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( obj%logwt, 1, MPI_DOUBLE_PRECISION, src, & 
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_SSEND( H, 1, MPI_DOUBLE_PRECISION, src, & 
                   rank, MPI_COMM_WORLD, ierr )
ENDIF
!!
!! Find global best model to test convergence
!!
ibestG = MAXLOC(logLM,1)

!!
!! Compute global evidence:
!!
IF(rank == src)THEN
   logZG = -1.0E308_RP       ! Initializes global Z to zero
   DO i = 1,NTHREAD
       logZG = LOGPLUS(logZG, logZM(i))
   ENDDO
ENDIF
CALL MPI_BCAST(logZG,1,MPI_DOUBLE_PRECISION,src,MPI_COMM_WORLD,ierr)

!!
!! Convergence
!!
IF(iconv == 0)THEN

   !! Max amount evidence can gain is small:
!   IF(ABS(LOGPLUS(logwtM(ibestG),logZG)-logZG) < logtol) iconv = 1

   !! Both seems much more stable in parallel code
   IF(ABS(LOGPLUS(logwtM(ibestG),logZG)-logZG) < logtol .AND. &
      inest > CEILING(endfac*NCPU*H)) THEN 
      iconv = 1
      IF(rank == src)WRITE(*,*) '**** CONVERGED ****'
   ENDIF
   CALL MPI_BCAST(iconv,1,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)

ENDIF

201 FORMAT(2x, 'Thread ',I3,' results:')
202 FORMAT(2x, 'Evidence ln(Z) = ',F12.6,'+-', F12.6)
203 FORMAT(2x, 'Information H  = ',F12.6,'nats = ', F12.6,' bits')
DEALLOCATE( logZM,HM,logLM,logwtM )
RETURN
END SUBROUTINE CONVERGENCE

!=======================================================================

SUBROUTINE LOGLHOOD(obj,dat)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifreq,iang
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP), DIMENSION(NAVEF,NANG):: Rpltry
REAL(KIND=RP), DIMENSION(NBAND)     :: Etmp
REAL(KIND=RP), DIMENSION(NAVEF)     :: fr
REAL(KIND=RP)                       :: flo,fhi,fstep

!!
!!  Compute plane wave refl. coeff. (band averaged)
!!
fr = 0._RP
DO ifreq = 1,NBAND
   flo = bands(ifreq) - bands(ifreq)*frbw
   fhi = bands(ifreq) + bands(ifreq)*frbw
   fstep = (fhi-flo)/(NAVEF-1)
   fr = flo + ((/ 1:NAVEF /)-1) * fstep
   CALL REF_NLAY3(dat%angobs,obj%par,fr,Rpltry,cw,rw,NAVEF,NANG,NFP)

   DO iang = 1,NANG
      dat%Rrep(ifreq,iang) = SUM(Rpltry(:,iang)*EXP(-(fr-bands(ifreq))**2._RP/ &
                             (frbw*bands(ifreq))**2._RP)*fstep,1)/ &
                             SUM(EXP(-(fr-bands(ifreq))**2._RP/ &
                             (frbw*bands(ifreq))**2._RP)*fstep,1);
   ENDDO
ENDDO

!!
!!  Compute log likelihood
!!
dat%res = (dat%Robs-dat%Rrep)*dat%Rex
DO ifreq = 1,NBAND
   Etmp(ifreq) = DOT_PRODUCT( MATMUL(dat%res(ifreq,:),&
                 dat%Cdi(:,:,ifreq)),dat%res(ifreq,:))/2._RP
ENDDO
obj%logL = dat%lognorm - SUM(Etmp)

RETURN
END SUBROUTINE LOGLHOOD

!=======================================================================

SUBROUTINE PRIOR(obj,dat)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=RP) :: j,ichck_cr
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP)    :: ran_uni(NFP)

DO
   CALL RANDOM_NUMBER(ran_uni)
   obj%par = minlim + maxpert * ran_uni
   CALL COUPLE_CR(obj,dat,ichck_cr)
   IF(ichck_cr == 1)CYCLE
   CALL LOGLHOOD(obj,dat)
   icount(rank+1) = icount(rank+1) + 1
   EXIT
ENDDO

128 FORMAT(60F16.6)
RETURN
END SUBROUTINE PRIOR

!=======================================================================

SUBROUTINE EXPLORE_SLICE_ROT(obj,dat,logLstar)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=IB)  :: accept,reject,i,j,k,jcount,kcount
INTEGER(KIND=IB)  :: ichck_cr = 0
TYPE (objstruc) :: obj,try,tryp
TYPE (datastruc):: dat
REAL(KIND=RP)  :: logLstar,logy
REAL(KIND=RP)  :: ran_uni,u
REAL(KIND=RP)  :: L,R,x0,x1,gx0


try = obj

DO i = 1,NFP

   tryp     = obj
   gx0      = obj%logL
   tryp%par = MATMUL(TRANSPOSE(dat%VV),obj%par)
   x0       = tryp%par(i)

   !!
   !!  Slice level = logLstar
   !!
   logy = logLstar
   !!
   !!  Fin initial sampling interval
   !!
   CALL RANDOM_NUMBER(ran_uni)
   u = ran_uni * wp(i)
   L = x0 - u
   R = x0 + (wp(i)-u)
   
   !!
   !!  Expand interval until outside slice
   !!
   jcount = 0
   DO
      tryp%par(i) = L
      try%par = MATMUL(dat%VV,tryp%par)
      IF(MINVAL(try%par - minlim) < 0._RP)EXIT
      IF(MINVAL(maxlim - try%par) < 0._RP)EXIT

      CALL LOGLHOOD(try,dat)
      icount(rank+1) = icount(rank+1) + 1
      jcount = jcount + 1
      IF(try%logL <= logy)EXIT
      L = L - wp(i)
      IF(jcount > 5000)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop1'
   ENDDO
   jcount = 0
   DO
      tryp%par(i) = R
      try%par = MATMUL(dat%VV,tryp%par)
      IF(MINVAL(try%par - minlim) < 0._RP)EXIT
      IF(MINVAL(maxlim - try%par) < 0._RP)EXIT

      CALL LOGLHOOD(try,dat)
      icount(rank+1) = icount(rank+1) + 1
      jcount = jcount + 1
      IF(try%logL <= logy)EXIT
      R = R + wp(i)
      IF(jcount > 5000)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop2'
   ENDDO

   !!
   !!  Sample from interval, shrinking it on each rejection
   !!
   jcount = 0
   kcount = 0
   DO
      CALL RANDOM_NUMBER(ran_uni)
      x1 = L+(R-L)*ran_uni
      tryp%par(i) = x1
      try%par = MATMUL(dat%VV,tryp%par)
      kcount = kcount + 1
      IF(kcount > 5000)THEN
         WRITE(*,*) 'WARNING: may be stuck in infinite cr-constrain loop',kcount
         WRITE(*,*) dc
         WRITE(*,*) dr
         WRITE(*,*) ''
         kcount = 0
      ENDIF
      !!
      !!  Make sure all models are from within phys. space  bounds.
      !!
      IF(MINVAL(try%par - minlim) < 0._RP)THEN
         IF(x1 > x0)THEN
            R = x1
         ELSE
            L = x1
         ENDIF
         CYCLE
      ENDIF
      IF(MINVAL(maxlim - try%par) < 0._RP)THEN
         IF(x1 > x0)THEN
            R = x1
         ELSE
            L = x1
         ENDIF
         CYCLE
      ENDIF

      CALL COUPLE_CR(try,dat,ichck_cr)
      IF(ichck_cr == 1)THEN
!         IF(x1 > x0)THEN
!            R = x1
!         ELSE
!            L = x1
!         ENDIF
!         CYCLE
      ENDIF

      CALL LOGLHOOD(try,dat)
      icount(rank+1) = icount(rank+1) + 1
      jcount = jcount + 1
      IF(try%logL >= logy)EXIT
      IF(x1 > x0)THEN
         R = x1
      ELSE
         L = x1
      ENDIF
      IF(jcount > 5000)WRITE(*,*) 'WARNING: may be stuck in infinite slice loop3'
   ENDDO
   obj = try

ENDDO  !  parameter loop

128 FORMAT(60F16.6)
RETURN
END SUBROUTINE EXPLORE_SLICE_ROT


!=======================================================================

SUBROUTINE COUPLE_CR(obj,dat,ichck_cr)
!=======================================================================
USE MPI
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ichck_cr
TYPE (objstruc) :: obj,try,tryp
TYPE (datastruc) :: dat

dc = obj%par(idxc(2:)) - obj%par(idxc(1:idxdim-1))
dr = obj%par(idxr(2:)) - obj%par(idxr(1:idxdim-1))

ichck_cr = 0
DO i=1,idxdim-1
   IF(dc(i) >= 0._RP)THEN
      IF(dr(i) <  0._RP) ichck_cr = 1
   ELSE
      IF(dr(i) >= 0._RP) ichck_cr = 1
   ENDIF
ENDDO

128 FORMAT(20F16.6)
END SUBROUTINE COUPLE_CR


!=======================================================================

SUBROUTINE REF_NLAY3(thd,m_rg,freq,ref,cw,rw,nfrq,NANG,NFP)
!
!
!     CANNOT HANDLE HALFSPACE AT THIS POINT
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
v=1._RP/CMPLX(1._RP/c,alf/dB2nep) ! force radiation condition to be satisfied
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

   zin = znlay*(znlay1-CMPLX(0._RP,1._RP)*znlay*TRANSPOSE(tankd))/ &
         (znlay - CMPLX(0._RP,1._RP)*znlay1*TRANSPOSE(tankd))

   DO i = NLAY,3,-1
      DO j = 1,nfrq
         thtmp(j,:) = sin(th(i-1,:))
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

      zin = zm1*(zin -CMPLX(0._RP,1._RP)*zm1*TRANSPOSE(tankd))/ &
           (zm1 - CMPLX(0._RP,1._RP)*zin*TRANSPOSE(tankd));
   ENDDO
   z1=z(1,:,:)

   reftmp = (zin-z1)/(zin+z1)
ENDIF

ref = ABS(TRANSPOSE(reftmp))

DEALLOCATE(c,alf,r,d,th1,v,v2,zz,th,z,freq2,k,reftmp,znlay, &
           znlay1,zin,zm1,z1,tankd,thtmp,ktmp)

RETURN
END SUBROUTINE REF_NLAY3

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
cosc = CMPLX(COS(x)*COSH(y),-SIN(x)*SINH(y))
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
sinc = CMPLX(SIN(x)*COSH(y),COS(x)*SINH(y))
RETURN
END FUNCTION

!==============================================================================

FUNCTION CACOS(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CACOS
COMPLEX(KIND=RP) :: z

CACOS = -CMPLX(0._RP,1._RP)*LOG(z+CMPLX(0._RP,1._RP)*SQRT(1._RP-z*z))
!CACOS = (PI/2._RP)+(CMPLX(0._RP,1._RP)* &
!        LOG(CMPLX(0._RP,1._RP)*z+SQRT(1._RP-(z*z))))

RETURN
END FUNCTION CACOS

!==============================================================================

FUNCTION CTAN(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CTAN
COMPLEX(KIND=RP) :: z

CTAN =  -CMPLX(0._RP,1._RP)*(EXP( CMPLX(0._RP,1._RP)*z) -EXP(-CMPLX(0._RP,1._RP)*z)) &
                            /(EXP( CMPLX(0._RP,1._RP)*z)+EXP(-CMPLX(0._RP,1._RP)*z))
RETURN
END FUNCTION CTAN

!=======================================================================

SUBROUTINE RESULTS(sample,logZ)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB)                   :: i,j
REAL(KIND=RP),DIMENSION(IMAX,NFP+5):: sample
REAL(KIND=RP)                      :: mean(NFP),tmp(NFP),std(NFP),wt,logZ

mean = 0._RP
tmp  = 0._RP
std  = 0._RP

Do i = 1,IMAX
    
   !
   ! Proportional weights
   !
   wt = exp(sample(i,3) - logZ); ! logwt - logZ

   !
   ! 1st and 2nd moments of parameters
   !
   Do j = 1,NFP
      mean(j) = mean(j) + wt * sample(i,5+j)
      tmp(j)  = tmp(j)  + wt * sample(i,5+j) * sample(i,5+j)
   ENDDO
ENDDO
Do j = 1,NFP
   std(j) = sqrt(tmp(j)-mean(j)*mean(j))
   WRITE(*,204) j,mean(j), std(j)
ENDDO

204 FORMAT(2x,'mean(',1I2,') =',F12.6,' std = ', F12.6)
END SUBROUTINE RESULTS
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

Subroutine CHOLESKY(A,n,np,L)
!=======================================================================
!  Cholesky decomposition of symmetric, positive-definite matix A
!  into lower-triangular matrix L. Modified from Numerical Recipes.
!-----------------------------------------------------------------------
USE DATA_TYPE
IMPLICIT NONE
INTEGER(KIND=IB):: n,np,i,j,k
REAL(KIND=RP)   :: A(np,np),A2(np,np),L(np,np),p(n),summ

A2 = A
DO i=1,n
   DO j=1,n
      IF (i .GT. j) A2(i,j) = 0._RP
   ENDDO
ENDDO

DO i=1,n
   DO j=i,n
      summ = A2(i,j)
      DO k=i-1,1,-1
         summ = summ-A2(i,k)*A2(j,k)
      ENDDO
      IF (i .EQ. j) THEN
         IF (summ .LE. 0.) THEN
            WRITE(6,*) 'Cholesky Decomp Failed ',summ
            STOP
         ENDIF
         p(i) = SQRT(summ)
      ELSE
         A2(j,i) = summ/p(i)
      ENDIF
   ENDDO
ENDDO

DO i=1,n
   DO j=1,n
      IF (i .GT. j) L(i,j) = A2(i,j)
      IF (i .EQ. j) L(i,i) = p(i)
   ENDDO
ENDDO

RETURN
END SUBROUTINE CHOLESKY
!=======================================================================

SUBROUTINE LINROT(dat,CdCi)
!!
!!  Estimates linearized model covariance matrix 
!!
!=======================================================================
USE MPI
USE NR
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ifq,ip,ibest,ithread,ichunk,ntot

TYPE (objstruc),DIMENSION(2):: obj          ! Objects in likelihood box
TYPE (objstruc)             :: try,tryp     ! Objects in likelihood box
TYPE (datastruc)            :: dat          ! Data

REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND)     :: CdCi

REAL(KIND=RP)                                :: best,test
REAL(KIND=RP)                                :: t1,t2
REAL(KIND=RP),DIMENSION(NFP)                 :: ran_uni2

REAL(KIND=RP)                                :: dm
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: dm_start
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: p1,p2
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: dat1,dat2
REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: dRdm

REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Jac,JactJac,Mpriorinv,Temp
REAL(KIND=RP), DIMENSION(:),   ALLOCATABLE   :: bufJac,WW,Jacscale
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: VV,Winv
REAL(KIND=RP)                                :: Wmax
REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE   :: Ctmp

IF(rank == src)WRITE(*,*)'Computing linear rotation estimate...'
ntot = NBAND*NANG

ALLOCATE( dat1(ntot),dat2(ntot),dRdm(NFP,ntot,NDM),       &
          dm_start(NFP),p1(NBAND,NANG),p2(NBAND,NANG),    &
          Jac(ntot,NFP),JactJac(ntot,NFP),bufJac(ntot),   &
          Jacscale(NFP),Mpriorinv(NFP,NFP),Temp(ntot,NFP),&
          Ctmp(NFP,NFP),VV(NFP,NFP),Winv(NFP,NFP),WW(NFP) )

!dm_start  = 1._RP
dm_start  = maxpert/4._RP


CHUNK_BG = FLOOR(FLOAT(NFP)/FLOAT(NTHREAD))+1
CHUNK_SM = FLOOR(FLOAT(NFP)/FLOAT(NTHREAD))
REST  = NFP - (CHUNK_SM*NTHREAD)

!-----------------------------------------------------------------------
!   IF preload Jac
!-----------------------------------------------------------------------
IF (IPREJAC .EQ. 0) THEN

!-----------------------------------------------------------------------
!   MPI IF
!-----------------------------------------------------------------------
IF (rank .EQ. src) THEN

WRITE(*,*)'SLAVES 1, SLAVES 2, MASTER: ',REST*CHUNK_BG,(NTHREAD-1-REST)*CHUNK_SM, CHUNK_SM
WRITE(*,*)'CHUNK_BG =',CHUNK_BG, 'CHUNK_SM =', CHUNK_SM, 'REST =',REST

DO ip=1,CHUNK_SM
   dm = dm_start(ip) ! Estimate deriv for range of dm values
   DO i=1,6!NDM 
      obj(1)%par = dat%start
      obj(1)%par(ip) = dat%start(ip)+dm
      CALL LOGLHOOD(obj(1),dat)
      p1 = dat%Rrep
      DO ifq=1,NBAND
         IF(ICOV == 0)THEN
            dat1((ifq-1)*NANG+1:(ifq)*NANG) = p1(ifq,:)/sdev(ifq)
         ELSE
            dat1((ifq-1)*NANG+1:(ifq)*NANG) = & 
            MATMUL(CdCi(:,:,ifq),p1(ifq,:))
         ENDIF
      ENDDO        
      obj(2)%par = dat%start
      obj(2)%par(ip) = dat%start(ip)-dm
      CALL LOGLHOOD(obj(2),dat)
      p2 = dat%Rrep
      DO ifq=1,NBAND
         IF(ICOV == 0)THEN
            dat2((ifq-1)*NANG+1:(ifq)*NANG) = p2(ifq,:)/sdev(ifq)
         ELSE
            dat2((ifq-1)*NANG+1:(ifq)*NANG) = & 
            MATMUL(CdCi(:,:,ifq),p2(ifq,:))
         ENDIF
      ENDDO       
      IF (abs(dm) .GT. 0._RP) THEN 
          dRdm(ip,:,i) = (dat1(:)-dat2(:))/(2._RP*dm)
      ELSE 
          dRdm(ip,:,i) = 0._RP
      ENDIF
      dm = dm/1.25_RP
   ENDDO 

   DO j=1,ntot   ! For each datum, choose best derivative estimate
      best = 1.e10_RP
      ibest = 1
      DO i=1,NDM-2 
         IF ((dRdm(ip,j,i)*dRdm(ip,j,i+1)*dRdm(ip,j,i+2)) &
              .NE. 0._RP) THEN
            t1 = abs( dRdm(ip,j,i+0)/ dRdm(ip,j,i+1)-1._RP)
            t2 = abs( dRdm(ip,j,i+1)/ dRdm(ip,j,i+2)-1._RP)  
            test = maxval((/t1,t2/))
         ELSE 
            test = 1.e10_RP
         ENDIF
         IF ((test < best) .and. (test /= 0._RP)) then        
            best  = test
            ibest = i
         ENDIF 
      ENDDO  
      Jac(j,ip) = dRdm(ip,j,ibest)  ! Best deriv into Jacobian               
   ENDDO  

ENDDO

IF(REST /= 0)THEN
   DO ithread = 1,REST
      DO ichunk = 1,CHUNK_BG
         from  = ithread
         tag   = MPI_ANY_TAG
         icount = SIZE(Jac,1)
         CALL MPI_RECV(bufJac, icount, MPI_DOUBLE_PRECISION, from, &
                       tag, MPI_COMM_WORLD, status, ierr )
         Jac(:,CHUNK_SM+((ithread-1)*CHUNK_BG+1)+(ichunk-1)) = bufJac
      ENDDO
   ENDDO
   DO ithread = REST+1,NTHREAD-1
      DO ichunk = 1,CHUNK_SM
         from  = ithread
         tag   = MPI_ANY_TAG
         icount = SIZE(Jac,1)
         CALL MPI_RECV(bufJac, icount, MPI_DOUBLE_PRECISION, from, &
                       tag, MPI_COMM_WORLD, status, ierr )
         Jac(:,CHUNK_SM+(REST*CHUNK_BG)+(ithread-rest-1)*CHUNK_SM+ichunk) = bufJac
      ENDDO
   ENDDO
ELSE
   DO ithread = 1,NTHREAD-1
      DO ichunk = 1,CHUNK_SM
         from  = ithread
         tag   = MPI_ANY_TAG
         icount = SIZE(Jac,1)
         CALL MPI_RECV(bufJac, icount, MPI_DOUBLE_PRECISION, from, &
                       tag, MPI_COMM_WORLD, status, ierr )
         Jac(:,CHUNK_SM+((ithread-1)*CHUNK_SM+1)+(ichunk-1)) = bufJac
      ENDDO
   ENDDO
ENDIF

OPEN(UNIT=20,FILE=jacfile,FORM='formatted',STATUS='REPLACE', &
     ACTION='WRITE',RECL=1024)
DO i = 1,ntot
    WRITE(20,*) Jac(i,:)
ENDDO
CLOSE(20)

!-----------------------------------------------------------------------
! MPI IF
!-----------------------------------------------------------------------
ELSEIF ((rank .NE. 0)) THEN

   IF(rank <= REST)THEN
      CHUNK_SL = CHUNK_BG
      CHUNK_BEG = CHUNK_SM+((rank-1)*CHUNK_SL+1)
      CHUNK_END = CHUNK_SM+(rank*CHUNK_SL)
   ELSE
      CHUNK_SL = CHUNK_SM
      CHUNK_BEG = CHUNK_SM+(REST*CHUNK_BG)+((rank-REST-1)*CHUNK_SM+1)
      CHUNK_END = CHUNK_SM+(REST*CHUNK_BG)+((rank-REST)*CHUNK_SL)
   ENDIF

   DO ip = CHUNK_BEG,CHUNK_END

   dm = dm_start(ip) ! Estimate deriv for range of dm values
   DO i=1,NDM 
      obj(1)%par = dat%start
      obj(1)%par(ip) = dat%start(ip)+dm
      CALL LOGLHOOD(obj(1),dat)
      p1 = dat%Rrep
      DO ifq=1,NBAND
         IF(ICOV == 0)THEN
            dat1((ifq-1)*NANG+1:(ifq)*NANG) = p1(ifq,:)/sdev(ifq)
         ELSE
            dat1((ifq-1)*NANG+1:(ifq)*NANG) = & 
            MATMUL(CdCi(:,:,ifq),p1(ifq,:))
         ENDIF
      ENDDO        
      obj(2)%par = dat%start
      obj(2)%par(ip) = dat%start(ip)-dm
      CALL LOGLHOOD(obj(2),dat)
      p2 = dat%Rrep
      DO ifq=1,NBAND
         IF(ICOV == 0)THEN
            dat2((ifq-1)*NANG+1:(ifq)*NANG) = p2(ifq,:)/sdev(ifq)
         ELSE
            dat2((ifq-1)*NANG+1:(ifq)*NANG) = & 
            MATMUL(CdCi(:,:,ifq),p2(ifq,:))
         ENDIF
      ENDDO       
      IF (abs(dm) .GT. 0._RP) THEN 
          dRdm(ip,:,i) = (dat1(:)-dat2(:))/(2._RP*dm)
      ELSE 
          dRdm(ip,:,i) = 0._RP
      ENDIF
      dm = dm/1.25_RP
   ENDDO 

   DO j=1,ntot   ! For each datum, choose best derivive estimate
      best = 1.e10_RP
      ibest = 1
      DO i=1,NDM-2 
         IF ((dRdm(ip,j,i)*dRdm(ip,j,i+1)*dRdm(ip,j,i+2)) &
              .NE. 0._RP) THEN
            t1 = ABS( dRdm(ip,j,i+0)/ dRdm(ip,j,i+1)-1._RP)
            t2 = ABS( dRdm(ip,j,i+1)/ dRdm(ip,j,i+2)-1._RP)  
            test = MAXVAL((/t1,t2/))
         ELSE 
            test = 1.e10_RP
         ENDIF
         IF ((test < best) .AND. (test /= 0._RP)) THEN
            best  = test
            ibest = i
         ENDIF 
      ENDDO  
      Jac(j,ip) = dRdm(ip,j,ibest)  ! Best deriv into Jacobian               
   ENDDO
   
!
!  ENDDO ip loop
!
   ENDDO


   DO ip = CHUNK_BEG,CHUNK_END
   
      icount  = SIZE(Jac,1)
      bufJac = Jac(:,ip)
      CALL MPI_SSEND( bufJac, icount, MPI_DOUBLE_PRECISION, src, & 
                      rank, MPI_COMM_WORLD, ierr )

   ENDDO

!------------------------------------------------------------------------  
!   MPI END IF
!------------------------------------------------------------------------  
ENDIF

CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
DO j=1,NFP
    CALL MPI_BCAST( Jac(:,j), ntot, MPI_DOUBLE_PRECISION, src, &
                    MPI_COMM_WORLD, ierr )
ENDDO
CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

ELSE

    IF(rank .EQ. src) WRITE(6,*) &
    'IPREJAC SET... LOADING PRE COMPUTED JACOBOIAN'
    OPEN(UNIT=20,FILE='preJac.txt',FORM='formatted',STATUS='OLD',&
         ACTION='READ')

    DO i = 1,ntot
        READ(20,*) (Jac(i,j),j=1,NFP)
    ENDDO
    CLOSE(20)

         !
ENDIF    ! ENDIF preload Jac
         !

DO i=1,NFP   ! Scale columns of Jacobian to unity for stability      
   Jacscale(i) = 0._RP
   DO j=1,ntot  
      Jacscale(i) = Jacscale(i)+(ABS(Jac(j,i))*ABS(Jac(j,i)))
   ENDDO
   Jacscale(i) = SQRT(Jacscale(i))
   DO j=1,ntot
      Jac(j,i) = Jac(j,i)/Jacscale(i)
   ENDDO
ENDDO

JactJac = MATMUL(TRANSPOSE(Jac),Jac)

Mpriorinv = 0._RP
DO i=1,NFP 
!   Mpriorinv(i,i) = 1./(scale(i)*(maxlim(i)-minlim(i))/2.)**2
   Mpriorinv(i,i) = 1._RP/(Jacscale(i)*(maxlim(i)-minlim(i))/4._RP)**2._RP
ENDDO
      
Ctmp = JactJac + Mpriorinv
CALL SVDCMP(Ctmp,WW,VV)

Wmax = maxval(WW,1)
Winv = 0._RP
DO i=1,NFP
   IF (WW(i)/Wmax > 1.e-7_RP) THEN
      Winv(i,i) = 1._RP/WW(i)
   ELSE
      WRITE(6,*) 'Warning: Small singular value'
      WRITE(6,*) WW/Wmax
   ENDIF
ENDDO
CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

dat%Cov0 = MATMUL(MATMUL(VV,Winv),TRANSPOSE(VV))

DO i=1,NFP 
   DO j=1,NFP
      dat%Cov0(i,j) = dat%Cov0(i,j)/(Jacscale(i)*Jacscale(j))
   ENDDO
ENDDO

Ctmp = dat%Cov0
CALL SVDCMP(Ctmp,WW,VV)
dat%Cov = dat%Cov0
dat%sdevm = SQRT(WW)
dat%VV = VV

!!
!! Compute approximate rotated bounds
!!
maxlimp=-1.e308_RP
minlimp= 1.e308_RP
                        
DO i=1,10000*NFP
   CALL RANDOM_NUMBER(ran_uni2)
   try%par = minlim + (ran_uni2*maxpert)
   tryp%par = MATMUL(TRANSPOSE(dat%VV),try%par)
   DO j=1,NFP
      IF(tryp%par(j) > maxlimp(j)) maxlimp(j)=tryp%par(j)
      IF(tryp%par(j) < minlimp(j)) minlimp(j)=tryp%par(j)
   ENDDO
ENDDO
maxpertp = maxlimp - minlimp
wp = maxpertp/wscale

IF(rank == src)THEN
   WRITE(6,*) 'Cov0 main diagonal:'
   DO i=1,NFP
      WRITE(6,FMT=204) i,dat%start(i),dat%sdevm(i),SQRT(dat%Cov(i,i))
   ENDDO
   OPEN(UNIT=41,FILE=lincovfile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
   DO i=1,NFP
      WRITE(41,FMT=205) dat%Cov0(i,:)
   ENDDO
   CLOSE(41)
   WRITE(*,*)'Done computing linear rotation estimate.'
ENDIF

203 FORMAT(i7,f16.5)
204 FORMAT(i7,6(f16.5))      
205 FORMAT(50F16.8)
 
DEALLOCATE( dat1,dat2,dRdm,p1,p2,Jac,JactJac,bufJac,Mpriorinv, &
            Ctmp,VV,Winv,WW )
CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
RETURN
END SUBROUTINE LINROT

!=======================================================================

SUBROUTINE NONLINROT(objold,dat)
!!
!!  Estimates nonlinear model covariance matrix by sampling
!!
!=======================================================================
USE MPI
USE NEST_COM
USE NR
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ikeep
TYPE (datastruc):: dat
TYPE (covstruc) :: msums
TYPE (objstruc) :: objold,obj1old,obj2old,obj1new,obj2new,try,tryp
REAL(KIND=RP)                   :: ran_uni,logy
REAL(KIND=RP),DIMENSION(NFP)    :: ran_uni2
REAL(KIND=RP),DIMENSION(NFP)    :: WW
REAL(KIND=RP),DIMENSION(NFP,NFP):: Ctmp,VV
REAL(KIND=RP),DIMENSION(NKEEP,NFP+1)        :: mkeep1,mkeep2
REAL(KIND=RP),DIMENSION(NKEEP*NTHREAD,NFP+1):: mkeep1M,mkeep2M

IF(rank == src)WRITE(*,*)'Computing non-linear rotation estimate (burn-in)...'


obj1old = objold
obj2old = objold

msums%jsamp = 0
msums%ksamp = 0
DO
   DO ikeep = 1,NKEEP
      obj1new = obj1old
      obj2new = obj2old

      !!
      !! Sample1 
      !!
!      IF(rank == src)PRINT*,ikeep
!      PRINT*,rank,'obj1',obj1old%par
!      PRINT*,rank,'obj2',obj1old%par
      CALL RANDOM_NUMBER(ran_uni)
      logy = obj1old%logL + LOG(ran_uni)
      CALL EXPLORE_SLICE_ROT(obj1new,dat,logy)

      !!
      !! Sample2 
      !!
      CALL RANDOM_NUMBER(ran_uni)
      logy = obj2old%logL + LOG(ran_uni)
      CALL EXPLORE_SLICE_ROT(obj2new,dat,logy)

      mkeep1(ikeep,:) = (/ obj1new%logL, obj1new%par /)
      mkeep2(ikeep,:) = (/ obj1new%logL, obj2new%par /)
      obj1old = obj1new
      obj2old = obj2new
   ENDDO

   CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
   CALL SR_SAMPLE(mkeep1,mkeep2,mkeep1M,mkeep2M)
   msums%jsamp = msums%jsamp + SIZE(mkeep1M,1)
   msums%ksamp = msums%ksamp + (2*SIZE(mkeep1M,1))

   IF(rank == src)THEN
      CALL COVCOR(mkeep1M,mkeep2M,msums,dat)
      WRITE(*,207) msums%jsamp,msums%dcor,obj1new%logL,obj1new%par
   ENDIF

   DO i=1,NFP
      CALL MPI_BCAST( dat%Cov(:,i), NFP, MPI_DOUBLE_PRECISION, src, &
                      MPI_COMM_WORLD, ierr )
   ENDDO
   CALL MPI_BCAST(msums%iburnin,1,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)

   !!
   !! Compute rotation matrix
   !!
   Ctmp = dat%Cov
   CALL SVDCMP(Ctmp,WW,VV)
   dat%VV = VV
   dat%sdevm = SQRT(WW)

   !!
   !! Compute approximate rotated bounds
   !!
   maxlimp = -1.e308_RP
   minlimp =  1.e308_RP
                        
   DO i=1,10000*NFP
      CALL RANDOM_NUMBER(ran_uni2)
      try%par = minlim + (ran_uni2*maxpert)
      tryp%par = MATMUL(TRANSPOSE(dat%VV),try%par)
      DO j=1,NFP
         IF(tryp%par(j) > maxlimp(j)) maxlimp(j)=tryp%par(j)
         IF(tryp%par(j) < minlimp(j)) minlimp(j)=tryp%par(j)
      ENDDO
   ENDDO
   maxpertp = maxlimp - minlimp
   wp = maxpertp/wscale

   IF(msums%iburnin == 1)EXIT
ENDDO
CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

IF(rank == src)THEN
   WRITE(6,*) ''
   WRITE(6,*) 'Burn-in completed.'
   WRITE(6,*) ''
   OPEN(UNIT=42,FILE=nonlincovfile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
   DO i=1,NFP
      WRITE(42,206) dat%Cov(i,:)
   ENDDO
   CLOSE(42)
ENDIF

206 FORMAT(50F16.8)
207 FORMAT(I4,50F12.4)
RETURN
END SUBROUTINE NONLINROT

!=======================================================================

SUBROUTINE COVCOR(mkeep1M,mkeep2M,msums,dat)
!!
!!  Compute covariances and correlations for two independent samples.
!!
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j
TYPE (datastruc):: dat
TYPE (covstruc) :: msums
REAL(KIND=RP),DIMENSION(NFP,NFP)  :: Cov1 = 0._RP
REAL(KIND=RP),DIMENSION(NFP,NFP)  :: Cov2 = 0._RP
REAL(KIND=RP),DIMENSION(NFP,NFP)  :: Cov3 = 0._RP
REAL(KIND=RP),DIMENSION(NFP,NFP)  :: Cor1 = 0._RP
REAL(KIND=RP),DIMENSION(NFP,NFP)  :: Cor2 = 0._RP
REAL(KIND=RP),DIMENSION(NKEEP*NTHREAD,NFP+1):: mkeep1M,mkeep2M

DO i=1,NFP 
   msums%msum1(i) = msums%msum1(i)+SUM(mkeep1M(:,i+1))         
   msums%msum2(i) = msums%msum2(i)+SUM(mkeep2M(:,i+1))
   msums%msum3(i) = msums%msum3(i)+SUM(mkeep1M(:,i+1))+SUM(mkeep2M(:,i+1))
   DO j=1,NFP 
      msums%mcsum1(i,j) = msums%mcsum1(i,j)+SUM(mkeep1M(:,i+1)*mkeep1M(:,j+1))
      msums%mcsum2(i,j) = msums%mcsum2(i,j)+SUM(mkeep2M(:,i+1)*mkeep2M(:,j+1))
      msums%mcsum3(i,j) = msums%mcsum3(i,j)+SUM(mkeep1M(:,i+1)*mkeep1M(:,j+1)) &
                                           +SUM(mkeep2M(:,i+1)*mkeep2M(:,j+1))
   ENDDO
ENDDO                       
msums%mave1   = msums%msum1 /REAL(msums%jsamp,KIND=RP)
msums%mave2   = msums%msum2 /REAL(msums%jsamp,KIND=RP)
msums%mave3   = msums%msum3 /REAL(2*msums%jsamp,KIND=RP)
msums%mcross1 = msums%mcsum1/REAL(msums%jsamp,KIND=RP)
msums%mcross2 = msums%mcsum2/REAL(msums%jsamp,KIND=RP) 
msums%mcross3 = msums%mcsum3/REAL(2*msums%jsamp,KIND=RP) 
DO i=1,NFP
   DO j=1,NFP 
      Cov1(i,j) = msums%mcross1(i,j)-msums%mave1(i)*msums%mave1(j)
      Cov2(i,j) = msums%mcross2(i,j)-msums%mave2(i)*msums%mave2(j) 
      Cov3(i,j) = msums%mcross3(i,j)-msums%mave3(i)*msums%mave3(j)
   ENDDO
ENDDO 

DO i=1,NFP 
   DO j=1,NFP 
      IF (Cov1(i,i)*Cov1(j,j) <= 0.) THEN 
         Cor1(i,j) = 0.
      ELSE 
         Cor1(i,j) = Cov1(i,j)/SQRT(Cov1(i,i)*Cov1(j,j))
      ENDIF 
      IF (Cov2(i,i)*Cov2(j,j) <= 0.) THEN
         Cor2(i,j) = 0.
      ELSE                  
         Cor2(i,j) = Cov2(i,j)/SQRT(Cov2(i,i)*Cov2(j,j))
      ENDIF         
   ENDDO
ENDDO

!-----------------------------------------------------------------------
!  If sample correlations agree, update Chol based on a weighted
!  average of linear and nonlinear covariance matrices.
!-----------------------------------------------------------------------    

msums%dcor = MAXVAL(ABS(Cor1-Cor2))  

IF (msums%iburnin == 1) THEN 
   dat%Cov = Cov3
ENDIF
IF ((msums%dcor < corconv1) .AND. (msums%dcor > corconv2) .AND. &
    (msums%jsamp >= NFP*NFP) .AND. msums%iburnin == 0) THEN

   dat%Cov = (corconv1-msums%dcor)/(1._RP-corconv2)*Cov3 + &
             (msums%dcor-corconv2)/(1._RP-corconv2)*dat%Cov0

ENDIF
IF (((msums%dcor <= corconv2) .AND. (msums%jsamp >= NFP*NFP) .OR. (msums%jsamp > 100000)) &
     .AND. msums%iburnin == 0) THEN

   dat%Cov = Cov3
   msums%ksamp = 0
   msums%iburnin = 1

ENDIF

RETURN
END SUBROUTINE COVCOR

!=======================================================================

SUBROUTINE SR_SAMPLE(mkeep1,mkeep2,mkeep1M,mkeep2M)
!!
!!  Send and receive samples
!!
!=======================================================================
USE MPI
USE NEST_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j
REAL(KIND=RP),DIMENSION(NKEEP,NFP+1)        :: mkeep1,mkeep2
REAL(KIND=RP),DIMENSION(NKEEP*NTHREAD,NFP+1):: mkeep1M,mkeep2M
REAL(KIND=RP),DIMENSION(NKEEP*(NFP+1)):: buf2

IF(rank == src)THEN
   mkeep1M(1:NKEEP,:) = mkeep1
   mkeep2M(1:NKEEP,:) = mkeep2
   DO i = 1,NTHREAD-1

      from  = i
      tag   = MPI_ANY_TAG
      icount = NKEEP*(NFP+1)
      CALL MPI_RECV(buf2, icount, MPI_DOUBLE_PRECISION, from,&
                    tag, MPI_COMM_WORLD, status, ierr )
      mkeep1M(i*NKEEP+1:(i+1)*(NKEEP),:) = RESHAPE(buf2,(/ NKEEP,NFP+1 /))

      CALL MPI_RECV(buf2, icount, MPI_DOUBLE_PRECISION, from,&
                    tag, MPI_COMM_WORLD, status, ierr )
      mkeep2M(i*NKEEP+1:(i+1)*(NKEEP),:) = RESHAPE(buf2,(/ NKEEP,NFP+1 /))

   ENDDO

ELSE

    icount = NKEEP*(NFP+1)
    buf2 = PACK(mkeep1,mkeep1 /= SQRT(-1._RP))
    CALL MPI_SSEND( buf2, icount, MPI_DOUBLE_PRECISION, src, &
                    rank, MPI_COMM_WORLD, ierr )
    buf2 = PACK(mkeep2,mkeep2 /= SQRT(-1._RP))
    CALL MPI_SSEND( buf2, icount, MPI_DOUBLE_PRECISION, src, &
                    rank, MPI_COMM_WORLD, ierr )

ENDIF



RETURN
END SUBROUTINE SR_SAMPLE
!=======================================================================
! This is the end my fiend...
! EOF
