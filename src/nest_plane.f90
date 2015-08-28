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
   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.0D0)
   REAL(KIND=RP),   PARAMETER :: PI  = 3.141592653589793238462643383279502884197_RP
END MODULE DATA_TYPE

!==============================================================================
MODULE NEST_COM
   USE DATA_TYPE
   IMPLICIT NONE

!
! General switches
!
   INTEGER(KIND=IB), PARAMETER :: ICOV    = 1            ! Data cov. switch
   INTEGER(KIND=IB), PARAMETER :: ISEEDST = 0            ! 1: Seed good starting model
   INTEGER(KIND=IB), PARAMETER :: NAVEF   = 12           ! # freq per band
   REAL(KIND=RP),    PARAMETER :: frbw    = 1._RP/20._RP ! Frac. bandwidth for freq ave.

!!
!!  SIM B
!!
   INTEGER(KIND=IB), PARAMETER :: NANG   = 131
   INTEGER(KIND=IB), PARAMETER :: NFP    = 23
   REAL(KIND=RP),    PARAMETER :: cw     = 1511.0_RP
   REAL(KIND=RP),    PARAMETER :: rw     = 1.029_RP
   CHARACTER(len=64) :: infile1    = 'sim_B_1_10_50_5lay.txt'
   CHARACTER(LEN=64) :: logfile    = 'sim_B_1_10_50_5lay_FGS.log'
   CHARACTER(len=64) :: covfile    = 'sim_B_1_10_50_4lay_cov.txt'
   CHARACTER(len=64) :: mapfile    = 'sim_B_1_10_50_5lay_map.dat'
   CHARACTER(len=64) :: sdfile     = 'sim_B_1_10_50_exstd.txt'
   CHARACTER(len=64) :: samplefile = 'sim_B_1_10_50_5lay_sample.txt'

   INTEGER(KIND=IB), PARAMETER :: NBAND  = 14
   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 1000._RP, 1200._RP, 1400._RP, 1683._RP, 1889._RP, 2121._RP, &
                                                  2381._RP, 2672._RP, 3000._RP, 3367._RP, 3779._RP, 4242._RP, &
                                                  4762._RP, 5345._RP /)

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
!   INTEGER(KIND=IB), PARAMETER :: NANG = 131
!   INTEGER(KIND=IB), PARAMETER :: NFP  = 3
!   REAL(KIND=RP),    PARAMETER :: cw   = 1510.78_RP
!   REAL(KIND=RP),    PARAMETER :: rw   = 1.029_RP
!   CHARACTER(len=64) :: infile1    = 'x_s02_1_3_16_0lay.txt'
!   CHARACTER(LEN=64) :: logfile    = 'x_s02_1_3_16_0lay_FGS.log'
!   CHARACTER(len=64) :: covfile    = 'x_s02_1_3_16_4lay_cov.txt'
!   CHARACTER(len=64) :: mapfile    = 'x_s02_1_3_16_0laycov_map.dat'
!   CHARACTER(len=64) :: sdfile     = 'x_s02_1_3_16_exstd.txt'
!   CHARACTER(len=64) :: samplefile = 'x_s02_1_3_16_0lay_sample.txt'
!
!   INTEGER(KIND=IB), PARAMETER :: NBAND  = 8
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/ 300._RP, 400._RP,  504._RP,  635._RP,  800._RP, 1008._RP, 1270._RP, & 
!                                                 1600._RP /)

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP), DIMENSION(NFP):: minlim  = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxlim  = 0._RP
   REAL(KIND=RP), DIMENSION(NFP):: maxpert = 0._RP

!!
!!  Nested sampling specific parameters
!!
   INTEGER(KIND=IB), PARAMETER :: N      = NFP*NFP*30_IB ! # objects in box
!   INTEGER(KIND=IB), PARAMETER :: N      = 10_IB      ! # objects in box

   REAL(KIND=RP),    PARAMETER :: endfac = 2._RP      ! Termination factor
   REAL(KIND=RP),    PARAMETER :: logtol = 1E-2_RP    ! Termination tolerance (log space)
   REAL(KIND=RP),DIMENSION(NFP):: step   = 0.1_RP     ! Initial guess suitable step-size in (0,1)
   INTEGER(KIND=IB), PARAMETER :: IMAX   = 1E6_IB     ! # iterations (# nested boxes)
   INTEGER(KIND=IB), PARAMETER :: JMAX   = 5_IB       ! Max # MCMC iterations
                                                      ! (must be large enough so
                                                      ! that correlation with
                                                      ! the original point is
                                                      ! entirely lost)

!!
!!  Structures for objects and data 
!!
   TYPE :: objstruc
      REAL(KIND=RP),DIMENSION(NFP):: par   = 0._RP   ! Forward parameters
      REAL(KIND=RP),DIMENSION(NFP):: prior = 0._RP   ! Prior parameters drawn from [0,1]
      REAL(KIND=RP)               :: logwt           ! log weights
      REAL(KIND=RP)               :: logL            ! log likelihood
   END TYPE objstruc

   TYPE :: datastruc
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Robs = 0._RP  ! Observed data
      REAL(KIND=RP),DIMENSION(NANG)           :: angobs = 0._RP! Observed angles
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rrep = 0._RP  ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: res  = 0._RP  ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rex  = 0._RP  ! Index for bad points/data gaps
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cdi  = 0._RP  ! Inverse data covariance matrices
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cd   = 0._RP  ! Data covariance matrices
      REAL(KIND=RP),DIMENSION(NFP)            :: start = 0._RP ! Good starting object
      REAL(KIND=RP),DIMENSION(NFP)            :: true  = 0._RP ! True object (only
                                                     ! meaningful for simulations)
      INTEGER(KIND=IB), DIMENSION(NBAND)      :: NDPF          ! # data per freq
      INTEGER(KIND=IB):: NANG  = NANG                          ! # data/angles (struc copy)
      INTEGER(KIND=IB):: NBAND = NBAND                         ! # freq bands (struc copy)
   END TYPE datastruc
   REAL(KIND=RP),DIMENSION(NBAND):: sdev = 0._RP               ! Standard devs

!!
!!  MPI Global Variables
!!
   INTEGER(KIND=IB) :: rank,NTHREAD,count,ierr
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB

END MODULE NEST_COM
  
!=======================================================================

PROGRAM  NESTED_PLANE

!=======================================================================
USE MPI
USE NEST_COM
USE NR
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ifreq,inest,ikeep,iprint1,iprint2,icopy2
INTEGER(KIND=IB),PARAMETER  :: NKEEP  = 400
INTEGER(KIND=IB),PARAMETER  :: NPRINT = 100
INTEGER(KIND=IB),PARAMETER  :: NCOPY  = 10
INTEGER(KIND=IB), DIMENSION(2), PARAMETER :: iseed = 5678

TYPE (objstruc), DIMENSION(N)         :: obj     ! Objects in likelihood box
TYPE (datastruc)                      :: dat     ! Data
REAL(KIND=RP),   DIMENSION(NKEEP+N,NFP+5):: sample  ! Unweighted posterior sample
REAL(KIND=RP)                         :: ran_uni ! Uniform random number

INTEGER(KIND=IB)  :: iworst,ibest,icopy
REAL(KIND=RP)     :: logLstar      ! Hard likelihood constraint
REAL(KIND=RP)     :: logwidth      ! Log width of integration interval
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
REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: CdC = 0._RP
REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: CdCi = 0._RP

!!
!!  MPI variables
!!
REAL(KIND=RP)                    :: tstart, tend 

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )
tstart = MPI_WTIME()

!!------------------------------------------------------------------------
!!
!!  Print nested sampling parameters to screen for logging
!!
WRITE(6,*) 'Max # iterations:',IMAX
WRITE(6,*) 'Max # MCMC steps:',JMAX
WRITE(6,*) '# objects in box:',N
WRITE(6,*) 'Convergence when inest > N*H*',endfac

WRITE(6,*) infile1
WRITE(6,*) covfile
WRITE(6,*) sdfile
WRITE(6,*) samplefile
WRITE(6,*) ''

!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
WRITE(6,*) 'Loading data...'
OPEN(UNIT=20,FILE=infile1,FORM='formatted',STATUS='OLD',ACTION='READ')

READ(20,*) (minlim(j),j=1,NFP)
READ(20,*) (maxlim(j),j=1,NFP)
READ(20,*) (dat%true(j),j=1,NFP)
maxpert = maxlim - minlim

DO ifreq = 1,NBAND
    READ(20,*) (dat%Robs(ifreq,j),j=1,NANG)
ENDDO

READ(20,*) (dat%angobs(j),j=1,NANG)

DO ifreq = 1,NBAND
    READ(20,*) (dat%Rex(ifreq,j),j=1,NANG)
ENDDO
CLOSE(20)

WRITE(6,*) 'Number of angles:         ',NANG
WRITE(6,*) 'Read data.'
CALL FLUSH(6)

IF(ICOV .EQ. 1)THEN
   WRITE(*,*) 'Using data covariance matrix estimate.'
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
   WRITE(*,*) 'Using experimental error (one stdev per band) estimate.'
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

   CdCi(:,:,ifreq) = SNGL(CdC(:,:,ifreq))
   CALL SVDCMP(CdCi(:,:,ifreq),WWW,VVV)
   Wmax = MAXVAL(WWW,1)
   WWWi = 0.  
   DO i=1,NFP
       IF (WWW(i)/Wmax .GT. 1.e-7) THEN
             WWWi(i,i) = 1./WWW(i)
       ELSE
             WRITE(6,*) 'Warning: Small singular value'
             WRITE(6,*) WWW/Wmax
       ENDIF
   ENDDO
   CdCi(:,:,ifreq) = MATMUL(MATMUL(VVV,WWWi),TRANSPOSE(VVV))
ENDDO

IF(ISEEDST == 1)THEN
   WRITE(6,*)'Reading map starting model:',mapfile
   OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   READ(20,*) (dat%start(j),j=1,NFP)
   CLOSE(20)

!!
!!  Constrain good starting object to within bounds 
!!
   DO i=1,NFP
     IF (dat%start(i) < minlim(i)) THEN
        WRITE(6,*) 'Starting Model Outside Bounds',i,dat%start(i),minlim(i),maxlim(i)
        dat%start(i) = minlim(i) + ABS(dat%start(i)-minlim(i))
     ENDIF
     IF (dat%start(i) > maxlim(i)) THEN
        WRITE(6,*) 'Starting Model Outside Bounds',i,dat%start(i),minlim(i),maxlim(i)
        dat%start(i) = maxlim(i) - ABS(dat%start(i)-maxlim(i))
!        CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
!        STOP
     ENDIF
   ENDDO
ENDIF

!!
!!  File to save posterior sample
!!
OPEN(UNIT=40,FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)

!!
!!  Initialize variables for running sums:
!!
logZ  = -1.0E308_RP      ! Initializes Z to zero
H     = 0._RP            ! Initializes H to zero

!CALL RANDOM_SEED(PUT=iseed)
CALL RANDOM_SEED()

!------------------------------------------------------------------------
!
!           ************Nested Sampling************
!
! -----------------------------------------------------------------------

!!
!!  Draw objects (obj) randomly from prior
!!
DO i=1,N

   CALL PRIOR(obj(i),dat)

ENDDO

!!
!!  Seed likely starting object
!!
IF(ISEEDST == 1)THEN
   obj(1)%par   = dat%start
   obj(1)%prior = (obj(1)%par - minlim) / maxpert
   CALL LOGLHOOD(obj(1),dat)
ENDIF

logwidth = LOG(1.0_RP - EXP(-1.0_RP / N)) ! Outermost interval of prior mass

ikeep   = 1
iprint1 = 0
iprint2 = 0
icopy2  = 0
DO inest = 1,IMAX

   iworst = 1
   ibest  = 1
   DO i=2,N
      IF(obj(i)%logL < obj(iworst)%logL) iworst = i
      IF(obj(i)%logL > obj(ibest)%logL)  ibest  = i
   ENDDO
   obj(iworst)%logwt = logwidth + obj(iworst)%logL
   obj(ibest)%logwt = logwidth + obj(ibest)%logL

   !
   ! Update evidence Z and information H
   !
   logZnew = LOGPLUS(logZ, obj(iworst)%logwt)
   H = EXP(obj(iworst)%logwt - logZnew) * obj(iworst)%logL &
     + EXP(logZ - logZnew) * (H + logZ) - logZnew
   logZ = logZnew

   !
   ! Save sample
   !
   sample(ikeep,:) = (/ obj(iworst)%logL,logwidth, obj(iworst)%logwt, H, logZ, obj(iworst)%par /)
   IF(ikeep == NKEEP)THEN
      DO i=1,NKEEP
         WRITE(40,207) sample(i,:)  
      ENDDO
      ikeep = 0
      sample = 0._RP
   ENDIF
   ikeep  = ikeep + 1
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
      WRITE(6,205) inest,N*H*endfac,logZ,H,obj(iworst)%logL,obj(ibest)%logL
      iprint2 = iprint2 + 1 
   ENDIF

   !
   ! Evolve worst object into better one
   !
   CALL RANDOM_NUMBER(ran_uni)
   logLstar  = obj(iworst)%logL
   IF(inest >= icopy2*NCOPY)THEN
     icopy = AINT(ran_uni*N)+1
     obj(iworst) = obj(icopy)
!     WRITE(*,*) 'Performed copy'
     icopy2 = icopy2 + 1
   ENDIF
   CALL EXPLORE(obj(iworst),dat,logLstar)

   !
   ! Shrink interval for prior mass
   !
   logwidth = logwidth - 1.0_RP/DBLE(N)
   
   !
   ! Convergence
   !
   ! Bulk of posterior mass passed:
   IF(CEILING(endfac*N*H) < inest) EXIT       
   ! Max amount evidence can gain is small:
   IF(ABS(LOGPLUS(obj(ibest)%logwt,logZ)-logZ) < logtol)EXIT 


ENDDO

!
! Final correction:
!
logwidth = -DBLE(inest)/DBLE(N) - LOG(DBLE(N))

DO i = 1,N
   obj(i)%logwt = logwidth + obj(i)%logL
   !
   ! Correct H and Z:
   !
   logZnew = LOGPLUS(logZ, obj(i)%logwt)
   H = EXP(obj(i)%logwt - logZnew) * obj(i)%logL &
     + EXP(logZ - logZnew) * (H + logZ) - logZnew
   logZ = logZnew
   !
   ! Additional last N samples
   !
   sample(ikeep-1+i,:) = (/ obj(i)%logL,logwidth, obj(i)%logwt, H, logZ, obj(i)%par /)

ENDDO

!!
!!  Total time
!!
IF (rank == src) THEN
   tend = MPI_WTIME()
   WRITE(*,*) 'Total time: ',tend-tstart,'on',NTHREAD,'CPUs'
ENDIF

!
! Print results:
!
WRITE(*,201) inest
WRITE(*,202) logZ, SQRT(H/N)
WRITE(*,203) H, H/LOG(2._RP)
201 FORMAT(2x, 'No. iterates = ',I6)
202 FORMAT(2x, 'Evidence ln(Z) = ',F12.6,' +-', F12.6)
203 FORMAT(2x, 'Information H  = ',F12.6,'\t nats = ', F12.6,' bits')
205 FORMAT(I6,5(f16.4))      
206 FORMAT(A86)

DO i=1,ikeep-1+N
   WRITE(40,207) sample(i,:)  
ENDDO
CLOSE(40)
207   FORMAT(50F16.4)
!CALL RESULTS(sample,logZ)

CALL MPI_FINALIZE( ierr ) 
END PROGRAM NESTED_PLANE

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
obj%logL = -SUM(Etmp)

RETURN
END SUBROUTINE LOGLHOOD

!=======================================================================

SUBROUTINE PRIOR(obj,dat)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=RP) :: j
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP)    :: ran_uni(NFP)

   CALL RANDOM_NUMBER(ran_uni)
   obj%prior = ran_uni
   obj%par = minlim + maxpert * obj%prior
   CALL LOGLHOOD(obj,dat)

RETURN
END SUBROUTINE PRIOR

!=======================================================================

SUBROUTINE EXPLORE2(obj,dat,logLstar)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=RP):: accept(NFP),reject(NFP),i,j
TYPE (objstruc) :: obj,try
TYPE (datastruc):: dat
REAL(KIND=RP)   :: logLstar
REAL(KIND=RP)   :: ran_uni(NFP)

accept = 0      ! # MCMC acceptances
reject = 0      ! # MCMC rejections

DO j = 1,JMAX  ! MCMC counter (pre-judged # steps)
   !!
   !!  Trial object
   !!

   !!
   !!  Gaussian proposal
   !!
   DO i = 1,NFP  
      CALL GASDEVJ(ran_uni(i))
      try%prior(i) = obj%prior(i) + step(i) * ran_uni(i)      ! |move| < step
      try%prior(i) = try%prior(i) - FLOOR(try%prior(i))    ! wraparound to stay within (0,1)
      try%par(i)   = minlim(i) + maxpert(i) * try%prior(i) ! map to x

      CALL LOGLHOOD(try,dat)
!      PRINT*,sngl(try%logL),sngl(logLstar),sngl(obj%logL)

      !
      ! Accept if and only if within hard likelihood constraint
      !
      IF( try%logL > logLstar )THEN
         obj = try
         accept(i) = accept(i) + 1
      ELSE
         reject(i) = reject(i) + 1
      ENDIF

      !
      ! Refine step-size to let acceptance ratio converge around 50%
      !
      IF( accept(i) > reject(i) ) step(i) = step(i) * exp(1.0_RP / accept(i));
      IF( accept(i) < reject(i) ) step(i) = step(i) / exp(1.0_RP / reject(i));
   ENDDO   

   !
   ! Exit if already accepted lots of models
   !
!   IF(accept > CEILING(dble(JMAX)/2._RP)) EXIT

ENDDO
!PRINT*,'L= ',sngl(obj%logL),sngl(logLstar),sngl(step),sngl(accept+1)/sngl(reject+1)

RETURN
END SUBROUTINE EXPLORE2

!=======================================================================

SUBROUTINE EXPLORE(obj,dat,logLstar)
!=======================================================================
USE NEST_COM
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER
INTEGER(KIND=RP):: accept, reject,i,j
TYPE (objstruc) :: obj,try
TYPE (datastruc):: dat
REAL(KIND=RP)   :: logLstar
REAL(KIND=RP)   :: ran_uni(NFP)
REAL(KIND=RP)   :: step1

accept = 0      ! # MCMC acceptances
reject = 0      ! # MCMC rejections

DO j = 1,JMAX  ! MCMC counter (pre-judged # steps)
   !
   ! Trial object
   !
   !
   !  Uniform 
   !
!   CALL RANDOM_NUMBER(ran_uni)
!   DO i = 1,NFP
!      try%prior(i) = obj%prior(i) + step1 * (2._RP*ran_uni(i) - 1._RP)  ! |move| < step
!      try%prior(i) = try%prior(i) - FLOOR(try%prior(i))                ! wraparound to stay within (0,1)
!      try%par(i)   = minlim(i) + maxpert(i) * try%prior(i)             ! map to x
!   ENDDO
   !
   !  Gaussian
   !
   DO i = 1,NFP  
      CALL GASDEVJ(ran_uni(i))
      try%prior(i) = obj%prior(i) + step1 * ran_uni(i)      ! |move| < step
      try%prior(i) = try%prior(i) - FLOOR(try%prior(i))    ! wraparound to stay within (0,1)
      try%par(i)   = minlim(i) + maxpert(i) * try%prior(i) ! map to x
   ENDDO   

   CALL LOGLHOOD(try,dat)
!   PRINT*,sngl(try%logL),sngl(logLstar),sngl(obj%logL)

   !
   ! Accept if and only if within hard likelihood constraint
   !
   IF( try%logL > logLstar )THEN
      obj = try
      accept = accept + 1
   ELSE
      reject = reject + 1
   ENDIF

   !
   ! Refine step-size to let acceptance ratio converge around 50%
   !
   IF( accept > reject ) step1 = step1 * exp(1.0_RP / accept);
   IF( accept < reject ) step1 = step1 / exp(1.0_RP / reject);

   !
   ! Exit if already accepted lots of models
   !
   IF(accept > CEILING(dble(JMAX)/2._RP)) EXIT

ENDDO
!PRINT*,'L= ',sngl(obj%logL),sngl(logLstar),sngl(step),sngl(accept+1)/sngl(reject+1)

RETURN
END SUBROUTINE EXPLORE

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
ALLOCATE(tankd(nfrq,NANG),tankdtmp(nfrq,NANG),thtmp(nfrq,NANG),ktmp(nfrq,NANG))

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
      thtmp(i,:) = sin(th(NLAY,:))
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
INTEGER :: i,j
REAL(KIND=RP),   DIMENSION(IMAX,NFP+5):: sample
REAL(KIND=RP) :: mean(NFP),tmp(NFP),std(NFP),wt,logZ

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
INTEGER      :: n,np,i,j,k
REAL(KIND=RP) :: A(np,np),A2(np,np),L(np,np),p(n),summ

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
! This is the end my fiend...
! EOF
