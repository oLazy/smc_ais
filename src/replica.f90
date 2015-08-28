!==============================================================================
!
!  Reversible Jump MCMC Sampling for relfection coefficient inversion
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, January 23 2010
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: January 25 2010
!
!  Based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
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
   INTEGER(KIND=IB), PARAMETER :: ISPHER     = 1     ! Spherical refl. coeff. 
   INTEGER(KIND=IB), PARAMETER :: IMAP       = 0     ! WRITE REPLICA AND EXIT
   INTEGER(KIND=IB), PARAMETER :: ICOUPLE_CR = 1     ! Constrain discontinuities of rho and v to the same sign
   INTEGER(KIND=IB), PARAMETER :: subsmp     = 1_IB  ! Subsample Sommerfeld integral plane wave part
   INTEGER(KIND=IB), PARAMETER :: ICOOL      = 1     ! Do cooling during burn-in

!!
!!  AUV DATA RJMCMC trial ping
!!
   INTEGER(KIND=IB), PARAMETER :: ICOV = 0    ! 0 = Sample implicit over sigma
                                              ! 1 = Sample over sigma
                                              ! 2 = Use sigma from ping ave
                                              ! 3 = Use sigma from ping ave but scale by factor (free parameter)
   INTEGER(KIND=IB), PARAMETER :: IdB  = 0    ! 1=Carry out computation in dB
   INTEGER(KIND=IB), PARAMETER :: NANG = 32
   INTEGER(KIND=IB), PARAMETER :: NLMX = 10
   INTEGER(KIND=IB), PARAMETER :: NAVEF= 1            ! # freq per band
   REAL(KIND=RP),    PARAMETER :: frbw = 1._RP/15._RP ! Frac. bandwidth for freq ave.

   CHARACTER(len=64) :: filebasefile      = 'filebase.txt'
   INTEGER(KIND=IB),PARAMETER  :: IGA     = 0    ! Type of averaging (0=intensity, 1=Gaussian, 2=1/3 octave intensity)
   INTEGER(KIND=IB), PARAMETER :: NBAND   = 3
   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 4
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP, 2400._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 6
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP, 2400._RP, &
!                                                 2800._RP, 3200._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 9
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/1000._RP, 1200._RP, 2000._RP, 2200._RP, 2400._RP, &
!                                                 2600._RP, 2800._RP, 3000._RP, 3200._RP/)
!   INTEGER(KIND=IB), PARAMETER :: NBAND   = 6
!   REAL(KIND=RP),   DIMENSION(NBAND):: bands = (/2100._RP, 2300._RP, 2500._RP, 2700._RP, 2900._RP, 3100._RP/)
   REAL(KIND=RP) :: FBW = 12.5_RP
!   INTEGER(KIND=IB), PARAMETER :: NSD  = NBAND
   INTEGER(KIND=IB), PARAMETER :: NSD  = 1

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP), DIMENSION(4+NSD):: minlim   = 0._RP
   REAL(KIND=RP), DIMENSION(4+NSD):: maxlim   = 0._RP
   INTEGER(KIND=IB)               :: kmin     = 0       ! Min number of layers
   INTEGER(KIND=IB)               :: kmax     = 0       ! Max number of layers
   REAL(KIND=RP), PARAMETER       :: hmin     = 0.05_RP ! Min allowed layer thickness
   REAL(KIND=RP), DIMENSION(4+NSD):: maxpert  = 0._RP
   REAL(KIND=RP), DIMENSION(4+NSD):: pertsd   = 0._RP 
   REAL(KIND=RP), DIMENSION(4+NSD):: pertsdsc = 40._RP 
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE:: fr,fstep   ! Total frequency array
   REAL(KIND=RP)               :: z_t
   REAL(KIND=RP)               :: cw
   REAL(KIND=RP)               :: rw
   REAL(KIND=RP)               :: hmx
   CHARACTER(len=64) :: filebase
   INTEGER(KIND=IB)  :: filebaselen
   CHARACTER(len=64) :: infile1
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(len=64) :: repfile
   CHARACTER(len=64) :: sdfile
   CHARACTER(len=64) :: samplefile

!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX = (NLMX * 4) + 3 + NSD
   INTEGER(KIND=IB)           :: ioutside = 0
   INTEGER(KIND=IB)           :: ireject = 0, iaccept = 0
   INTEGER(KIND=IB)           :: ireject_bd = 0, iaccept_bd = 0
   INTEGER(KIND=IB)           :: i_bd    ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_zpert ! Z perturb track (0=nothing, 1=birth-death, 2=perturb 1 z)

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: NCHAIN     = 1E8_IB  ! # iterations (max # MCMC steps)
   INTEGER(KIND=IB),PARAMETER    :: ICHAINTHIN = 500_IB  ! Chain thinning interval
   INTEGER(KIND=IB),PARAMETER    :: NKEEP      = 20      ! Number models to keep before writing
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 10      ! Misc parameters in sample (for bookeeping)

!!
!! Annealing burn-in parameters (sets beta schedule)
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: beta                 ! Inverse Temprerature
   INTEGER(KIND=IB),PARAMETER            :: NTEMP    = 5e4       ! Number of beta values
   INTEGER(KIND=IB),PARAMETER            :: NBAL     = 5e0       ! Number of balanching steps at each T
   REAL(KIND=RP)                         :: beta1gl  = 0.0001_RP ! Global inverse T to define annealing schedule (start)
   REAL(KIND=RP)                         :: beta4gl  = 1.00_RP   ! Global inverse T to define annealing schedule (end for 3 exp legs)

!!
!!  Structures for objects and data 
!!
   TYPE :: objstruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NLMX*4+3+NSD) :: par        ! Forward parameters
      REAL(KIND=RP),DIMENSION(NLMX)         :: z          ! Forward parameters
      REAL(KIND=RP),DIMENSION(3)            :: g          ! Acoustic parameters birth-death layer
      REAL(KIND=RP),DIMENSION(3)            :: gp         ! Acoustic parameters birth-death layer perturbed
      REAL(KIND=RP),DIMENSION(3,3)          :: Chat,Chati ! Covariance matrix for perturbing one BD layer
      REAL(KIND=RP)                         :: detChat
      INTEGER(KIND=IB)                      :: k          ! Layer dimension
      INTEGER(KIND=IB)                      :: NFP        ! Number forward parameters
      REAL(KIND=RP)                         :: logL       ! log likelihood
      REAL(KIND=RP)                         :: logPr      ! log Prior probability ratio
      REAL(KIND=RP)                         :: lognorm    ! Data covariance matrices
   END TYPE objstruc

   TYPE :: datastruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Robs   = 0._RP ! Observed data
      REAL(KIND=RP),DIMENSION(NANG)           :: angobs = 0._RP ! Observed angles
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rrep   = 0._RP ! Replica data for trial models
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: res    = 0._RP ! Replica data for trial models
!      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: Rex    = 0._RP ! Index for bad points/data gaps
      REAL(KIND=RP),DIMENSION(NBAND,NANG)     :: sigma  = 0._RP ! Index for bad points/data gaps
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cdi    = 0._RP ! Inverse data covariance matrices
      REAL(KIND=RP),DIMENSION(NANG,NANG,NBAND):: Cd     = 0._RP ! Data covariance matrices
                                                                ! meaningful for simulations)
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
!!  Global variables
!!
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE:: sample  ! Posterior sample
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE  :: buf_save_snd,buf_save_rcv
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE::buffer1
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE   ::buffer2,buffer3

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
   INTEGER(KIND=IB)            :: rank,NTHREAD,ncount,ierr
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

PROGRAM  RJMCMC_PLANE

!=======================================================================
USE MPI
USE RJMCMC_COM
USE NR
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ipar,ifreq,imcmc,ithin,isource,ikeep

TYPE (objstruc)                            :: obj,objmax ! Objects in likelihood box
TYPE (datastruc)                           :: dat        ! Data

REAL(KIND=RP)                              :: ran_uni    ! Uniform random number
REAL(KIND=RP)                              :: ran_nor    ! Normal random number
REAL(KIND=RP)                              :: logLG      ! Global best log likelihood
REAL(KIND=RP)                              :: flo,fhi

REAL(KIND=RP)                              :: LOGPLUS    ! Function: Addition carried out in log space

INTEGER(KIND=IB),DIMENSION(ICHAINTHIN)     :: idxchain   !! Chain thinning array (use random 
                                                         !! perturbation to randomize chain thinning)
INTEGER(KIND=IB):: NTHIN

!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

REAL(KIND=RP)               :: tstart, tend              ! Overall time 
REAL(KIND=RP)               :: tstart2, tend2            ! Time for one forward model computation
REAL(KIND=RP)               :: tstartsnd, tendsnd        ! Communication time
REAL(KIND=RP)               :: tstartcmp, tendcmp, tcmp  ! Forward computation time

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )

OPEN(UNIT=20,FILE=filebasefile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)

infile1        = filebase(1:filebaselen) // '.txt'
logfile        = filebase(1:filebaselen) // '_RJMH.log'
seedfile       = filebase(1:filebaselen) // '_seeds.log'
mapfile        = filebase(1:filebaselen) // '_map.dat'
repfile        = filebase(1:filebaselen) // '_rep.dat'
sdfile         = filebase(1:filebaselen) // '_sigma.txt'
samplefile     = filebase(1:filebaselen) // '_sample.txt'

IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '...running on ',NTHREAD,' cores'

ALLOCATE( sample(NKEEP,NFPMX+NAP) )
ncount = NKEEP*(NFPMX+NAP)
ALLOCATE( buffer1(1+MPI_BSEND_OVERHEAD),buffer2(ncount+MPI_BSEND_OVERHEAD),buffer3(1+MPI_BSEND_OVERHEAD) )
ALLOCATE( buf_save_snd(NKEEP*(NFPMX+NAP)),buf_save_rcv(NKEEP*(NFPMX+NAP)) )
buf_save_snd= 0._RP
buf_save_rcv= 0._RP
sample      = 0._RP
!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
IF(rank == src)WRITE(6,209) 'Loading data from file...',infile1
OPEN(UNIT=20,FILE=infile1,FORM='formatted',STATUS='OLD',ACTION='READ')
!IF(ISPHER == 1)THEN
READ(20,*) z_t
!ENDIF
READ(20,*) cw
READ(20,*) rw
READ(20,*) hmx
DO ifreq = 1,NBAND
    READ(20,*) (dat%Robs(ifreq,j),j=1,NANG)
ENDDO
READ(20,*) (dat%angobs(j),j=1,NANG)
!DO ifreq = 1,NBAND
!    READ(20,*) (dat%Rex(ifreq,j),j=1,NANG)
!ENDDO
CLOSE(20)
IF(ICOV >= 2)THEN
   OPEN(UNIT=20,FILE=sdfile,FORM='formatted',STATUS='OLD',ACTION='READ')
   DO ifreq = 1,NBAND
      READ(20,*) (dat%sigma(ifreq,j),j=1,NANG)
      dat%logdet(ifreq) = SUM(2._RP*LOG(dat%sigma(ifreq,:)))
      dat%lognorm(ifreq) = -0.5_RP*REAL(NANG,KIND=RP)*LOG(2._RP*PI) - 0.5_RP*dat%logdet(ifreq)
   ENDDO
   CLOSE(20)
ENDIF

!!
!!  Prior bounds
!!
!! AUV DATA:
minlim = (/ hmin, 1450._RP, 1.20_RP, 0.001_RP, 0.001_RP /)
maxlim = (/ hmx,  1700._RP, 2.10_RP, 1.000_RP, 1.000_RP /)
!minlim = (/ hmin, 1450._RP, 1.20_RP, 0.001_RP, 0.001_RP, 0.001_RP, 0.001_RP /)
!maxlim = (/ hmx,  1750._RP, 2.20_RP, 1.000_RP, 6.000_RP, 6.000_RP, 6.000_RP /)

kmin = 1
kmax = NLMX
maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc
!!------------------------------------------------------------------------
!!
!!  Print sampling parameters to screen for logging
!!
IF(rank == src)THEN
   WRITE(6,210) 'ICOV                 :   ',ICOV
   WRITE(6,210) 'Number of angles     :   ',NANG
   WRITE(6,210) 'Number of frequencies:   ',NBAND
   CALL FLUSH(6)
   WRITE(6,209) 'Sample file:             ',samplefile
   WRITE(6,*) ''
   WRITE(6,*) 'NFPMX    = ',NFPMX
   WRITE(6,*) 'ang(min) = ',dat%angobs(1)
   WRITE(6,*) 'ang(max) = ',dat%angobs(NANG)
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

!!
!! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)
!!
!CALL RANDOM_SEED
CALL PARALLEL_SEED()
!!
!! Make cooling schedule
!!
IF(ICOOL == 1) THEN
   IF(rank == src)THEN
      WRITE(6,*) ''
      WRITE(6,*) '  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
      WRITE(6,*) 'Cooling turned on:'
      WRITE(6,*) 'NTEMP = ',NTEMP
      WRITE(6,*) 'NBAL  = ',NBAL
      WRITE(6,*) 'beta1 = ',beta1gl
      WRITE(6,*) 'beta4 = ',beta4gl
      WRITE(6,*) ''
   ENDIF
   ALLOCATE( beta(NTEMP) )
   CALL MAKE_BETA(beta1gl,beta4gl,NTEMP)
ENDIF

!!
!! some factors for the acceptance ratio in birth/death case:
!!
obj%Chat    = 0._RP
obj%Chati   = 0._RP
obj%detChat = 1._RP
DO ipar = 1,3
   obj%Chat(ipar,ipar)  = pertsd(ipar+1)**2._RP
   obj%Chati(ipar,ipar) = 1._RP/pertsd(ipar+1)**2._RP
   obj%detChat = obj%detChat * obj%Chat(ipar,ipar)
ENDDO
obj%logPr = SUM(LOG(maxlim(2:4)-minlim(2:4))) ! prior ratio for change in dimension of 1
obj%lognorm = LOG(SQRT((2._RP*PI)**3._RP*obj%detChat))

!!
!! Test model (first element gives number of layers):
!!
obj%k   = 3

obj%NFP = (obj%k * 4) + 3 + NSD
!ALLOCATE( obj%par(NFPMX),obj%z(kmax) )
obj%par = 0._RP
obj%z   = 0._RP

!! SIMULATION:
obj%par(1:obj%NFP) = (/1.00_RP, 1480._RP, 1.30_RP, 0.01_RP,&
                       0.50_RP, 1600._RP, 1.75_RP, 0.20_RP,&
                       1.60_RP, 1550._RP, 1.50_RP, 0.20_RP,&
                                1680._RP, 1.95_RP, 0.10_RP,0.033_RP/)
!obj%par(1:obj%NFP) = (/1.76997149_RP,  1561.0642_RP,  1.71335923_RP,  0.977785179_RP,&
!                       5.07706022_RP,  1527.2039_RP,  1.46706705_RP,  0.758054738_RP,&
!                                   1644.7873_RP,  2.07418183_RP,  0.0419632979_RP,  0.0688066987_RP /)
!CALL PRIOR(obj,dat)
CALL CALCZ(obj)

ALLOCATE( icount(NTHREAD) )
icount = 0

tstart = MPI_WTIME()

!!
!! 
!!
ALLOCATE( fr(NBAND*NAVEF),fstep(NBAND) )
!bands = 304._RP
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

IF(IMAP == 1)THEN
   IF(rank == src)WRITE(6,*) 'IMAP activated, exiting after computing replica for MAP.'
   tstart2 = MPI_WTIME()
   CALL LOGLHOOD(obj,dat)
   tend2 = MPI_WTIME()
   IF(rank == src)CALL SAVEREPLICA(obj,dat,repfile)
   IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
   STOP
ELSE
   IF(rank == src)WRITE(6,*) 'Starting model:'
   IF(rank == src)CALL PRINTPAR(obj)
   tstart2 = MPI_WTIME()
   CALL LOGLHOOD(obj,dat)
   tend2 = MPI_WTIME()
   IF(rank == src)WRITE(6,*) 'logL = ',obj%logL
   IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
ENDIF
CALL FLUSH(6)


!CALL LOGLHOOD(obj,dat)
icount(rank+1) = icount(rank+1) + 1

!!
!! Make simulated data:
!!
!DO i = 1,NBAND
!   DO j = 1,NANG
!      CALL GASDEVJ(ran_nor)
!      dat%Rrep(i,j) = dat%Rrep(i,j)+ran_nor*obj%sd
!   ENDDO
!   WRITE(110,201) dat%Rrep(i,:)
!ENDDO
!WRITE(110,201) dat%angobs

!------------------------------------------------------------------------
!
!          ************ RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
IF(rank == src)THEN
   WRITE(6,*) 'Starting RJMCMC sampling...'
   !!
   !!  File to save posterior sample
   !!
   OPEN(UNIT=40,FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
   OPEN(UNIT=44,FILE=logfile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN

DO imcmc = 1,NCHAIN

   tstartsnd = MPI_WTIME() 
   CALL SAVESAMPLE(logLG,isource,tcmp)
   tendsnd = MPI_WTIME()
   IF(MOD(imcmc,30+1)==0)THEN
      WRITE(*,*) ''
      WRITE(6,203) '   imcmc,           logL,          logPr,        lognorm,      k,       iacc/irej, iacc_bd,  time(send), source,  time(comp)'
      WRITE(6,203) '----------------------------------------------------------------------------------------------------------------------------'
   ENDIF
!   sample(ikeep,:) =  (/ obj%logL,obj%logPr,REAL(obj%k,RP),obj%par,& 
!                         REAL(iaccept,RP)/REAL(ireject,RP),REAL(iaccept_bd,RP),&
!                         REAL(ireject_bd,RP),REAL(i_bd,RP),REAL(i_zpert,RP) /)
   WRITE(6,202) imcmc,sample(1,(/1,2,3/)),INT(sample(1,4)),sample(1,5+NFPMX),INT(sample(1,6+NFPMX)),tendsnd-tstartsnd,isource,tcmp
   WRITE(44,202) imcmc,sample(1,(/1,2,3/)),INT(sample(1,4)),sample(1,5+NFPMX),INT(sample(1,6+NFPMX)),tendsnd-tstartsnd,isource,tcmp
   CALL FLUSH(6)
   CALL FLUSH(44)

ENDDO
IF(rank == src)THEN
   CLOSE(40)
ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    SLAVE PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE
IF(ICOOL == 1)THEN
!   IF(rank==2) WRITE(*,211) beta
   objmax = obj
   DO imcmc = 1,NTEMP
      DO ithin = 1,NBAL
         CALL EXPLORE_MH(obj,dat,beta(imcmc))
         IF(obj%logL > objmax%logL) objmax=obj
      ENDDO
      IF(MOD(imcmc,500) == 0)WRITE(*,213) rank,imcmc,'  logL=',obj%logL,'  k=',obj%k,'  objmax_logL=',objmax%logL,'  objmax_k=',objmax%k
   ENDDO
   obj = objmax
   DO ithin = 1,500
      CALL EXPLORE_MH(obj,dat,1._RP)
   ENDDO
   WRITE(*,212)rank,'logL=',obj%logL,'k=',obj%k
ENDIF
ikeep   = 1_IB
idxchain = 0
DO imcmc = 1,NCHAIN

   IF(ikeep == 1)tstartcmp = MPI_WTIME()
   !! Apply radom chain thinning between ICHAINTHIN and 2*ICHAINTHIN steps
!   idxchain = RANDPERM(ICHAINTHIN)
!   NTHIN = idxchain(1)+ICHAINTHIN
   NTHIN = ICHAINTHIN
   DO ithin = 1,NTHIN
      CALL EXPLORE_MH(obj,dat,1._RP)
   ENDDO
   sample(ikeep,:) =  (/ obj%logL, obj%logPr, obj%lognorm, REAL(obj%k,RP), & ! 4 parameters
                         obj%par , & 
                         REAL(iaccept,RP)/REAL(ireject,RP),REAL(iaccept_bd,RP),&
                         REAL(ireject_bd,RP),REAL(i_bd,RP),REAL(i_zpert,RP),REAL(rank,RP) /)

   IF(ikeep == NKEEP)THEN
      tendcmp = MPI_WTIME()
      tcmp = tendcmp-tstartcmp
      CALL SAVESAMPLE(logLG,isource,tcmp)
      ikeep  = 0
      sample = 0._RP
!      IF(rank == 2)PRINT*,iaccept_bd
   ENDIF

   ikeep  = ikeep + 1_IB
ENDDO

ENDIF !! MPI ENDIF

201 FORMAT(200F12.4)
202 FORMAT(I8,3F16.6,I8,F16.6,I8,1F13.3,I8,1F13.3)
203 FORMAT(A124)
!204 FORMAT(3I12,2F16.6,I6,44F16.6,4I12)
205 FORMAT(10F10.4)
209 FORMAT(A26,A40)
210 FORMAT(A26,I4)
211 FORMAT(20000ES12.4)
212 FORMAT(I4,A8,F12.4,A5,I3)
213 FORMAT(I4,I6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
CALL MPI_FINALIZE( ierr ) 

END PROGRAM RJMCMC_PLANE

!=======================================================================

SUBROUTINE LOGLHOOD(obj,dat)
!=======================================================================
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB) :: ifreq,ifr,iang,ncra,ipar,iidx
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP), DIMENSION(obj%NFP)   :: m_inr
REAL(KIND=RP),DIMENSION(NAVEF,obj%k+1):: vp,alf1dB
REAL(KIND=RP), DIMENSION(NAVEF,NANG):: Rpltry
REAL(KIND=RP), DIMENSION(NBAND)     :: Etmp
INTEGER(KIND=IB),DIMENSION((kmax+1)*3+NSD) :: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+NSD) :: idxc,idxa
REAL(KIND=RP)                              :: fref,cref,alfrdB,y

!!
!!  Compute plane wave refl. coeff. (band averaged)
!!

DO ifreq = 1,NBAND
   IF(ISPHER == 0)THEN
      IF(NAVEF > 1)THEN
         CALL REF_NLAY3(dat%angobs,obj%par(1:obj%NFP-NSD),fr((ifreq-1)*NAVEF+1:ifreq*NAVEF),&
                        Rpltry,cw,rw,NAVEF,NANG,obj%NFP-NSD)
      ELSE
         CALL REF_NLAY3(dat%angobs,obj%par(1:obj%NFP-NSD),fr(ifreq),&
                        Rpltry,cw,rw,NAVEF,NANG,obj%NFP-NSD)
      ENDIF
   ELSE
      IF(NAVEF > 1)THEN
         CALL SPH_REF_NLAY(dat%angobs,obj%par(1:obj%NFP-NSD),fr((ifreq-1)*NAVEF+1:ifreq*NAVEF),&
                           Rpltry,cw,rw,z_t,NAVEF,NANG,obj%NFP-NSD,ifreq)
      ELSE
         CALL SPH_REF_NLAY(dat%angobs,obj%par(1:obj%NFP-NSD),fr(ifreq),&
                           Rpltry,cw,rw,z_t,NAVEF,NANG,obj%NFP-NSD,ifreq)
      ENDIF
   ENDIF

!   PRINT*,fr
!   DO ifr=1,NAVEF
!      WRITE(77,207) Rpltry(ifr,:)
!   ENDDO
!   STOP

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

!!
!!  Compute log likelihood
!!
dat%res = (dat%Robs-dat%Rrep)!*dat%Rex
IF(ICOV == 0)THEN
   !!
   !! Sample over sigma (one for all freqs)
   !!
   Etmp = 0._RP
   DO ifr = 1,NBAND
      Etmp(ifr) = REAL(NANG,RP)/2._RP*LOG(SUM(dat%res(:,ifreq)*dat%res(:,ifreq)))
   ENDDO
   Etmp = -Etmp
ELSEIF(ICOV == 1)THEN
   !!
   !! Sample over sigma (one for all freqs)
   !!
   DO ifreq = 1,NBAND
!      Etmp(ifreq) = LOG(1._RP/(2._RP*PI)**(REAL(NANG,RP)/2._RP)) &
!          -(SUM(dat%res(ifreq,:)**2._RP)/(2._RP*obj%par(obj%NFP-NSD+ifreq)**2._RP)&
!          +REAL(NANG,RP)*LOG(obj%par(obj%NFP-NSD+ifreq)))
      Etmp(ifreq) = LOG(1._RP/(2._RP*PI)**(REAL(NANG,RP)/2._RP)) &
          -(SUM(dat%res(ifreq,:)**2._RP)/(2._RP*obj%par(obj%NFP)**2._RP)&
          +REAL(NANG,RP)*LOG(obj%par(obj%NFP)))
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
      Etmp(ifreq) = LOG(1._RP/(2._RP*PI)**(REAL(NANG,RP)/2._RP)) &
          -(SUM(dat%res(ifreq,:)**2._RP/(2._RP*(obj%par(obj%NFP-NSD+ifreq)*dat%sigma(ifreq,:))**2._RP))&
          +REAL(NANG,RP)*LOG(obj%par(obj%NFP-NSD+ifreq))+SUM(LOG(dat%sigma(ifreq,:))))
   ENDDO
ENDIF
obj%logL = SUM(Etmp)

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
INTEGER(KIND=IB) :: iwhich,NFPnew,idel,ipar,ilay,idxz
INTEGER(KIND=IB) :: ncra
TYPE(objstruc) :: obj,objnew
TYPE(datastruc):: dat
REAL(KIND=RP)  :: logPLratio,logPQ,ran_uni,ran_uni_BD, ran_unik
REAL(KIND=RP)  :: znew,beta_mh
INTEGER(KIND=IB),DIMENSION(NFPMX) :: idxrand
INTEGER(KIND=IB),DIMENSION((kmax+1)*3+NSD) :: idxcra
REAL(KIND=RP),DIMENSION(obj%k) :: ztmp

!ALLOCATE( objnew%par(NFPMX),objnew%z(kmax) )

objnew = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)

!! Do normal MCMC with 0.5 probability
IF(ran_uni_BD <= 0.5_RP)THEN
!! Do nothing to k and z
   i_zpert = 0
   i_bd = 0
ENDIF

!! Do BIRTH-DEATH MCMC with 0.25 probability
IF((ran_uni_BD > 0.5_RP) .AND. (ran_uni_BD <= 0.75))THEN 
   i_zpert = 1
   !! Perturbing k:
   CALL RANDOM_NUMBER(ran_unik)
   i_bd = 0
   IF(obj%k == kmax)THEN  !! If k == kmax, no birth allowed
      IF(ran_unik>=0.5_RP) i_bd = 2
   ELSEIF(obj%k == kmin)THEN  !! If k == kmin, no death allowed
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
!! Perturb z with 0.25 probability
IF((ran_uni_BD > 0.75_RP))THEN
   i_zpert = 2
   i_bd = 0
   idxrand = 0
   idxrand(1:obj%k) = RANDPERM(obj%k)
   idxz = idxrand(1)
   iwhich = (idxz-1) * 4 + 1
   DO
      CALL PROPOSALZ(obj,objnew,dat,idxz,iwhich)
      IF(objnew%z(objnew%k) > maxlim(1))CYCLE
      ztmp = objnew%z(1:objnew%k) - objnew%z(idxz)
      IF(hmin>MINVAL(ABS(ztmp)).AND.(hmin>objnew%z(1)))CYCLE
      EXIT
   ENDDO
ENDIF

!!
!! Do Metropolis-Hastings
!!
IF(obj%k /= objnew%k)THEN
   IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)
   IF(ioutside == 0)THEN
      CALL LOGLHOOD(objnew,dat)
      IF(i_bd == 1)THEN
         !! BIRTH ACCEPTANCE RATIO
         logPQ = -objnew%logPr + objnew%lognorm + &
                 (0.5_RP*DOT_PRODUCT(MATMUL((objnew%gp-objnew%g),objnew%Chati),(objnew%gp-objnew%g)))
         logPLratio =  logPQ + (objnew%logL - obj%logL)*beta_mh
      ELSEIF(i_bd == 2)THEN
         !! DEATH ACCEPTANCE RATIO
         logPQ = objnew%logPr - objnew%lognorm - &
                 (0.5_RP*DOT_PRODUCT(MATMUL((objnew%gp-objnew%g),objnew%Chati),(objnew%gp-objnew%g)))
         logPLratio = logPQ + (objnew%logL - obj%logL)*beta_mh
      ENDIF
      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
         ireject_bd = ireject_bd + 1
      ELSE
         obj = objnew
         iaccept_bd = iaccept_bd + 1
      ENDIF 
   ELSE
      objnew = obj
      ireject_bd = ireject_bd + 1
      ioutside = 0
   ENDIF
ELSE
   !! Perturb only c, rho, alpha, and sigma
   IF(ICOV == 0)THEN
      CALL GETIDXCRA(obj,idxcra,ncra)
   ELSEIF(ICOV == 1)THEN
      CALL GETIDXCRA(obj,idxcra,ncra)
   ELSEIF(ICOV == 2)THEN
      CALL GETIDXCRA2(obj,idxcra,ncra)
   ELSEIF(ICOV == 3)THEN
      CALL GETIDXCRA(obj,idxcra,ncra)
   ENDIF
   idxrand = 0
   idxrand(1:ncra) = RANDPERM(ncra)
   DO ipar = 1,ncra
      iwhich = idxcra(idxrand(ipar))
      CALL PROPOSAL(obj,objnew,dat,iwhich)
      IF(ICOUPLE_CR == 1)CALL COUPLE_CR(objnew)

      IF(ioutside == 0)THEN
         CALL LOGLHOOD(objnew,dat)
         logPLratio = (objnew%logL - obj%logL)*beta_mh
         CALL RANDOM_NUMBER(ran_uni)
         IF(ran_uni >= EXP(logPLratio))THEN
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
ENDIF

!DEALLOCATE( objnew%par,objnew%z )
END SUBROUTINE EXPLORE_MH
!=======================================================================

SUBROUTINE PROPOSALZ(obj,objnew,dat,idxz,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,iwhich,idxz,iloop
TYPE(objstruc) :: obj,objnew
TYPE(datastruc):: dat
REAL(KIND=RP)  :: ran_uni, ran_nor

!! Perturb only one z with Gaussian proposal
CALL GASDEVJ(ran_nor)
objnew%z(idxz) = obj%z(idxz) + pertsd(1)*ran_nor
CALL CALCPAR(objnew)

!! Check if outside prior, if so, reject
IF(((objnew%z(idxz) - minlim(1)) < 0._RP).OR.((maxlim(1) - objnew%z(idxz)) < 0._RP))ioutside = 1
!! Check that z is not perturbed into another layer; reject if it happened
IF(idxz ==1)THEN
   IF(objnew%z(idxz) > obj%z(idxz+1))ioutside = 1
ELSEIF(idxz == objnew%k)THEN
   IF(objnew%z(idxz) < objnew%z(idxz-1))ioutside = 1
ELSE
   IF((objnew%z(idxz) < objnew%z(idxz-1)).OR.(objnew%z(idxz) > obj%z(idxz+1)))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSALZ
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,dat,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,iwhich,ipar,iloop,ilay
TYPE(objstruc) :: obj,objnew
TYPE(datastruc):: dat
REAL(KIND=RP)  :: ran_uni, ran_nor

ilay = CEILING(REAL(iwhich,RP)/4._RP)
ipar = iwhich-(ilay-1)*4
IF(ilay > obj%k) ipar = ipar + 1
!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%par(iwhich) = obj%par(iwhich) + pertsd(ipar)*ran_nor
IF(((objnew%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - objnew%par(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL
!=======================================================================

SUBROUTINE GETIDXCRA(obj,idxcra,ncra)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,ilay,isd,ncra
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*3+NSD) :: idxcra

idxcra = 0
ipar = 1
DO ilay = 1,obj%k
   idxcra(ipar) = (ilay-1)*4+2
   ipar = ipar + 1
   idxcra(ipar) = (ilay-1)*4+3
   ipar = ipar + 1
   idxcra(ipar) = (ilay-1)*4+4
   ipar = ipar + 1
ENDDO
idxcra(ipar) = obj%k*4+1
ipar = ipar + 1
idxcra(ipar) = obj%k*4+2
ipar = ipar + 1
idxcra(ipar) = obj%k*4+3
DO isd = 1,NSD
   ipar = ipar + 1
   idxcra(ipar) = obj%k*4+3+isd !! standard deviation
ENDDO

ncra = (obj%k+1)*3+NSD

END SUBROUTINE GETIDXCRA
!=======================================================================

SUBROUTINE GETIDXCRA2(obj,idxcra,ncra)
!=======================================================================
!!
!! This gives index without standard deviation
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,ilay,ncra
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*3) :: idxcra

idxcra = 0
ipar = 1
DO ilay = 1,obj%k
   idxcra(ipar) = (ilay-1)*4+2
   ipar = ipar + 1
   idxcra(ipar) = (ilay-1)*4+3
   ipar = ipar + 1
   idxcra(ipar) = (ilay-1)*4+4
   ipar = ipar + 1
ENDDO
idxcra(ipar) = obj%k*4+1
ipar = ipar + 1
idxcra(ipar) = obj%k*4+2
ipar = ipar + 1
idxcra(ipar) = obj%k*4+3
!ipar = ipar + 1
!idxcra(ipar) = obj%k*4+4 !! standard deviation

ncra = (obj%k+1)*3

END SUBROUTINE GETIDXCRA2
!=======================================================================

SUBROUTINE CALCPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxh
REAL(KIND=IB),DIMENSION(obj%k) :: h

idxh = 0
idxh = (/1:obj%k*4:4/)
h(1) = obj%z(1)
DO i = 2,obj%k
   h(i) = obj%z(i)-obj%z(i-1)
ENDDO
obj%par(idxh) = h

END SUBROUTINE CALCPAR
!=======================================================================

SUBROUTINE CALCZ(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxz

idxz = 0
idxz = (/1:obj%k*4:4/)
obj%z = 0._RP
!obj%z(1:obj%k) = obj%par(idxz)
DO i = 1,obj%k
   obj%z(i) = SUM(obj%par(idxz(1:i)))
ENDDO

END SUBROUTINE CALCZ
!=======================================================================

SUBROUTINE DEATH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: idel
TYPE(objstruc)   :: obj,objnew
TYPE(datastruc):: dat
INTEGER(KIND=IB),DIMENSION(NFPMX) :: idxdeath
REAL(KIND=RP)  :: ran_uni
REAL(KIND=RP),DIMENSION(3) :: cra_ave

objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * 4) + 3 + NSD
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
   cra_ave = (obj%par((idel-1)*4+2:(idel-1)*4+4) + obj%par(idel*4+2:idel*4+4))/2._RP
   objnew%par(1:objnew%NFP) = obj%par(5:obj%NFP)
   objnew%par(2:4) = cra_ave
   objnew%par(1) = obj%par(1) + obj%par(5)
   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par(2:4)
   objnew%gp = objnew%par(2:4)
ELSEIF(idel >= obj%k)THEN
   cra_ave = (obj%par(obj%NFP-NSD-5:obj%NFP-NSD-3) + obj%par(obj%NFP-NSD-2:obj%NFP-NSD))/2._RP
   objnew%par(1:objnew%NFP) = (/ obj%par(1:obj%NFP-NSD-7),obj%par(obj%NFP-NSD-5:obj%NFP-NSD-3),obj%par(obj%NFP-NSD+1:obj%NFP) /)
   objnew%par(objnew%NFP-NSD-2:objnew%NFP-NSD) = cra_ave
   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par(obj%NFP-NSD-2:obj%NFP-NSD)
   objnew%gp = objnew%par(objnew%NFP-NSD-2:objnew%NFP-NSD)
ELSE
   cra_ave = (obj%par((idel-1)*4+2:(idel-1)*4+4) + obj%par(idel*4+2:idel*4+4))/2._RP
   objnew%par(1:objnew%NFP) = (/ obj%par(1:(idel-1)*4), obj%par(idel*4+1:obj%NFP) /)
   objnew%par((idel-1)*4+2:(idel-1)*4+4) = cra_ave
   objnew%par((idel-1)*4+1) = obj%par((idel-1)*4+1) + obj%par((idel-1)*4+5)

   !! This records the perturbation for the bd acceptance ratio
   objnew%g  = obj%par((idel-1)*4+2:(idel-1)*4+4)
   objnew%gp = objnew%par((idel-1)*4+2:(idel-1)*4+4)
ENDIF
CALL CALCZ(objnew)

END SUBROUTINE DEATH
!=======================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i
TYPE(objstruc) :: obj

DO i=1,obj%k
   WRITE(6,201) obj%par((i-1)*4+1:i*4)
ENDDO
WRITE(6,202) '            ',obj%par(obj%k*4+1:obj%k*4+3+NSD)

201 FORMAT(4F12.4)
202 FORMAT(A12,8F12.4)
END SUBROUTINE PRINTPAR
!=======================================================================

SUBROUTINE BIRTH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i, iznew, ipert, iloop, iwhich, ipar
TYPE(objstruc) :: obj,objnew,objnew2
REAL(KIND=RP)  :: znew,hnew1,hnew2,ran_uni,ran_nor
REAL(KIND=RP),DIMENSION(obj%k) :: ztmp

!ALLOCATE( objnew2%par(NFPMX),objnew2%z(kmax) )
objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * 4) + 3 + NSD
!!
!! Draw new z until new layer > hmin
!!
DO
   CALL RANDOM_NUMBER(ran_uni)
   znew = maxpert(1)*ran_uni
   ztmp = obj%z(1:obj%k) - znew
   IF(hmin<MINVAL(ABS(ztmp)).AND.(hmin<znew)) EXIT
ENDDO

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
   hnew1 = znew
   hnew2 = obj%z(iznew)-znew
   objnew%par(1:objnew%NFP) = (/ obj%par(1:4),obj%par /)
   objnew%par(1)            = hnew1
   objnew%par(5)            = hnew2
ELSEIF(iznew > obj%k)THEN
   hnew1 = znew - obj%z(iznew-1)
   objnew%par(1:objnew%NFP) = (/ obj%par(1:obj%NFP-NSD-3),hnew1, &
                                 obj%par(obj%NFP-NSD-2:obj%NFP-NSD), &
                                 obj%par(obj%NFP-NSD-2:obj%NFP) /)
   objnew%par(objnew%NFP-NSD-6) = hnew1
ELSE
   hnew1 = znew - obj%z(iznew-1)
   hnew2 = obj%z(iznew)-znew
   objnew%par(1:objnew%NFP) = (/ obj%par(1:(iznew)*4), &
                                 obj%par(((iznew-1)*4)+1:iznew*4), &
                                 obj%par(iznew*4+1:obj%NFP) /)
   objnew%par(((iznew-1)*4)+1) = hnew1
   objnew%par(((iznew)*4)+1)   = hnew2
ENDIF
CALL CALCZ(objnew)

!!
!! Pick one of the twe new layers at random and perturb:
!!
CALL RANDOM_NUMBER(ran_uni)
IF(ran_uni <= 0.5_RP)THEN
   ipert = iznew
ELSE
   ipert = iznew+1
ENDIF

!! Perturb only one parameter (determined in EXPLORE_MH call)
objnew2 = objnew
DO ipar = 2,4
   iloop = 0
   iwhich = (ipert-1)*4+ipar
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

!DEALLOCATE( objnew2%par,objnew2%z )
END SUBROUTINE BIRTH
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

DO i=1,obj%NFP-1
   ilay = CEILING(REAL(i,RP)/4._RP)
   ipar = i-(ilay-1)*4
   IF(ilay > obj%k) ipar = ipar + 1
   CALL RANDOM_NUMBER(ran_uni)
   obj%par(i) = minlim(ipar) + maxpert(ipar) * ran_uni
ENDDO
CALL RANDOM_NUMBER(ran_uni)
obj%par(obj%NFP) = minlim(5) + maxpert(5) * ran_uni
CALL CALCZ(obj)
CALL LOGLHOOD(obj,dat)

RETURN
END SUBROUTINE PRIOR
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
INTEGER(KIND=IB) :: kf,lf,i,j,l,idx,Ni,Nbuf
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
!      IF(MOD(Nr(lf),subsmp) /= 0._RP)THEN
!         i = i+1
!         DO j=1,MOD(Nr(lf),subsmp)
!!            ref2(l) = ref(i)+(ref(i+1)-ref(i))*REAL(j-1,RP)/REAL(MOD(N(lf)-Ni,subsmp),RP)
!            ref2(l) = ref(i)+(ref(i+1)-ref(i))*REAL(j-1,RP)/REAL(MOD(Nr(lf),subsmp),RP)
!            l = l + 1
!         ENDDO
!      ENDIF
      ref2(N(lf)-Ni:N(lf)) = ref(N_tmp-Ni:N_tmp)

!      DO i = 1,N_tmp
!         WRITE(68,*) 90._RP-REAL(rTh_tmp(i))*180._RP/PI,REAL(ref(i))
!         WRITE(69,*) 90._RP-AIMAG(rTh_tmp(i))*180._RP/PI,AIMAG(ref(i))
!         WRITE(70,*) 90._RP-ABS(rTh_tmp(i))*180._RP/PI,ABS(ref(i))
!      ENDDO
!
!
!      DO i = 1,N(lf)
!         WRITE(78,*) 90._RP-REAL(rTh(i,lf))*180._RP/PI,REAL(ref2(i))
!         WRITE(79,*) 90._RP-AIMAG(rTh(i,lf))*180._RP/PI,AIMAG(ref2(i))
!         WRITE(80,*) 90._RP-ABS(rTh(i,lf))*180._RP/PI,ABS(ref2(i))
!      ENDDO
!
!      !!!!!
!      ref2 = 0._RP
!      CALL REF_NLAY4(90._RP-rTh(:,lf)*180._RP/PI,m_rg,freq(kf),ref2,cw,rw,1,N(lf),NFP)
!      DO i = 1,N(lf)
!         WRITE(88,*) 90._RP-REAL(rTh(i,lf))*180._RP/PI,REAL(ref2(i))
!         WRITE(89,*) 90._RP-AIMAG(rTh(i,lf))*180._RP/PI,AIMAG(ref2(i))
!         WRITE(90,*) 90._RP-ABS(rTh(i,lf))*180._RP/PI,ABS(ref2(i))
!      ENDDO
!
!      STOP
!      !!!!!
   ELSE
      ref2 = 0._RP
      CALL REF_NLAY4(90._RP-rTh(1:N(lf),lf)*180._RP/PI,m_rg,freq(kf),ref2,cw,rw,1,N(lf),NFP)
   ENDIF

   rRP = drTh(lf)/3._RP*MATMUL(ref2(1:Nr(lf)),btR(1:Nr(lf),:,lf))
   iRP = diTh(lf)/3._RP*MATMUL(ref2(Nr(lf):N(lf)),btR(Nr(lf):N(lf),:,lf))

   ref_sph(kf,:) = ABS(CMPLX(0._RP,1._RP,RP)*(rRP+iRP)/RPs(:,lf))
   DEALLOCATE( ref,ref2,rTh_tmp )
ENDDO
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SPH_REF_NLAY
!=======================================================================

SUBROUTINE SPH_REF_NLAY_BKUP(thd,m_rg,freq,ref_sph,cw,rw,z_t,nfrq,NANG,NFP,idx)
!
!     Compute spherical wave reflection coeff
!     CANNOT HANDLE HALFSPACE AT THIS POINT
!     Based on code by Charles W. Holland and John Camin ARL PSU
!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: kf,lf,i,j,l,idx
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
   IF(MOD(N(lf),subsmp) == 0._RP)THEN
      N_tmp = CEILING(REAL(N(lf)-1,RP)/REAL(subsmp,RP))
   ELSE
      N_tmp = CEILING(REAL(N(lf)-1,RP)/REAL(subsmp,RP))+1
   ENDIF
   ALLOCATE( ref(N_tmp),ref2(N(lf)),rTh_tmp(N_tmp) )
   ref = 0._RP
   ref2= 0._RP
   rTh_tmp = 0._RP
   j = 1
   DO i=1,N(lf),subsmp
      rTh_tmp(j) = rTh(i,lf)
      j = j + 1
   ENDDO
   rTh_tmp(N_tmp) = rTh(N(lf),lf)

   CALL REF_NLAY4(90._RP-rTh_tmp*180._RP/PI,m_rg,freq(kf),ref,cw,rw,1,N_tmp,NFP)

   l = 1
   DO i=1,N_tmp-2
      DO j=1,subsmp
         ref2(l) = ref(i)+(ref(i+1)-ref(i))*REAL(j-1,RP)/REAL(subsmp,RP)
         l = l + 1
      ENDDO
   ENDDO
   IF(MOD(N(lf),subsmp) /= 0._RP)THEN
      i = i+1
      DO j=1,MOD(N(lf),subsmp)
         ref2(l) = ref(i)+(ref(i+1)-ref(i))*REAL(j-1,RP)/REAL(MOD(N(lf),subsmp),RP)
         l = l + 1
      ENDDO
   ENDIF
   ref2(N(lf)) = ref(N_tmp)
   rRP = drTh(lf)/3._RP*MATMUL(ref2(1:Nr(lf)),btR(1:Nr(lf),:,lf))
   iRP = diTh(lf)/3._RP*MATMUL(ref2(Nr(lf):N(lf)),btR(Nr(lf):N(lf),:,lf))

   ref_sph(kf,:) = ABS(CMPLX(0._RP,1._RP,RP)*(rRP+iRP)/RPs(:,lf))
   DEALLOCATE( ref,ref2,rTh_tmp )
ENDDO
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SPH_REF_NLAY_BKUP
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
   STOP
ELSE
    vp=1._RP/(1._RP/cr + Q2*alfr/(2._RP*PI*fr)**y)
    alf1dB= alfrdB*(f1/fr)**(y-1._RP)
ENDIF

RETURN
END SUBROUTINE KK_WATERS_dBmkHz
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
      7717070,      141180,     3783525,     3087889,     4812786,     3028075, &
      3712062,     6316731,      436800,     7957708,     2055697,     1944360, &
      1222992,     7537775,     7769874,     5588112,     7590383,     1426393, &
      1753301,     7681841,     2842400,     4411488,     7304010,      497639, &
      4978920,     5345495,      754842,     7360599,     5776102/)
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

201   FORMAT(50I10)
END SUBROUTINE PARALLEL_SEED
!=======================================================================

SUBROUTINE SAVESAMPLE(logLG,isource,tcmp)
!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB):: i,j,ikeep,isource
REAL(KIND=RP):: logLG,tcmp
REAL(KIND=RP):: LOGPLUS
REAL(KIND=RP):: t1,t2

!!
!!  Sending samples to master
!!
IF(rank == src)THEN
   ikeep = NKEEP
   IF(iconv == 1)THEN
      iconv2 = 1
      iconv3 = iconv3 + 1
   ENDIF
   tag   = MPI_ANY_TAG
   buf_save_rcv = 0._RP

   CALL MPI_RECV(buf_save_rcv, ncount, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
                 tag, MPI_COMM_WORLD, status, ierr )
   isource = status(MPI_SOURCE)
   CALL MPI_RECV(tcmp, 1, MPI_DOUBLE_PRECISION, isource,&
                 tag, MPI_COMM_WORLD, status, ierr )

   isize1 = SIZE(buffer1,1)*IB
   CALL MPI_BUFFER_ATTACH(buffer1,isize1,ierr)
   CALL MPI_BSEND(iconv,1,MPI_INTEGER,isource,rank,MPI_COMM_WORLD,ierr)
   CALL MPI_BUFFER_DETACH(buffer1,isize1,ierr)

   sample = RESHAPE(buf_save_rcv,(/ NKEEP,NFPMX+NAP /))
   !!
   !! Master writes sample to file
   !!
   DO j=1,NKEEP
      WRITE(40,207) sample(j,:)
   ENDDO
   CALL FLUSH(40)
   ikeep = ikeep + NKEEP
   logLG = MAXVAL(sample(:,1))
ELSE
   buf_save_snd = 0._RP
   buf_save_snd = PACK(sample(1:NKEEP,:),.true.)

   isize2 = SIZE(buffer2,1)*RP
   CALL MPI_BUFFER_ATTACH(buffer2,isize2,ierr)
   CALL MPI_BSEND( buf_save_snd, ncount, MPI_DOUBLE_PRECISION, src, &
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_BUFFER_DETACH(buffer2,isize2,ierr)

   isize3 = size(buffer3,1)*RP
   CALL MPI_BUFFER_ATTACH(buffer3,isize3,ierr)
   CALL MPI_BSEND( tcmp, 1, MPI_DOUBLE_PRECISION, src, &
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_BUFFER_DETACH(buffer3,isize3,ierr)

   CALL MPI_RECV(iconv,1,MPI_INTEGER,src,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
ENDIF !  MPI
CALL FLUSH(6)
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SAVESAMPLE
!=======================================================================

SUBROUTINE COUPLE_CR(obj)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ncra
TYPE (objstruc)  :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*3+NSD) :: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+NSD) :: idxc,idxr
REAL(KIND=RP),DIMENSION(kmax):: dc,dr

idxcra= 0
idxc  = 0
idxr  = 0
CALL GETIDXCRA(obj,idxcra,ncra)
idxc = idxcra(1:ncra-NSD:3)
idxr = idxcra(2:ncra-NSD:3)

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
!==============================================================================

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
!=======================================================================

FUNCTION RANDPERM(num)
!==============================================================================
USE data_type, ONLY : IB, RP
IMPLICIT NONE
INTEGER(KIND=IB), INTENT(IN) :: num
INTEGER(KIND=IB) :: number, i, j, k
INTEGER(KIND=IB), DIMENSION(num) :: randperm
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
   randperm(i)=number
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
