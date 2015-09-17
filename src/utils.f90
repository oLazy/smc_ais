module utils
  use data_type
  use rjmcmc_com
  implicit none
contains
  function resample(oldSample, oldSize, newSize, mode, k_mode) result (newSample)
    ! this function resamples keeping the k-distribution.
    ! modes
    ! 'n' : no weight on the particles: all the particles with the same k have the same probability to be drawn for the sample
    ! 'w' : Likelihood weight on the particles: within the set of particles with the same k, the probability of a particle to be drawn is proportional to its likelyhood value

    ! developer note: it is easy to extend this code by adding modes and defining new cdfs in the select case block.
    ! it works as expected on test data
    ! Eric Mandolesi 20 August 2015
    
    implicit none
    integer (kind=ib), intent(in) :: oldSize, newSize
    type(objstruc), dimension(oldSize), intent(in) :: oldSample
    character :: mode
    character(len=5) :: k_mode

    type(objstruc), dimension(:),allocatable :: newSample
    !! dummy arguments

    integer(kind=ib), dimension(oldSize) :: kk
    
    integer(kind=ib) :: i,j,k 
    integer(kind=ib), dimension(:), allocatable :: v
    real(kind=rp), dimension(:,:), allocatable :: findit
    real(kind=rp), dimension(:), allocatable :: cs
    real(kind=rp), dimension(:,:), allocatable :: lkcum
    real(kind=rp) :: sumVr , lknorm
    integer(kind=ib) :: maxK
    real(kind=rp) :: r1, r2
    
    if(.not.allocated(newSample))allocate(newSample(newSize))

    !    print*, 'resample: in module'    
    do i = 1,oldSize
       kk(i) = oldSample(i)%k
    end do
    maxK = maxval(kk)
    !write(*,'(I3)') maxK
    allocate(v(maxK))
    allocate(findit(oldSize,maxK))
    allocate(cs(maxK))
    allocate(lkcum(oldSize,maxK))    
    print*, 'resample: arrays allocated!'
    
    sumVr =  1._rp/oldSize
    v = 0;
    findit = 0._rp
    
    do i = 1,maxK
       do j = 1,oldSize
          if(oldSample(j)%k .eq. i )then
             v(i) = v(i) + 1
             
             select case (mode)
             case ('n','N')
                findit(j,i) = 1._rp
             case ('w','W')
!                print*, 'resample: w case'    
                findit(j,i) = exp(oldSample(j)%logL)
             case('t')
                findit(j,i) = oldSample(j)%logL
             case default
                write(*,*) 'wrong usage function resample in module utils.f90 - unknown mode used. Quitting.'
                stop
             end select
             
          end if
       end do
    end do

    if(sum(v).ne.oldSize)then
       print*, 'error in module utils.f90 - wrong cumulative sum for v(:). Quiting.'
       print*, sum(v)
       stop
    end if

    ! compute the cpd for k
    cs = 0;
    do i = 1,maxK
       if (i.eq.1)then
          cs(i) = v(i)*sumVr
       else
          cs(i) = cs(i-1) + v(i)*sumVr
       end if
    end do

    do i = 1,maxK
       lknorm = sum(findit(1:oldSize,i));
       if (lknorm == 0)then
          lkcum(1:oldSize,i) = 1
       else
          select case (mode)
          case ('t')
             ! test case for the logCumsum function
             lkcum(1:oldSize,i) = logCumsum(findit(1:oldSize,i))
          case default
             lkcum(1:oldSize,i) = cumsum(findit(1:oldSize,i))/lknorm
          end select
       end if
    end do

    
    !!! CRASH AFTER THIS POINT
    do i = 1,newSize
       print*, 'i = ', i
       call random_number(r1)
       print*, 'r1 = ', r1
!       r1 = rand();
       jdo: do j = 1,maxK
          print*, 'j = ', j
          print*, 'cs(j) = ', cs(j)
          if(r1 <= cs(j)) then
!             print*, i+1,' test'
 !            newSample(i+1)%k = j
  !           print*, i+1,' test: passed!'
             print*, 'in the if construct!'
   !          print*, 'i = ', i
    !         print*, 'size(newSample) = ', size(newSample)
             newSample(i)%k = j; !%newSample containst the number of particle per k
             exit jdo
             print*, 'out of jdo cycle'
          end if
       end do jdo
    end do
    !!! CRASH BEFORE THIS POINT
    select case (k_mode)


    case('fixed')
       if (newSize.ne.oldSize) then
          print*, 'error in module utils.f90 - k_mode = fixed'
          print*, 'the two samples have different dimension. Quitting.'
          stop
       end if
       
       i_do1: do i = 1,newSize
          call random_number(r2)
          j = oldSample(i)%k
          k_do1: do k = 1,oldSize
             if(r2 <= lkcum(k,j))then
                newSample(i) = oldSample(k)
                exit k_do1
             end if
          end do k_do1
       end do i_do1
          
       
    case ('kdist')
       print*, 'resample: in kdist'    
       i_do: do i = 1,newSize
          call random_number(r1)
          j_do: do j = 1,maxK
             if(r1 <= cs(j))then
                call random_number(r2)
                k_do: do k = 1,oldSize
                   if(r2 <= lkcum(k,j)) then
                      newSample(i) = oldSample(k);
                      exit k_do
                   end if
                end do k_do
                exit j_do
             end if
          end do j_do
       end do i_do
    end select

    print*, 'resample: before deallocate calls'    
    deallocate(v)
    deallocate(findit)
    deallocate(cs)
    deallocate(lkcum)  
    print*, 'resample: exit'
    return
  end function resample


  
  function cumsum(array,shift)
    ! performs the cumulative sum of the vector "array" shifting each element by a constant "shift"
    use data_type
    !use rjmcmc_com
    implicit none
    real(kind=rp), dimension(:), intent(in) :: array
    real(kind=rp), optional, intent(in) :: shift
    real(kind=rp), dimension(size(array)) :: cumsum
    integer (kind=ib) :: n
    
    integer(kind=ib) :: i
    real(kind=rp) :: temp

    n = size(array)
    
    if (n.eq.0) return
    temp = 0._rp
    if(present(shift)) temp = shift
    cumsum(1) = array(1) + temp
    do i = 2,n
       cumsum(i) = cumsum(i-1) + array (i)
    end do
    return
  end function cumsum

  function logsum(a,b) result(c)
    use data_type
    implicit none
    real(kind=rp), intent(in) :: a,b
    real(kind=rp) :: c

    if(a.gt.b) then
       c = b + log (1 + exp (a - b) );
    else
       c = a + log (1 + exp (b - a) );
    end if
    
  end function logsum

  function logCumsum(array)
    use data_type
    implicit none
    real(kind=rp), dimension(:), intent(in) :: array
    
    real(kind=rp), dimension(size(array)) :: logCumsum
    real(kind=rp), dimension(size(array)) :: t
    integer(kind=ib) :: i,n

    t(1) = exp(log(array(1)))
    do i = 2,size(array)
       t(i) = logsum(t(i-1),array(i))
    end do
    logCumsum = t - t(size(t));
  end function logCumsum

end module utils
