program test
  use data_type
  use rjmcmc_com
  use utils
  implicit none

  type(objstruc), dimension(10000) :: oldSample
  type(objstruc), dimension(10000) :: newSample
  character(len=8) :: filename = 'test.txt'

  integer(ib) :: i
  integer(ib) :: handle = 1979
  integer :: info = 0
 
  open (handle,file=trim(filename),action='read',status='old',iostat=info)
  if(info.ne.0)then
     print*, 'error opening file ', trim(filename), '.'
     stop
  end if
  
  do i = 1,10000
     read(handle,*) oldSample(i)%k, oldSample(i)%logL
  end do
  
  close(handle)
  
  newSample = resample(oldSample,10000,10000,'w')
  

  open (handle,file='test.out',action='write',status='unknown',iostat=info)
  if(info.ne.0)then
     print*, 'error opening file ', 'test.out', '.'
     stop
  end if
  
  do i = 1,10000
     write(handle,*) newSample(i)%k, newSample(i)%logL
  end do
  
  close(handle)


  newSample = resample(oldSample,10000,10000,'t')
  

  open (handle,file='lc_test.out',action='write',status='unknown',iostat=info)
  if(info.ne.0)then
     print*, 'error opening file ', 'lc_test.out', '.'
     stop
  end if
  
  do i = 1,10000
     write(handle,*) newSample(i)%k, newSample(i)%logL
  end do
  
  close(handle)
end program test
