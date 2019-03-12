!
! Remove deadend from experiment list
!

subroutine remove_deadend(experiment,irm)
  
  use topolink_data
  implicit none
  integer :: i, irm, j
  type(experiment_data) :: experiment
  type(observed_deadend), allocatable :: newset(:)

  allocate( newset(experiment%ndeadends-1) )

  j = 0
  do i = 1, experiment%ndeadends
    if ( i == irm ) cycle
    j = j + 1
    newset(j) = experiment%deadend(i)
  end do

  deallocate(experiment%deadend)
  experiment%ndeadends = experiment%ndeadends-1
  allocate( experiment%deadend(experiment%ndeadends) )
  experiment%deadend  = newset

  deallocate(newset)

end subroutine remove_deadend


