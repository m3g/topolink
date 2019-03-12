!
! Remove observation from experiment list
!

subroutine remove_observed(experiment,irm)
  
  use topolink_data
  implicit none
  integer :: i, irm, j
  type(experiment_data) :: experiment
  type(observed_link), allocatable :: newset(:)

  allocate( newset(experiment%nobs-1) )

  j = 0
  do i = 1, experiment%nobs
    if ( i == irm ) cycle
    j = j + 1
    newset(j) = experiment%observed(i)
  end do

  deallocate(experiment%observed)
  experiment%nobs = experiment%nobs-1
  allocate( experiment%observed(experiment%nobs) )
  experiment%observed = newset

  deallocate( newset )

end subroutine remove_observed


