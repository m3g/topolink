!
! Subroutine that checks the conistency of the link with the 
! experimental observations 
!

subroutine linkconsistency(link,nexp,experiment)

  use topolink_data
  use topolink_operations
  implicit none
  integer :: nexp, iexp, j 
  type(specific_link), intent(inout) :: link
  type(experiment_data), intent(in) :: experiment(nexp)

  ! If the link was not found, it is consistent with experiments
  ! that did not observed it

  if ( .not. link%found ) then
    do iexp = 1, nexp
      if ( .not. link%exp(iexp)%observed ) then
        do j = 1, experiment(iexp)%ntypes
          if ( .not. ( link .matches. experiment(iexp)%linktype(j) ) ) cycle
          if ( link%exp(iexp)%type_reactive ) then
            link%exp(iexp)%type_consistent = .true.
            link%n_type_consistent = link%n_type_consistent + 1
            if ( link%exp(iexp)%obs_reactive ) then
              link%exp(iexp)%obs_consistent = .true.
              link%n_obs_consistent = link%n_obs_consistent + 1
            end if
          end if
        end do
      end if
    end do
  end if

  ! If the link was found, it is consistent only with experiments
  ! that reported having observed it

  if ( link%found ) then
    do iexp = 1, nexp
      do j = 1, experiment(iexp)%ntypes
        if ( .not. ( link .matches. experiment(iexp)%linktype(j) ) ) cycle

        ! If the distance is good and the link was observed, the structure is consistent

        if ( link%topodist <= experiment(iexp)%linktype(j)%dist .and. &
             link%exp(iexp)%observed ) then
          if ( link%exp(iexp)%type_reactive ) then
            link%exp(iexp)%type_consistent = .true.
            link%n_type_consistent = link%n_type_consistent + 1
            if ( link%exp(iexp)%obs_reactive ) then
              link%exp(iexp)%obs_consistent = .true.
              link%n_obs_consistent = link%n_obs_consistent + 1
            end if
          end if
        end if

        ! If the link was not observed, and the distance is too large, the structure
        ! is also consistent

        if ( link%topodist > experiment(iexp)%linktype(j)%dist .and. &
             .not. link%exp(iexp)%observed ) then
          if ( link%exp(iexp)%type_reactive ) then
            link%exp(iexp)%type_consistent = .true.
            link%n_type_consistent = link%n_type_consistent + 1
            if ( link%exp(iexp)%obs_reactive ) then
              link%exp(iexp)%obs_consistent = .true.
              link%n_obs_consistent = link%n_obs_consistent + 1
            end if
          end if
        end if

      end do
    end do
  end if

end subroutine linkconsistency
