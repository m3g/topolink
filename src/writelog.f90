
subroutine writelog(string)

  use inputoptions

  character(len=*) :: string

  if ( screen_log ) then
    write(*,*) trim(adjustl(string))
  end if

  if ( output_log ) then
    write(output_log_unit,*) trim(adjustl(string))
  end if

  return

end subroutine writelog
