
subroutine writelog(string)

  use inputoptions

  character(len=*) :: string

  write(*,*) " "//trim(adjustl(string))

  if ( output_log ) then
    write(output_log_unit,*) " "//trim(adjustl(string))
  end if

  return

end subroutine writelog
