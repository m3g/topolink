!
! Module with functions to operate on file names and strings
!

module file_operations

  use ioformat, only : max_string_length 

  contains

    !
    ! Function that determines the basename of a file,
    ! removing the path and the extension
    !
    
    character(len=max_string_length) function basename(filename)
     
      use ioformat, only : max_string_length 
      implicit none
      character(len=max_string_length) :: filename
    
      basename = remove_path(filename)
      basename = remove_extension(basename)
    
    end function basename
    
    !
    ! Function that removes the extension of a file name
    !
    
    character(len=max_string_length) function remove_extension(filename)
     
      use ioformat, only : max_string_length 
      implicit none
      integer :: i, idot
      character(len=max_string_length) :: filename
    
      remove_extension = filename
      i = length(remove_extension)
      idot = i+1
      do while(i > 0)
        if ( remove_extension(i:i) == "." ) then
          idot = i
        end if
        i = i - 1
      end do
      i = i + 1
      remove_extension = remove_extension(1:idot-1)
      do i = idot, max_string_length
        remove_extension(i:i) = achar(32)
      end do
    
    end function remove_extension

    !
    ! Function that removes the path from a file name
    !
    
    character(len=max_string_length) function remove_path(filename)
    
      use ioformat, only : max_string_length 
      implicit none
      integer :: i, ilength
      character(len=max_string_length) :: filename
    
      remove_path = trim(adjustl(filename))
      ilength = length(remove_path)
      i = ilength
      do while(remove_path(i:i) /= "/")
        i = i - 1
        if ( i == 0 ) exit
      end do
      i = i + 1
      remove_path(1:ilength-i+1) = remove_path(i:ilength)
      do i = ilength-i+2, max_string_length
        remove_path(i:i) = achar(32)
      end do
    
    end function remove_path
    
    !
    ! Function that determines the length of a string
    !
    
    integer function length(string)
    
      use ioformat, only : max_string_length 
      implicit none
      character(len=max_string_length) :: string
      length = max_string_length
      do while( empty_char(string(length:length)) ) 
        length = length - 1
      end do
    
    end function length
    
    !
    ! Function that determines if a character is empty
    !
    
    logical function empty_char(char)
    
      implicit none
      character :: char
      empty_char = .false.
      if ( char == achar(9) .or. &
           char == achar(32) .or. &
           char == '' ) then
        empty_char = .true.
      end if 
    
    end function empty_char
    
    !
    ! Function that checks if a line is a comment line
    !
    
    logical function comment(string)
      
      use ioformat, only : max_string_length 
      implicit none
      integer :: i
      character(len=max_string_length) :: string
      i = 1
      do while( empty_char(string(i:i)) .and. i < max_string_length ) 
        i = i + 1
      end do
      comment = .false.
      if ( string(i:i) == "#" .or. i == max_string_length ) comment = .true.
    
    end function comment

    !
    ! Subroutine that tests if file exists, and tries to open it.
    ! Asks the user if he/she wants the file to be overwritten 
    !

    subroutine checkfile(file)
    
      use ioformat, only : max_string_length 
      implicit none
      integer :: ioerr
      character(len=max_string_length ) :: file
      character(len=1) :: char
    
      open(10,file=file,status='new',action='write',iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Trying to create file: ', trim(adjustl(file)),' but file already exists '
        write(*,"(a,$)") '  Overwrite it? (Y/N): '
        read(*,*) char
        if ( char == "Y" ) then
          open(10,file=file,action='write',iostat=ioerr)
          if ( ioerr /= 0 ) then
            write(*,*) ' Could not open file. Quitting. '
            stop
          end if
          close(10)
        else
          write(*,*) ' Quitting. '
          stop
        end if
      end if
      close(10)
    
    end subroutine checkfile

end module file_operations









