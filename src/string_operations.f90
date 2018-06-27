!
! Module with functions to operate on file names and strings
!

module string_operations

  use ioformat, only : max_string_length

  contains

    !
    ! Function that determines the basename of a file,
    ! removing the path and the extension
    !
    
    character(len=max_string_length) function basename(filename)
    
      use ioformat, only : max_string_length
      integer :: i
      character(len=max_string_length) :: filename
    
      basename = trim(adjustl(filename))
      i = length(basename)
      idot = i+1
      do while(basename(i:i) /= "/")
        if ( basename(i:i) == "." ) then
          idot = i
        end if
        i = i - 1
        if ( i == 0 ) exit
      end do
      i = i + 1
      basename = basename(i:idot-1)
      do i = idot, max_string_length
        basename(i:i) = achar(32)
      end do
    
    end function basename
    
    !
    ! Subroutine that removes the path and trailing blanks from a file name
    !
    
    subroutine cleanname(string)
    
      use ioformat, only : max_string_length
      integer :: i, j
      character(len=max_string_length) :: string
    
      ! find last '/' character
    
      j = 1
      do i = 1, length(string)
        if ( string(i:i) == '/' ) j = i
      end do
      string = string(j:length(string))
      string = trim(adjustl(string))
    
    end subroutine cleanname

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
        if ( length == 0 ) exit
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
      character(len=max_string_length) :: file
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

    !
    ! Gets keyword from input file
    !
    
    character(len=max_string_length) function keyword(string)
    
      use ioformat, only : max_string_length
      implicit none
      integer :: if, il
      character(len=max_string_length) :: string
    
      if = 1
      do while( empty_char(string(if:if)) .and. if < max_string_length)
        if = if + 1
      end do
      il = if
      do while( .not. empty_char(string(il:il)) .and. il < max_string_length)
        il = il + 1
      end do
      il = il - 1
      keyword = string(if:il)
    
    end function keyword

    !
    ! Gets keyword value from input file
    !
    
    character(len=max_string_length) function keyvalue(string,ivalue)
    
      use ioformat, only : max_string_length
      implicit none
      integer :: if, il, ivalue, i
      character(len=max_string_length) :: string
    
      ! Jump keyword 
    
      if = 1
      do while( empty_char(string(if:if)) .and. if < max_string_length )
        if = if + 1
      end do
      il = if
      do while( .not. empty_char(string(il:il)) .and. il < max_string_length )
        il = il + 1
      end do
    
      ! The keyword ended, now reading values
    
      do i = 1, ivalue
        il = il - 1
        if = il + 1
        do while( empty_char(string(if:if)) .and. if < max_string_length )
          if = if + 1
        end do
        il = if
        do while( .not. empty_char(string(il:il)) .and. il < max_string_length )
          il = il + 1
        end do
      end do
    
      keyvalue = string(if:il)
      if(length(keyvalue) == 0) then
        write(*,*) ' ERROR: Some keyword without value: '
        write(*,*) string(1:length(string))
        stop
      end if
    
    end function keyvalue

    !
    ! Gets file name value from keyword value
    !
    
    character(len=max_string_length) function filename(string)
    
      use ioformat, only : max_string_length
      implicit none
      integer :: if, il, i
      character(len=max_string_length) :: string
    
      ! Jump keyword 
    
      if = 1
      do while( empty_char(string(if:if)) .and. if < max_string_length )
        if = if + 1
      end do
      il = if
      do while( .not. empty_char(string(il:il)) .and. il < max_string_length )
        il = il + 1
      end do
    
      ! The keyword ended, now reading the file name (with spaces)
    
      il = il - 1
      if = il + 1
      do while( empty_char(string(if:if)) .and. if < max_string_length )
        if = if + 1
      end do

      ! If the first character is a quote, read until the next quote
      if ( string(if:if) == '"' ) then
        if = if + 1
        il = if + 1
        do while( string(il:il) /= '"' .and. il < max_string_length )
          il = il + 1
        end do
        if ( il == max_string_length ) then
          write(*,*) ' ERROR: file name defined starting with quote but end of name not found. '
          write(*,*) '        Names with paths of up to ',max_string_length,' characters are accepted. Too long? '
          write(*,*) '        If this length is a problem, change max_string_length in ioformat.f90 and recompile. '
          stop
        end if
        il = il - 1

      ! If the name does not start with a quote, read until next blank character

      else
        il = if
        do while( .not. empty_char(string(il:il)) .and. il < max_string_length  )
          il = il + 1
        end do
        il = il - 1
      end if

      ! Check if there was a comment on this line

      do i = if, il
        if ( string(i:i) == "#" ) then
          il = i - 1
          exit
        end if
      end do
    
      filename = string(if:il)
      if(length(filename) == 0) then
        write(*,*) ' ERROR: Could not read file name: '
        write(*,*) string(1:length(string))
        stop
      end if
    
    end function filename

    !
    ! Count the number of input data words in an input line
    !
    
    integer function countwords(record)
    
      use ioformat, only : max_string_length
      implicit none
      integer :: ioerr, i
      character(len=max_string_length) :: record
      character(len=1) :: string
    
      countwords = 0
      do
        countwords = countwords + 1
        read(record,*,iostat=ioerr) (string, i = 1, countwords)
        if ( ioerr /= 0 ) then
          countwords = countwords - 1
          exit
        end if
      end do
    
    end function countwords

    !
    ! Subroutine that removes tab characters from string
    !
    
    subroutine strclean(record)
    
      use ioformat, only : max_string_length
      character(len=max_string_length) :: record
      record = trim(adjustl(record))
      do i = 1, max_string_length
        if ( record(i:i) == achar(9) ) record(i:i) = achar(32)
      end do
    
    end subroutine strclean

end module string_operations









