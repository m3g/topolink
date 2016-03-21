!
! L. Martinez
! Institute of Chemistry
! University of Campinas
! http://leandro.iqm.unicamp.br
!

!
! Gets keyword from input file
!

function keyword(string)

  implicit none
  integer :: if, il
  character(len=200) :: keyword, string

  if = 1
  do while(string(if:if) <= ' '.and. if < 200)
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do
  il = il - 1
  keyword = string(if:il)

return
end function keyword

!
! Gets keyword value from input file
!

function keyvalue(string,ivalue)

  implicit none
  integer :: if, il, length, ivalue, i
  character(len=200) :: keyvalue, string

  ! Jump keyword 

  if = 1
  do while(string(if:if) <= ' '.and.if < 200)
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do

  ! The keyword ended, now reading values

  do i = 1, ivalue
    il = il - 1
    if = il + 1
    do while(string(if:if) <= ' '.and.if < 200)
      if = if + 1
    end do
    il = if
    do while(string(il:il) > ' '.and.il < 200)
      il = il + 1
    end do
  end do

  keyvalue = string(if:il)
  if(length(keyvalue) == 0) then
    write(*,*) ' ERROR: Some keyword without value: '
    write(*,*) string(1:length(string))
    stop
  end if

return
end function keyvalue

! Count the number of input data words in an input line

function countwords(record)

  implicit none
  integer :: countwords, ioerr, i
  character(len=200) :: record
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
! Function that sets the length of a string
!

function length(string)

  implicit none
  integer :: length
  character(len=200) :: string

  length = 200
  do while(string(length:length) <= ' ')
    length = length - 1
    if ( length == 0 ) exit
  end do

return
end function length

!
! Subroutine that removes the path and trailing blanks from a file name
!

subroutine cleanname(string)

  integer :: i, j, length
  character(len=200) :: string

  ! find last '/' character

  j = 1
  do i = 1, length(string)
    if ( string(i:i) == '/' ) j = i
  end do
  string = string(j:length(string))
  string = trim(adjustl(string))

end subroutine cleanname









