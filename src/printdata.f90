!
! Subroutine that prints the output data of each link search
!

subroutine printdata(print,link)

  use ioformat
  use topolink_data
  implicit none
  integer :: print, ib
  type(specific_link) :: link
  character(len=1) :: ra1, ra2, aa1, aa2
  character(len=3) :: charobs
  character(len=9) :: charmax, chardist
  character(len=13) :: statuschar
  character(len=800) :: lineformat

  ! Print title

  if ( print == -1 ) then
    write(str,dashes) ; call writelog(str)
    write(str,"( '>     RESIDUE1   ATOM1 RESIDUE2   ATOM2 EUCLDIST  TOPODIST OBSERVED&
                 &    DMIN      DMAX        RESULT  OBSRES REACRES RA AA')") ; call writelog(str)
    write(str,dashes) ; call writelog(str)
    return
  end if

  ! Output line format for each link
  !  -------------------------------------------------------------------------------------------------------------------------
  !        RESIDUE1   ATOM1 RESIDUE2   ATOM2 EUCLDIST  TOPODIST OBSERVED    DMIN      DMAX        RESULT  OBSRES REACRES RA AA
  !  LINK: LYSX A 1000 XXCB LYSX A 1008 XXCB 0014.000 >0013.756   YES   0000.000 >0034.000 NOTFOUND GOOD   00/00   11/11 YN YN

  lineformat = "( t3,"//&
                 &"'LINK:',"//&           ! LINK:
                 &"t9,a4,"//&             ! LYS
                 &"t14,a1,"//&            ! A
                 &"t16,i4,"//&            ! 6
                 &"t21,a4,"//&            ! CB
                 &"t26,a4,"//&            ! LYS
                 &"t31,a1,"//&            ! A
                 &"t33,i4,"//&            ! 8
                 &"t38,a4,"//&            ! CB
                 &"t43,f8.3,"//&          ! 14.000 
                 &"t52,a9,"//&            ! 13.756
                 &"t64,a3,"//&            ! YES
                 &"t70,f8.3,"//&          ! 0.000
                 &"t79,a9,"//&            ! 34.000
                 &"t89,a13,"//&           ! FOUND GOOD
                 &"t105,i2,'/',i2,"//&    ! 00/00
                 &"t113,i2,'/',i2,"//&    ! 11/11
                 &"t119,a1,a1,"//&        ! YN
                 &"t122,a1,a1,"//&        ! YN
                 &")" 

  if ( link%observed ) then
    charobs = 'YES'
    write(charmax,"( f9.3 )") link%dmax
    ib = 9 - len(trim(adjustl(charmax)))
    charmax(ib:ib) = " "
  else
    charobs = 'NO'
    write(charmax,"( f9.3 )") link%dmin
    ib = 9 - len(trim(adjustl(charmax)))
    charmax(ib:ib) = ">"
  end if

  ! If a topological distance was found, print it

  if ( link%status == 0 .or. &
       link%status == 1 .or. &
       link%status == 2 .or. &
       link%status == 5 .or. &
       link%status == 6 ) then
    write(chardist,"( f9.3 )") link%topodist
    ib = 9 - len(trim(adjustl(chardist)))
    chardist(ib:ib) = " "

  else

  ! Otherwise, report that it is greater than dsearch, or than 
  ! the euclidean distance, the one which is greater

    write(chardist,"( f9.3 )") max(link%dsearch,link%euclidean)
    ib = 9 - len(trim(adjustl(chardist)))
    chardist(ib:ib) = ">"
  end if

  ra1 = "Y" ; if ( .not. link%atom1%residue%accessible ) ra1 = "N"
  ra2 = "Y" ; if ( .not. link%atom2%residue%accessible ) ra2 = "N"
  aa1 = "Y" ; if ( .not. link%atom1%accessible ) aa1 = "N"
  aa2 = "Y" ; if ( .not. link%atom2%accessible ) aa2 = "N"

  write(str,lineformat) link%atom1%residue%name, link%atom1%residue%chain, &
                        link%atom1%residue%index, link%atom1%name, &
                        link%atom2%residue%name, link%atom2%residue%chain, &
                        link%atom2%residue%index, link%atom2%name, &
                        link%euclidean, chardist, charobs, &
                        link%dmin, charmax, statuschar(link%status), &
                        link%n_obs_consistent, link%n_obs_expected, &
                        link%n_type_consistent, link%n_type_expected, &
                        ra1, ra2, aa1, aa2
  call writelog(str)

end subroutine printdata

function statuschar(status)

  implicit none
  integer :: status
  character(len=13) :: statuschar

  if ( status == 0 ) statuschar ="    OK: FOUND"
  if ( status == 1 ) statuschar ="   BAD: SHORT"
  if ( status == 2 ) statuschar ="    BAD: LONG"
  if ( status == 3 ) statuschar ="    BAD: EUCL"
  if ( status == 4 ) statuschar ="BAD: NOTFOUND"
  if ( status == 5 ) statuschar =" BAD: MISSING"
  if ( status == 6 ) statuschar ="     OK: LONG"
  if ( status == 7 ) statuschar ="     OK: EUCL"
  if ( status == 8 ) statuschar =" OK: NOTFOUND"

end function statuschar








