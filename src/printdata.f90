!
! Subroutine that prints the output data of each link search
!

subroutine printdata(print,link)

  use ioformat
  use topolink_data
  implicit none
  integer :: print, ib
  type(specific_link) :: link
  character(len=3) :: charobs
  character(len=9) :: charmax, chardist
  character(len=13) :: charresult
  character(len=800) :: lineformat

  ! Print title

  if ( print == -1 ) then
    write(*,dashes) 
    write(*,"( '        RESIDUE1   ATOM1 RESIDUE2   ATOM2 EUCLDIST  TOPODIST OBSERVED&
               &    DMIN      DMAX        RESULT  OBSRES REACRES')")
    write(*,dashes)
    return
  end if

  ! Output line format for each link
  !  -------------------------------------------------------------------------------------------------------------------
  !        RESIDUE1   ATOM1 RESIDUE2   ATOM2 EUCLDIST  TOPODIST OBSERVED    DMIN      DMAX        RESULT  OBSRES REACRES
  !  LINK: LYSX A 1000 XXCB LYSX A 1008 XXCB 0014.000 >0013.756   YES   0000.000 >0034.000 NOTFOUND GOOD   00/00   11/11

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

  if ( link%status == 0 .or. &
       link%status == 1 .or. &
       link%status == 4 ) then
    write(chardist,"( f9.3 )") link%topodist
    ib = 9 - len(trim(adjustl(chardist)))
    chardist(ib:ib) = " "
  else
    if ( link%observed ) then
      write(chardist,"( f9.3 )") link%dmax
      ib = 9 - len(trim(adjustl(chardist)))
      chardist(ib:ib) = ">"
    else
      write(chardist,"( f9.3 )") link%dmin
      ib = 9 - len(trim(adjustl(chardist)))
      chardist(ib:ib) = ">"
    end if
  end if

  if ( link%status == 0 ) charresult ="    OK: FOUND"
  if ( link%status == 1 ) charresult ="   BAD: SHORT"
  if ( link%status == 2 ) charresult ="    BAD: LONG"
  if ( link%status == 3 ) charresult ="    BAD: EUCL"
  if ( link%status == 4 ) charresult ="BAD: NOTFOUND"
  if ( link%status == 5 ) charresult =" BAD: MISSING"
  if ( link%status == 6 ) charresult ="     OK: LONG"
  if ( link%status == 7 ) charresult ="     OK: EUCL"
  if ( link%status == 8 ) charresult =" OK: NOTFOUND"

  write(*,lineformat) link%atom1%residue%name, link%atom1%residue%chain, &
                      link%atom1%residue%index, link%atom1%name, &
                      link%atom2%residue%name, link%atom2%residue%chain, &
                      link%atom2%residue%index, link%atom2%name, &
                      link%euclidean, chardist, charobs, &
                      link%dmin, charmax, charresult, &
                      link%n_obs_consistent, link%n_obs_expected, &
                      link%n_type_consistent, link%n_type_expected

end subroutine printdata



