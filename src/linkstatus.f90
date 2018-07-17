!
! Function that sets the value of the status property of the computed links,
! according to the result of the computations and the observations on that link
!
! Link status:
!
!                                                          It is consistent with   --YES: Status: 0 
!                                                            dmin and dmax         |      (OK: FOUND)
!                                                    --YES-------------------------|    
!                                                    |                             |      --Too short: Status: 1  
!                                                    |                             ---NO--|            (BAD: SHORT)
!                                                    |                                    |           
!                                                    |                                    --Too long: Status: 2
!                                 The link was found |                                                (BAD: LONG)
!                                  (on structure)    |                                                      
!                           --YES--------------------|                                
!                           |                        |        Euclidean violation     --YES: Status: 3
!                           |                        |             of dmaxlink        |     (BAD: EUCL)
!  The link was observed -- |                        --NO-----------------------------|
!   (experimentally)        |                                                         |
!                           |                                                         --NO: Status: 4 
!                           |                                                              (BAD: NOTFOUND)
!                           |
!                           |                                                                         
!                           |                              The distance is         --YES: Status: 5   
!                           |                              shorter than dmax       |      (BAD: MISSING)
!                           |                       --YES:-------------------------|                   
!                           |                       |                              --NO: Status: 6     
!                           |    The link was found |                                    (OK: LONG)                  
!                           |    (on structure)     |                   
!                           --NO--------------------|                                                                   
!                                                   |                               --YES: Status: 7    
!                                                   |   Euclidean violation of dcut |      (OK: EUCL)
!                                                   --NO----------------------------|                 
!                                                                                   |                   
!                                                                                   --NO: Status: 8     
!                                                                                        (OK: NOTFOUND)
!

integer function linkstatus(link)

  use topolink_data
  implicit none
  type(specific_link), intent(in) :: link

  if ( link%observed ) then
    if ( link%found ) then
      if ( link%topodist >= link%dmin .and. &
           link%topodist <= link%dmax ) then
        linkstatus = 0 ! OK: FOUND
      else
        if ( link%topodist < link%dmin ) then
          linkstatus = 1 ! BAD: SHORT
        end if
        if ( link%topodist > link%dmax .and. &
             link%euclidean <= link%dsearch ) then
          linkstatus = 2 ! BAD: LONG
        end if
        if ( link%euclidean > link%dsearch ) then
          linkstatus = 3 ! BAD: EUCL
        end if
      end if
    else
      if ( link%euclidean > link%dmaxlink ) then
        linkstatus = 3 ! BAD: EUCL
      else
        linkstatus = 4 ! BAD: NOTFOUND
      end if
    end if
  else
    if ( link%found )then
      if ( link%topodist <= link%dmaxlink ) then
        linkstatus = 5 ! BAD: MISSING
      else
        linkstatus = 6 ! OK: LONG
      end if
    else
      if ( link%euclidean > link%dmaxlink ) then
        linkstatus = 7 ! OK: EUCL
      else
        linkstatus = 8 ! OK: NOTFOUND
      end if
    end if
  end if

end function linkstatus




