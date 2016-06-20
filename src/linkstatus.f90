!
! Function that sets the value of the status property of the computed links,
! according to the result of the computations and the observations on that link
!
! Link status:
!
!                                                          It is consistent with   --YES: Status: 0 
!                                                            dmin and dmax         |      (FOUND GOOD)
!                                                    --YES-------------------------|
!                                 The link was found |                             --NO: Status: 1 
!                                 (on structure)     |                                   (FOUND VIOL)
!                           --YES--------------------|                                
!                           |                        |                                --YES: Status: 2
!                           |                        |    Euclidean violation of dcut |      (EUCL VIOL)
!  The link was observed -- |                        --NO-----------------------------|
!   (experimentally)        |                                                         |
!                           |                                                         --NO: Status: 3
!                           |                                                              (NOTFOUND VIOL)
!                           |
!                           |                                                                         
!                           |                              The distance is         --YES: Status: 5   
!                           |                              shorter than dmax       |      (NOTFOUND GOOD)
!                           |                       --YES:-------------------------|                   
!                           |                       |                              --NO: Status: 1     
!                           |    The link was found |                                    (FOUND VIOL)                  
!                           |    (on structure)     |                   
!                           --NO--------------------|                                                                   
!                                                   |                               --YES: Status: 5    
!                                                   |   Euclidean violation of dcut |      (EUCL GOOD)
!                                                   --NO----------------------------|                 
!                                                                                   |                   
!                                                                                   --NO: Status: 6     
!                                                                                        (NOTFOUND GOOD)
!

integer function linkstatus(link)

  use topolink_data
  implicit none
  type(specific_link), intent(in) :: link

  if ( link%observed ) then
    if ( link%found ) then
      if ( link%topodist >= link%dmin .and. &
           link%topodist <= link%dmax ) then
        linkstatus = 0 ! FOUND GOOD
      else
        linkstatus = 1 ! FOUND VIOL
      end if
    else
      if ( link%euclidean > link%dcut ) then
        linkstatus = 2 ! EUCL VIOL
      else
        linkstatus = 3 ! NOTFOUND VIOL
      end if
    end if
  else
    if ( link%found )then
      if ( link%topodist <= link%dmax ) then
        linkstatus = 5 ! NOTFOUND GOOD
      else
        linkstatus = 1 ! FOUND VIOL
      end if
    else
      if ( link%euclidean > link%dcut ) then
        linkstatus = 5 ! EUCL GOOD
      else
        linkstatus = 6 ! NOTFOUND GOOD
      end if
    end if
  end if

end function linkstatus
