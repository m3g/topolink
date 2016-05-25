!
! Program TOPOLINK
!
! L. Martinez, 
! Institute of Chemistry, University of Campinas
! http://leandro.iqm.unicamp.br
! Out 5, 2015
!
! Reference:
!
! L. Martinez, A. Ferrari, F. C. Gozzo,
! TopoLink: A package to compute the likelyhood of structural models
! based on surface accessible topological distances
! 2015
!

! Subroutine that prints the title

subroutine title()

  write(*,*) 
  write(*,"(a)") "  ########################################### TOPOLINK #################################################"
  write(*,"(a)") '        L. Martinez, Institute of Chemistry - University of Campinas. http://leandro.iqm.unicamp.br'
  write(*,*) 
  write(*,"( t90, a )") " Version 16.146 "
  write(*,*) 

end subroutine title

