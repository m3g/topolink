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
! TopoLink: A package to compute the likelihood of structural models
! based on surface accessible topological distances
! 2015
!

! Subroutine that prints the title

subroutine title()

  write(*,*) 
  write(*,"(t3,56('#'),' TOPOLINK ',55('#'))") 
  write(*,"(t15,a)") 'L. Martinez, Institute of Chemistry - University of Campinas. http://leandro.iqm.unicamp.br'
  write(*,*) 
  write(*,"(t105,a)") " Version 18.229 "
  write(*,*) 

end subroutine title
