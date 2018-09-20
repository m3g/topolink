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

subroutine title(iout)

  integer :: iout
  write(iout,*) 
  write(iout,"(t3,56('#'),' TOPOLINK ',55('#'))") 
  write(iout,"(t15,a)") 'L. Martinez, Institute of Chemistry - University of Campinas. http://leandro.iqm.unicamp.br'
  write(iout,*) 
  write(iout,"(t105,a)") " Version 18.263 "
  write(iout,*) 

end subroutine title
