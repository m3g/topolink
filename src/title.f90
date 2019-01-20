!
! Program TOPOLINK
!
! L. Martinez, 
! Institute of Chemistry, University of Campinas
! http://leandro.iqm.unicamp.br
! Out 5, 2015
!

! Subroutine that prints the title

subroutine title(iout)

  integer :: iout
  write(iout,*) 
  write(iout,"(t3,56('#'),' TOPOLINK ',55('#'))") 
  write(iout,"(t18,a)") 'Institute of Chemistry - University of Campinas. http://m3g.iqm.unicamp.br/topolink'
  write(iout,*) 
  write(iout,"(t105,a)") " Version 19.020 "
  write(iout,*) 

end subroutine title
