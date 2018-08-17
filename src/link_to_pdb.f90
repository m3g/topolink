!
! Subroutine that prints the output data of each link search
!

subroutine link_to_pdb(link,nlinkatoms,linkdir,pdbfile,natoms,atom,n,x)

  use ioformat
  use topolink_data
  use string_operations
  implicit none

  integer :: ioerr
  integer :: j, ix, iy, iz
  integer :: natoms, nlinkatoms

  type(specific_link) :: link
  type(pdbatom) :: writeatom, atom(natoms)

  character(len=max_string_length) :: linkfile, pdbfile, linkdir
  character(len=4) :: char1, char2

  integer :: n
  double precision :: x(n)
  double precision :: overlap, stretch
  character(len=13) :: statuschar

  linkfile = pdbfile
  call cleanname(linkfile)
  write(char1,"( i4 )") link%atom1%residue%index
  write(char2,"( i4 )") link%atom2%residue%index
  linkfile=trim(adjustl(linkdir(1:length(linkdir))))//linkfile(1:length(linkfile)-4)//&
    '_'//trim(adjustl(link%atom1%residue%name))//trim(adjustl(link%atom1%residue%chain))//&
    trim(adjustl(char1))//trim(adjustl(link%atom1%name))//&
    '-'//trim(adjustl(link%atom2%residue%name))//trim(adjustl(link%atom2%residue%chain))//&
    trim(adjustl(char2))//trim(adjustl(link%atom2%name))//'.pdb'
  open(10,file=linkfile,iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not create link PDB file: ', trim(adjustl(linkfile))
    write(*,*) '        Perhaps the output directory does not exist.'
    write(*,*) '        Output directory: ', trim(adjustl(linkdir))
    stop
  end if
  write(10,"( 'REMARK EUCLIDEAN DISTANCE: ', (tr2,f12.5) )") link%euclidean
  write(10,"( 'REMARK TOPOLOGICAL DISTANCE: ', (tr2,f12.5) )") link%topodist
  write(10,"( 'REMARK OVERLAP and STRETCH: ', 2(tr2,f12.5) )") overlap(n,x), stretch(n,x)
  write(10,"( 'REMARK LINK STATUS: ', a13 )") statuschar(link%status)
  writeatom = atom(link%atom1%index)
  writeatom%residue%name = "LINK"
  writeatom%residue%chain = "A"
  writeatom%residue%index = 1
  writeatom%index = 1
  write(10,"(a)") trim(print_pdbhetatm(writeatom))
  do j = 1, nlinkatoms
    ix = (j-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    writeatom%index = j + 1
    writeatom%name = "O"
    writeatom%residue%name = "LINK"
    writeatom%residue%chain = "A"
    writeatom%residue%index = 1
    writeatom%x = x(ix)
    writeatom%y = x(iy)
    writeatom%z = x(iz)
    write(10,"(a)") trim(print_pdbhetatm(writeatom))
  end do
  writeatom = atom(link%atom2%index)
  writeatom%residue%name = "LINK"
  writeatom%residue%chain = "A"
  writeatom%residue%index = 1
  writeatom%index = nlinkatoms + 2
  write(10,"(a)") trim(print_pdbhetatm(writeatom))
  close(10)

end subroutine link_to_pdb

