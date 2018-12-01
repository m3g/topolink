!
! Subroutine that prints the links as PDB files
!

subroutine link_to_pdb(link,nlinkatoms,linkdir,pdbfile,natoms,atom,n,x)

  use ioformat
  use topolink_data
  use string_operations
  implicit none

  integer :: ioerr
  integer :: j, ix, iy, iz, index
  integer :: natoms, nlinkatoms, nsteps

  type(specific_link) :: link
  type(pdbatom) :: writeatom, atom(natoms)

  character(len=max_string_length) :: linkfile, pdbfile, linkdir
  character(len=4) :: char1, char2
  character(len=1) :: color

  integer :: n
  double precision :: x(n), stepx, stepy, stepz
  double precision :: overlap, stretch
  character(len=13) :: statuschar

  ! Define the color of the linker (using standard atom colors)
  if ( link%status == 0 ) then
    color = "N" ! Blue for links that are completely ok
  else
    if ( link%found ) then
      color = "S" ! Yellow for links that are too long, too short, or not observed 
    else 
      color = "O" ! Red for links that were not found at all
    end if
  end if

  linkfile = pdbfile
  call cleanname(linkfile)
  write(char1,"( i4 )") link%atom1%residue%index
  write(char2,"( i4 )") link%atom2%residue%index
  linkfile=trim(adjustl(linkdir(1:length(linkdir))))//linkfile(1:length(linkfile)-4)//&
    '_'//trim(adjustl(link%atom1%residue%name))//trim(adjustl(link%atom1%residue%chain))//&
    trim(adjustl(char1))//trim(adjustl(link%atom1%name))//&
    '-'//trim(adjustl(link%atom2%residue%name))//trim(adjustl(link%atom2%residue%chain))//&
    trim(adjustl(char2))//trim(adjustl(link%atom2%name))//'.pdb'
  open(11,file=linkfile,iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(str,*) ' ERROR: Could not create link PDB file: ', trim(adjustl(linkfile)) ; call writelog(str)
    write(str,*) '        Perhaps the output directory does not exist.' ; call writelog(str)
    write(str,*) '        Output directory: ', trim(adjustl(linkdir)) ; call writelog(str)
    stop
  end if
  write(11,"( 'REMARK EUCLIDEAN DISTANCE: ', (tr2,f12.5) )") link%euclidean
  write(11,"( 'REMARK LINK STATUS: ', a13 )") statuschar(link%status)
  if ( link%found ) then
    write(11,"( 'REMARK TOPOLOGICAL DISTANCE: ', (tr2,f12.5) )") link%topodist
    write(11,"( 'REMARK OVERLAP and STRETCH: ', 2(tr2,f12.5) )") overlap(n,x), stretch(n,x)
  end if
  index = 1
  writeatom = atom(link%atom1%index)
  writeatom%residue%name = "LINK"
  writeatom%residue%chain = "A"
  writeatom%residue%index = 1
  writeatom%index = index
  write(11,"(a)") trim(print_pdbhetatm(writeatom))
  if ( link%found ) then
    do j = 1, nlinkatoms
      index = index + 1
      ix = (j-1)*3 + 1
      iy = ix + 1
      iz = ix + 2
      writeatom%index = index
      writeatom%name = color
      writeatom%residue%name = "LINK"
      writeatom%residue%chain = "A"
      writeatom%residue%index = 1
      writeatom%x = x(ix)
      writeatom%y = x(iy)
      writeatom%z = x(iz)
      write(11,"(a)") trim(print_pdbhetatm(writeatom))
    end do
  else
    nsteps = int(link%euclidean)-2
    stepx = (link%atom2%x - link%atom1%x)/nsteps
    stepy = (link%atom2%y - link%atom1%y)/nsteps
    stepz = (link%atom2%z - link%atom1%z)/nsteps
    do j = 1, nsteps
      index = index + 1
      writeatom%index = index
      writeatom%name = color
      writeatom%residue%name = "LINK"
      writeatom%residue%chain = "A"
      writeatom%residue%index = 1
      writeatom%x = link%atom1%x + j*stepx
      writeatom%y = link%atom1%y + j*stepy
      writeatom%z = link%atom1%z + j*stepz
      write(11,"(a)") trim(print_pdbhetatm(writeatom))
    end do
  end if
  index = index + 1
  writeatom = atom(link%atom2%index)
  writeatom%residue%name = "LINK"
  writeatom%residue%chain = "A"
  writeatom%residue%index = index
  writeatom%index = nlinkatoms + 2
  write(11,"(a)") trim(print_pdbhetatm(writeatom))
  close(11)

end subroutine link_to_pdb

