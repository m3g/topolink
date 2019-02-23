
function linkname( link :: LinkData )

  atom1 = "$(link.atom1.residue.name) $(link.atom1.residue.chain) $(link.atom1.residue.number) $(link.atom1.name)" 
  atom2 = "$(link.atom2.residue.name) $(link.atom2.residue.chain) $(link.atom2.residue.number) $(link.atom2.name)" 
  return "$atom1 - $atom2"

end
