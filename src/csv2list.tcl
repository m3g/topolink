#!/usr/bin/tclsh

set file [ open $argv ]
set filedata [ read $file ]
set filedata [ split $filedata "\n" ] 
close $file

set i 0
foreach line $filedata {
  incr i
  set dataline [ split $line "," ]
  set j 0
  foreach data $dataline {
    incr j
    set table($i,$j) $data
  }
}
set n [ expr $i - 1 ]

for { set i 1 } { $i <= $n } { incr i } {
  set temp [ string map {E GLU D ASP \
                         K LYS S SER \
                         C CYS \
                         ( "-" ) " "} $table($i,1) ]
  set table($i,1) $temp
  set temp [ string map {E GLU D ASP ( "-" ) " "} $table(1,$i) ]
  set table(1,$i) $temp
}

for { set i 2 } { $i <= [ expr $n - 1 ] } { incr i } {
  for { set j [ expr $i + 1 ] } { $j <= $n } { incr j } {
    set res1 [ split $table($i,1) "-" ]
    set res2 [ split $table(1,$j) "-" ]
    if { $table($i,$j) > " " } { 
      puts [ format "observed %3s A %3i %3s A %3i %3i" \
             [ lindex $res1 0 ] [ lindex $res1 1 ] \
             [ lindex $res2 0 ] [ lindex $res2 1 ] $table($i,$j) ]
    } elseif { $table($j,$i) > " " } {
      puts [ format "observed %3s A %3i %3s A %3i %3i" \
             [ lindex $res1 0 ] [ lindex $res1 1 ] \
             [ lindex $res2 0 ] [ lindex $res2 1 ] $table($j,$i) ]
    }
  }
}


