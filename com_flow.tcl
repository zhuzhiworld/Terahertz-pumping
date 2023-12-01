###This is an analysis program for calculating water flow rate.###

mol addfile [glob *.dcd] type dcd first 5000 last 55000 step 1 waitfor all
mol addfile [glob *.psf]

set upperEnd  35
set lowerEnd -35
set dcdfreq  1000
set step     1.0
puts "Compute water net flux, upperEnd=$upperEnd, lowerEnd=$lowerEnd, dcdFreq=$dcdfreq, step=$step"
puts "Please Wait...."
set outf [open "ftraj.dat" w]
proc set_status {} {
  global wat statusList upperEnd lowerEnd
  set statusList {}
  foreach z [$wat get z] {
    if {$z < $lowerEnd} {
      lappend statusList -1
    } elseif {$z > $upperEnd} {
      lappend statusList 1
    } else {
      lappend statusList 0
    }
  }
}

set wat [atomselect top "name OH2"]
set numFrame [molinfo top get numframes]
set total 0

molinfo top set frame 0
set_status


for {set fr 1} {$fr < $numFrame} {incr fr} {
  molinfo top set frame $fr
  set oldList $statusList
  set_status
  foreach oldSt $oldList newSt $statusList {
    if {$oldSt!=$newSt && $oldSt+$newSt!=0} {
      incr total [expr $newSt - $oldSt]
    }
  }
  puts $outf "$fr\t$total"
}


if {$total > 0} {
  puts "The net flow is [expr $total/2.0] water molecules along +z"
} elseif {$total < 0} {
  puts "The net flow is [expr -$total/2.0] water molecules along -z"
} else {
  puts "The net flow is 0"
}
set flow_average [expr ($total/2.0)*1E6/($numFrame*$dcdfreq*$step)]
close $outf
puts "flow_ave=$flow_average"
puts "GOOD LUCY, BABY!"
exit