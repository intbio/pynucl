exec mkdir -p tmp/dat

set render 0
set render [lindex $argv 0]

set timesh 0
set timesh [lindex $argv 1]

set ssflag 0
set ssflag [lindex $argv 2]

set step 1
set step [lindex $argv 3]

mol delete 0    
    
set first 0
set last 1000
set nframes [expr  [molinfo 1 get numframes] - 1 ]
#set nframes 20


#mol top 0
add_text_layer TIME

#animate goto 0
#render snapshot ../../analysis_data/initial.tga
set filen 1
for {set i $first} {$i <= $nframes} {incr i $step} {
animate goto $i

if {$ssflag == 1 } {
mol ssrecalc 1
}



mol delete top
add_text_layer TIME
draw color 0
set time [format "Time: %5.1f ns" [expr $i * 1.0]]
if {$timesh == 1 } {
draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3
}

display update

if {$render == 1 } {
render TachyonInternal tmp/dat/$filen.dat.tga
} else {

render snapshot tmp/dat/$filen.dat.tga
}
set filen [expr $filen + 1]
}
