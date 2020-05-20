#to make movie run vmd -e view_nucl.tcl -args big_data/h3-h4_xray.pdb big_data/h3-h4.xtc title 1 1 1 1 1 1 1 180 0
#first number - smoothing window \
#second number - 0/1 do movie of preview
#third number - 0/1 render with tachyon or not (tachyon allows commandline rendering)
#forth number -0/1 display time or not
#fifth number -0/1 update or not ssecondary structure during movie
#sixth number movie step in frames
#senventh number - timestep
#eighth numer - rotation angle around y
#ninth numer - rotation angle around x
set mov [lindex $argv 4] 

set title [lindex $argv 2] 

set sm 1 
set sm [lindex $argv 3] 

set init big_data/h3-h4_xray.pdb
set init [lindex $argv 0]

set trj big_data/h3-h4.xtc
set trj [lindex $argv 1]

set render 0
set render [lindex $argv 5]

set timesh 0
set timesh [lindex $argv 6]

set ssflag 0
set ssflag [lindex $argv 7]

set movstep 1
set movstep [lindex $argv 8]

set timestep 1
set timestep [lindex $argv 9]

set rotby 180
set rotby [lindex $argv 10]

set rotbyx 0
set rotbyx [lindex $argv 11]

source VMD_scripts/input_param.tcl
source VMD_scripts/add_text_layer.tcl

#display rendermode GLSL

#set title "MD simulations of H3-H4, tails truncated"

####
#Display the reference crystal structure
mol load pdb $init

mol modstyle 0 0 NewCartoon 0.840000 20.000000 2.630000 0

###

mol load pdb $init
mol addfile $trj waitfor all

mol ssrecalc top

set molid [ molinfo top ]
 

mol delrep 0 top

#  H3 H3 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Chain
mol selection {chain A E}
mol material AOShiny
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0
mol smoothrep top 0 $sm

# H4 H4 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Chain
mol selection {chain B F}
mol material AOShiny
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0
mol smoothrep top 1 $sm



# H2a H2a Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Chain
mol selection {chain C G}
mol material AOShiny
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0
mol smoothrep top 2 $sm



# H2b H2b Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Chain
mol selection {chain D H}
mol material AOShiny
mol addrep top
mol selupdate 3 top 0
mol colupdate 3 top 0
mol smoothrep top 3 $sm



# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol color ColorID 6
mol selection {nucleic and backbone}
mol material AOEdgy
mol addrep top
mol selupdate 4 top 0
mol colupdate 4 top 0
mol smoothrep top 4 $sm



# Dna nucleobases representation
mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 6
mol selection { nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
mol material AOEdgy
mol addrep top
mol selupdate 5 top 0
mol colupdate 5 top 0
mol smoothrep top 5 $sm



mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 3
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
mol material AOEdgy
mol addrep top
mol selupdate 6 top 0
mol colupdate 6 top 0
mol smoothrep top 6 $sm


#Key argininges
mol representation VDW
mol color ColorID 3
mol selection {(chain A E and resid 83 63 49) or (chain B F and resid 45) or (chain C G and resid 42 77) or (chain D H and resid 30)}
mol material AOShiny
mol addrep top
mol selupdate 7 top 0
mol colupdate 7 top 0
mol smoothrep top 7 $sm

# 5' end of DNA chain I

set ai [atomselect top "chain I and name O5'"]
set na [lindex [$ai get resid] 0]
set selchI [format "chain I and name O5' and resid \"%s\"" $na]

mol representation VDW 3 12
mol color ColorID 11
mol selection $selchI
mol material AOShiny
mol addrep top
mol selupdate 8 top 0
mol colupdate 8 top 0
mol smoothrep top 8 $sm

#H3F104C-H4V43C sCx1
#H3L82C-H4V81C sCx2

mol representation VDW 1 12
mol color ColorID 21
mol selection {(chain A E and resid 104) or (chain B F and resid 43)}
mol material AOShiny
mol addrep top
mol selupdate 9 top 0
mol colupdate 9 top 0
mol smoothrep top 9 $sm

#color scheme  for histones and DNA
#old scheme
# color chain A orange2
# color chain B yellow3
# color chain C mauve
# color chain D cyan3
# color chain E orange2
# color chain F yellow3
# color chain G mauve
# color chain H cyan3

color Chain A blue3
color Chain B green
color Chain C yellow2
color Chain D red3
color Chain E blue3
color Chain F green
color Chain G yellow2
color Chain H red3

color Highlight Nonback 6
color Highlight Nucback 2
color Display Background white
color scale method RWB
color change rgb  0 0.1 0.2 0.7 ;# blue
color change rgb  1 0.7 0.2 0.1 ;# red
color change rgb  3 0.7 0.4 0.0 ;# orange
color change rgb  4 0.8 0.7 0.1 ;# yellow
color change rgb  7 0.1 0.7 0.2 ;# green
color change rgb 10 0.1 0.7 0.8 ;# cyan
color change rgb 11 0.6 0.1 0.6 ;# purple
color change rgb 23 0.550000011920929 0.7099999785423279 0.9800000190734863
color change rgb 24 0.15 0.25 0.93
color change rgb 30 0.81 0.18 0.18

# ================= movie starts ==================== 
# sets a variable for knowing were real molecules ends if there are more than one
 set molnum [expr [molinfo num]]

# #--------set start txt position in px-----
 set txtx [expr -($dispw/2) + 20 ]
 set txty [expr ($disph/2) - 20 ]

 # set step of txt lines
 set txtstep  [expr $disph/30 ]

# # set counter for lines
 set txtlncount 0

# #prepare scene
scale by $scale
translate by $transx $transy $transz
#added to show dimer from interesting side
rotate y by $rotby
rotate x by $rotbyx
axes location off
display update ui

# #add some text on new layer
 add_text_layer HEADNAME
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " $title size 1.5 thickness 3
 incr txtlncount 2

 add_text_layer HISTONES

# #Turn off all
set num_of_reps [molinfo $molid get numreps]
#for {set i 0} { $i < $num_of_reps } { incr i 1 } {
#mol showrep $molid $i off }

set matnum 22

set text(0) "Histones H3"
set text(1) "Histones H4"
set text(2) "Histones H2A"
set text(3) "Histones H2B"
set text(4) "Min groove ARG"
set color_text(0) blue3
set color_text(1) green
set color_text(2) yellow2
set color_text(3) red3
set color_text(4) orange

# #First four pair of histones
for {set r 0} {$r < 5} {incr r 1} {
	#text
	draw color $color_text($r)
	draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "$text($r)" size 1.5 thickness 3
	incr txtlncount 1
	#molecule
	#mol showrep $molid $r on
	#copy a material for them
	#material add copy AOShiny
	# mol modmaterial $r top Material$matnum
	#mol modmaterial $r top AOShiny
	
	
# 	#make a smooth transparent appearance for 30 frames with rotation
	
	#	material change opacity Material$matnum 1
		#render Tachyon dat/$fc.dat
		
	

# 	

	#incr matnum 1
}

# # Adding DNA
# #text
add_text_layer DNA
draw color 6
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "DNA" size 1.5 thickness 3
incr txtlncount 1

draw color 3
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "AT pairs" size 1.5 thickness 3
incr txtlncount 1

draw color 6
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "GC pairs" size 1.5 thickness 3
incr txtlncount 1

draw color 11
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "O5'DNA chain I" size 1 thickness 3
incr txtlncount 1

draw color 21
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "sCX1 crosslinks" size 1 thickness 3
incr txtlncount 1
#ion lables
# draw color 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Sodium ions" size 1.5 thickness 3
# incr txtlncount 1

# draw color 7
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chloride ions" size 1.5 thickness 3
# incr txtlncount 1

#molecule
mol showrep $molid 4 on
mol showrep $molid 5 on
mol showrep $molid 6 on
mol showrep $molid 7 on
mol showrep $molid 8 on
mol modmaterial 4 $molid AOEdgy
mol modmaterial 5 $molid AOEdgy
mol modmaterial 6 $molid AOEdgy
mol modmaterial 7 $molid AOEdgy
mol modmaterial 8 $molid AOEdgy

	display update
	#material change opacity AOEdgy 0.8
	
display update


#mol top 1

proc src {file args} {
  set argv $::argv
  set argc $::argc
  set ::argv $args
  set ::argc [llength $args]
  set code [catch {uplevel [list source $file]} return]
  set ::argv $argv
  set ::argc $argc
  return -code $code $return
}


if {$mov == 1 } {
puts "Will make images for movie"

src VMD_scripts/movie.tcl $render $timesh $ssflag $movstep $timestep

exit
}
