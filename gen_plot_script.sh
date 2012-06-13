#!/bin/zsh
# vim: ft=zsh

[[ -r output/frames.dat && -r output/frames_info.txt ]] || exit
# get frame count as $framecount
source output/frames_info.txt

outfile="plot_particle_movement.gpl"

cat << PREAMBLE > $outfile
#!/usr/bin/env gnuplot
set term gif size 800, 600 animate delay 2 loop 1
set output "particle_movement.gif"
set xrange [0:5]
set yrange [0:5]
unset key

PREAMBLE

for ((i=0; i < framecount; i+=2)); do
    print "plot 'output/frames.dat' index $i using 2:3:(1+\$0) with points pt 6 lc variable" 
done >> $outfile

