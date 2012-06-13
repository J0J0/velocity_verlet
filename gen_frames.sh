#!/bin/zsh
# vim: ft=zsh

# particle count n
n=$(print output/particle_*.dat(N) | wc -w)
[[ $n -ge 1 ]] || exit
# frame count m
m=0

fids=()
for ((i=1; i <= n; i++)); do
    tmpname=output/particle_$(printf '%04i' $i).dat
    exec {tmpid}<$tmpname
    fids+=$tmpid
done

while read -t -u $fids[1]
do
    print $REPLY 
    for ((i=2; i <= n; i++)); do
        read -t -u $fids[$i] || exit 2
        print $REPLY 
    done
    print
    print
    let m++
done > output/frames.dat

for ((i=1; i <= n; i++)) do 
    tmpid=$fids[$i]
    exec {tmpid}<&-
done

: > output/frames_info.txt
print "framecount=$m" >> output/frames_info.txt
