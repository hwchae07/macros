#!/bin/bash

rm AnaBeamPID.C

ln -s AnaBeamPID_for_gamma.C AnaBeamPID.C

make clean
make -j8


for ((i=146;i<=154;i++));do
#for j in 1 2 3 4 5 6 7 8 9
#for j in 7 
#do
#./AnaTree ./sdaq02/run0$i.ridf <<EOF
#$j
#EOF
#done
for j in 16
do
./AnaTree $i<<EOF
$j
EOF
done
done

