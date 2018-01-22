#!/bin/bash

rm AnaBeamPID.C

ln -s AnaBeamPID_original.C AnaBeamPID.C

make clean
make -j8



#for ((i=133;i<=138;i++)); do
#for j in 2 3 4 5
#for j in 4
#do
#./AnaTree ./sdaq02/run0$i.ridf <<EOF
#$j
#EOF
#done
#./AnaTree $i <<EOF
#12
#EOF
#done

#for ((i=146;i<=154;i++));do
#for j in 1 2 3 4 5 6 7 8 9
#for j in 7 8 
#do
#./AnaTree ./sdaq02/run0$i.ridf <<EOF
#$j
#EOF
#done
#./AnaTree $i<<EOF
#12
#EOF
#done


for ((i=275;i<=294;i++)); do
#for ((i=275;i<=301;i++)); do
#for ((i=295;i<=301;i++)); do
#for ((i=240;i<=243;i++)); do
#for ((i=302;i<=302;i++));do
#for j in 1 2 3 4 5 6 7 8 9 10
for j in 6
do
./AnaTree ./sdaq02/run0$i.ridf <<EOF
$j
EOF
done
./AnaTree $i <<EOF
15
EOF
done
