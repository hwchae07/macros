#!/bin/bash



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

make clean
make -j8

#for i in 26 27 28 36 37 38
#for i in 36 37 38
#do
#./AnaTree ./sdaq14/nebula00$i.ridf <<EOF
#7
#EOF
#done

#for i in 26 27 28 36 37 38
for i in 36 37 38
do
./AnaTree ./sdaq14/nebula00$i.ridf <<EOF
11
EOF
done


#for ((i=275;i<=294;i++)); do
#for ((i=275;i<=301;i++)); do
#for ((i=295;i<=301;i++)); do
#for ((i=240;i<=243;i++)); do
#for ((i=302;i<=302;i++));do
#for j in 1 2 3 4 5 6 7 8 9 10
#for j in 7
#do
#./AnaTree ./sdaq02/run0$i.ridf <<EOF
#$j
#EOF
#done
#./AnaTree $i <<EOF
#16
#EOF
#done
