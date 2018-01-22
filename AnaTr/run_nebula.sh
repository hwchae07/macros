#!/bin/bash


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
for ((i=26 ; i<=44 ; i++));
do
./AnaTree ./sdaq14/nebula00$i.ridf <<EOF
7
EOF
done

