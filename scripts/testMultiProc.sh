#!/bin/bash

nProc=(1 2 3 4 5 6 7 8 9 10 50 51)

for nP in ${nProc[@]}; do
    output_directory="MyImplementation_${nP}"
    mpiexec -n "$nP" ./mainrespa -a eauta -i 5 -r 0 -d 0.002 -gx 0 -gy 0 -gz 0 -s 1.0 -e 1.0 -nu 0.073 -csv tools/uniform.csv -o "$output_directory"/out.csv
done