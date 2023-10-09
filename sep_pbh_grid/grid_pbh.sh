#!/bin/bash

export OMP_NUM_THREADS=8

mkdir grid

run () {
    M=$1
    MBHf=$2
    
    MBH=$(echo "10^$MBHf * $M" | bc -l)
    
    echo $M, $MBH
    dirname="grid/M=${M}_MBHf=${MBHf}"
    cp -R template $dirname
    cd $dirname
    shmesa change inlist_project \
        initial_mass $M \
        new_core_mass $MBH
    ./star inlist_project
    
    rm -rf .mesa_temp_cache photos star 
    
    cd -
}

for M in 1; do
    for MBHf in -12 -11 -10 -4; do
        run $M $MBHf
    done
done
