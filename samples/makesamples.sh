#!/bin/bash

EXE=felix.INT64Nifort.d

for dir in GaAs SrTiO3
do 
    echo $dir
    pushd $dir
    pwd
    ../../src/$EXE
    rm -rf sample_outputs/*
    mv $dir_Sim_*/* sample_outputs
    rmdir $dir_Sim_*
    popd
done


