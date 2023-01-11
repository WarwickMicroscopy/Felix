#!/bin/bash

#EXE=felix.INT64Nifort.d
EXE=felix.INT64NGNU.d

for dir in GaAs SrTiO3
do 
    echo $dir
    pushd $dir
    pwd
    ../../src/$EXE
    diff -s -W132 $dir_Sim_*/ 
    popd
done


