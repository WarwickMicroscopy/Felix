#!/bin/bash

#EXE=felix.INT64Nifort.d
EXE=felix.OPT64NGNU.d

for dir in GaAs SrTiO3
do 
    echo $dir
    pushd $dir
    pwd
    ../../src/$EXE
    diff -s -W132 GaAs_Sim_0930A_080x080/ sample_outputs/
    popd
done


