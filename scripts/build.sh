#!/bin/bash

echo "build script for FelixSim/Draw/Refine"

version=${1:-0.0}
status=${2:-0}
build=${3:-0.0}
author=${4:-0}
date=${5:-`date -R | cut -b 1-16`}
time=${6:-`date -R | cut -b 18-31`}

echo "attempting to build version" ${version} "with status" ${status} "build" ${build} "for author" ${author} "on" ${date} "at time" ${time} 

srcdirs="src samples docs"

for dir in src samples docs; do
    echo "--- working on files in directory" $dir

    cd $dir

    for file in *.f90 makefile*.GF *.txt *.inp; do

	echo $file " updating!"
	sed "s/VERSION/${version}/g" $file | sed "s/DATE/${date}/g" | sed "s/TIME/${time}/g" | sed "s/STATUS/${status}/g" | sed "s/BUILD/${build}/g" | sed "s/AUTHOR/${author}/g" > $file.tmp 

	#mv $file.tmp $file

    done

    cd ..

done
