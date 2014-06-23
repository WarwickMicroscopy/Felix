#!/bin/bash
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# FelixSim
#
# Richard Beanland, Keith Evans and Rudolf A Roemer
#
# (C) 2013/14, all right reserved
#
# Version: VERSION
# Date:    DATE
# Time:    TIME
# Rls:     RLSTATUS
# Build:   BUILD
# Author:  AUTHOR
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of FelixSim.
#
#  FelixSim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  FelixSim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo "build script for FelixSim/Draw/Refine"

version=${1:-0.0}
rls=${2:-0}
build=${3:-0.0}
author=${4:-0}
date=${5:-`date -R | cut -b 1-16`}
time=${6:-`date -R | cut -b 18-31`}

echo "attempting to build version" ${version} "with rls" ${rls} "build" ${build} "for author" ${author} "on" ${date} "at time" ${time} 

sourcedir=`pwd`

cd ..
targetdir=`pwd`/felix-${version}
[ -d ${targetdir} ] || mkdir ${targetdir}
cd ${sourcedir}

tarball=felix-${version}.tar.bz2

# tag the files with the right version/rls/build/author information

cd ${sourcedir}
cp -vr * ${targetdir}
cd ${targetdir}

echo "--- working on files in directory" $dir

for file in `find . \( -name "*.f90" -o -name "*.inp" -o -name "makefile*.GF" -o -name "README.txt" \) -print`; do

    echo $file " updating!"
    sed "s/:VERSION:/${version}/g" $file | sed "s/:DATE:/${date}/g" | sed "s/:TIME:/${time}/g" | sed "s/:RLSTATUS:/${rls}/g" | sed "s/:BUILD:/${build}/g" | sed "s/:AUTHOR:/${author}/g" > $file.tmp 

	mv $file.tmp $file
	#ls $file.tmp $file

done

# create the tarball for OBS deployment

pwd
cd ..
#cp -vr ${sourcedir}/* ${targetdir}

echo "--- creating tarball" ${tarball} "from files in" ${targetdir}
tar -cjf ${tarball} `basename ${targetdir}`
ls
